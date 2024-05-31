#!/usr/bin/env Rscript

library(hyprcoloc)
library(argparse)
library(data.table)
library(arrow)
library(stringr)
library(IGUtilityPackage)
library(ggplot2)
library(tidyr)
library(rtracklayer)

parser <- ArgumentParser(description = 'Run HyprColoc for every eQTL locus to detect colocalisation between eQTL locus and GWAS locus.')

parser$add_argument('--loci', metavar = 'file', type = 'character', 
help = 'File which contains variant list for each locus.')
parser$add_argument('--eqtl_folder', metavar = 'file', type = 'character', 
help = 'eQTLGen parquet folder format with per gene output files.')
parser$add_argument('--gwas_folder', metavar = 'file', type = 'character', 
help = 'GWAS parquet folder format with per phenotype GWAS sumstats files.')
parser$add_argument('--gtf', metavar = 'file', type = 'character',
                    help = "ENSEMBL .gtf file, needs to be hg38.")
parser$add_argument('--pheno_manifest', metavar = 'file', type = 'character',
                    help = "UKBB phenotype description file.")
parser$add_argument('--i2_thresh', type = 'numeric', default = 40,
                    help = 'Heterogeneity threshold. Defaults to <40%.')
parser$add_argument('--maxN_thresh', type = 'numeric', default = 0.8,
                    help = 'Per gene maximal sample size threshold. Defaults to 0.8 (SNPs with >=0.8*max(N))')
parser$add_argument('--minN_thresh', type = 'numeric', default = 0,
                    help = 'Minimal sample size threshold. Defaults to 0 (no filtering)')
parser$add_argument('--output', metavar = 'file', type = 'character', 
help = 'Output file.')
#parser$add_argument('--locusplot', type = 'logical', action = 'store_false')

args <- parser$parse_args()

#TEMP
# args <- list(loci = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenEqtlGwasColoc/work/82/7c0de0b3c8b918cd8f9ceb5038137b/ENSG00000140932.txt.gz", 
# eqtl_folder = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenEqtlGwasColoc/work/82/7c0de0b3c8b918cd8f9ceb5038137b/eqtls/", 
# gwas_folder = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenEqtlGwasColoc/work/82/7c0de0b3c8b918cd8f9ceb5038137b/neale_ukb_sumstats_parquet",
# gtf = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenEqtlGwasColoc/work/82/7c0de0b3c8b918cd8f9ceb5038137b/Homo_sapiens.GRCh38.108.gtf.gz", 
# i2_thresh = 40, maxN_thresh = 0.8, minN_thresh = 6000, 
# pheno_manifest = "")

# Functions
ParseInput <- function(eqtl_database, locus_input, gwas_database){

  message("Reading locus info...")
  loci2 <- locus_input
  gene <- unique(loci2$phenotype)
  message("Reading locus info...done!")

  # Parse SNP list
  message("Parsing SNP list...")
  snps <- unlist(strsplit(loci2$SNPs, "|", fixed = TRUE))
  message(paste(length(snps), "SNPs in the filter."))

  message("Parsing SNP list...done!")

  # Read in eQTL gene
  message("Reading eQTL parquet file...")

  eqtls <- eqtl_database %>% 
  filter(variant %in% snps) %>% 
  collect() %>%
  filter((i_squared < args$i2_thresh | is.na(i_squared)) & 
  sample_size >= args$maxN_thresh * max(sample_size) & 
  sample_size >= args$minN_thresh) %>% 
  as.data.table()

  message("Reading eQTL parquet file...done!")
  message(paste(nrow(eqtls), "variants in the data!"))
  
  # Make beta and se matrix
  message("Make beta and se matrices...")
  eqtl_beta <- data.table(variant = eqtls$variant, eqtls$beta)
  eqtl_se <- data.table(variant = eqtls$variant, eqtls$standard_error)

  if (nrow(eqtl_beta[is.na(V2)]) > 0){message("NAs in eQTL beta table!")}
  if (nrow(eqtl_se[is.na(V2)]) > 0){message("NAs in eQTL se table!")}

  rm(eqtls)
  gc()

  # Read in GWAS sumstats
  # TODO: these are hardcoded, make those explicit
  message("Reading GWAS parquet file...")
  gwas <- gwas_database %>% 
    filter(ID %in% snps & info > 0.4) %>% 
    select(gwas_id, ID, beta, se) %>% 
    collect()

  gwas <- gwas %>% 
    group_by(gwas_id) %>% 
    distinct(ID, .keep_all = TRUE) %>%
    ungroup() %>% 
    as.data.table()

  if (nrow(gwas[is.na(beta)]) > 0){message("NAs in GWAS beta table!")}
  if (nrow(gwas[is.na(se)]) > 0){message("NAs in GWAS se table!")}

  message("Reading GWAS parquet file...done!")
  # beta file
  gwas_beta <- gwas[, c(1:3), with = FALSE]
  gwas_beta <- pivot_wider(gwas_beta, names_from = "gwas_id", values_from = "beta")
  colnames(gwas_beta)[1] <- "variant"
  # se file
  gwas_se <- gwas[, c(1:2, 4), with = FALSE]
  gwas_se <- pivot_wider(gwas_se, names_from = "gwas_id", values_from = "se")
  colnames(gwas_se)[1] <- "variant"

  gwas_beta <- gwas_beta[!duplicated(gwas_beta$variant), ]
  gwas_se <- gwas_se[!duplicated(gwas_se$variant), ]

  if (any(is.na(gwas_beta))){ message("NAs in GWAS beta table after pivoting!") }
  if (any(is.na(gwas_se))){ message("NAs in GWAS se table after pivoting!") }

  if(nrow(gwas_beta) > 100){ #So that there are at least some variants reasonable to run.
  
  eqtl_beta <- merge(eqtl_beta, gwas_beta, by = "variant")
  eqtl_se <- merge(eqtl_se, gwas_se, by = "variant")
  message(paste(nrow(eqtl_beta), "variants in the combined matrix"))

  } else {
    message("Does not make sense to include phenotype: <100 variants after filters!")
    }
  gc()

  colnames(eqtl_beta)[1:2] <- c("SNP", gene)
  colnames(eqtl_se)[1:2] <- c("SNP", gene)

  snplist <- eqtl_beta$SNP
  
  eqtl_beta <- as.matrix(eqtl_beta[, -1, with = FALSE])
  rownames(eqtl_beta) <- snplist
  eqtl_se <- as.matrix(eqtl_se[, -1, with = FALSE])
  rownames(eqtl_se) <- snplist

  return(list(betas = eqtl_beta, standard_errors = eqtl_se))
  message("Make beta and se matrices...done!")
}


# iterate over loci
# Prepare inputs
temp_loci <- fread(args$loci)
gene <- unique(temp_loci$phenotype)

message(paste("Analysing", gene))

message("Reading GTF...")
hgnc <- readGFF(args$gtf)
hgnc <- unique(hgnc[, c(9, 11)])
message("Reading GTF...done!")

message("Reading eQTL database...")
eqtls_input <- arrow::open_dataset(list.files(paste0(args$eqtl_folder, "/phenotype=", gene), full.names = TRUE))
message("Reading eQTL database...done!")

message("Reading GWAS database...")
gwas_input <- arrow::open_dataset(args$gwas_folder)
message("Reading GWAS database...done!")

res_final <- data.table(
eQTL_gene = NA,
eQTL_gene_name = NA,
type = NA,
lead_SNP = NA,
lead_SNP_chr = NA,
lead_SNP_pos = NA,
lead_SNP_ea = NA,
lead_SNP_nea = NA,
lead_SNP_beta = NA,
lead_SNP_se = NA,
iteration = NA,
traits = NA,
posterior_prob = NA,
regional_prob = NA,
candidate_snp = NA,
posterior_explained_by_snp = NA,
dropped_trait = NA,
nr_snps_included = NA
)[-1]


for (locus in 1:nrow(temp_loci)){
  message(paste0("Analysing ", locus, "/", nrow(temp_loci), " locus..."))
  message("Parsing inputs...")

print(temp_loci[locus])

inputs <- ParseInput(eqtl_database = eqtls_input, locus_input = temp_loci[locus], gwas_database = gwas_input)
  message("Parsing inputs...done!")

if (locus == 1){
  message("Reading pheno manifest file...")

  pheno_manifest <- fread(args$pheno_manifest)
  pheno_manifest2 <- pheno_manifest
  pheno_manifest <- pheno_manifest[, c(1, 3), with = FALSE]
  pheno_manifest$num_type <- NA
  pheno_manifest <- as.data.frame(pheno_manifest)
  pheno_manifest[pheno_manifest$variable_type %in% c("continuous_raw", "continuous_irnt", "categorical", "ordinal"), ]$num_type <- 0
  pheno_manifest[pheno_manifest$variable_type %in% c("binary"), ]$num_type <- 1
  pheno_manifest <- pheno_manifest[pheno_manifest$phenotype %in% colnames(inputs$betas)[-1], -2]

  pheno_manifest <- pheno_manifest[match(colnames(inputs$betas)[-1], pheno_manifest$phenotype), ]$num_type
  message(paste(length(pheno_manifest[pheno_manifest == 0]), "continuous traits and", length(pheno_manifest[pheno_manifest == 1]), "binary traits."))
  message("Reading pheno manifest file...done!")
}

gene <- colnames(inputs$betas)[1]

betas <- inputs$betas
ses <- inputs$standard_errors

rm(inputs)
gc()

if (ncol(betas) > 1) {

############################
# Analyse eQTL-GWAS regions#
############################

trait_names <- colnames(betas)

# Remove rows where SE is 0
# TODO: check why such variants are in the results
any_row_contains_zero <- apply(ses, 1, function(row) !any(row == 0 | is.na(row)))

if(any(is.na(betas))){message("NAs in the beta table!")}
if(any(is.na(ses))){message("NAs in the se table!")}

message("Running hyprcoloc analysis...")

res <- hyprcoloc(betas[any_row_contains_zero, ], 
                 ses[any_row_contains_zero, ], 
                 trait.names = colnames(betas), 
                 snp.id = rownames(betas[any_row_contains_zero, ]),
                 binary.outcomes = pheno_manifest
                 )  

message("Running hyprcoloc analysis...done!")
res <- as.data.table(res$results)

res <- data.table(eQTL_gene = gene, 
                  eQTL_gene_name = hgnc[hgnc$gene_id == gene, ]$gene_name,
                  type = temp_loci$type[locus],
                  lead_SNP = temp_loci$SNP[locus],
                  lead_SNP_chr = temp_loci$chr[locus],
                  lead_SNP_pos = temp_loci$pos[locus],
                  lead_SNP_ea = temp_loci$ea[locus],
                  lead_SNP_nea = temp_loci$nea[locus],
                  lead_SNP_beta = temp_loci$beta[locus],
                  lead_SNP_se = temp_loci$se[locus],
                  as.data.table(res)
                  )

res <- res[res$traits != "None", ]

message("Splitting traits to separate rows...")

res[, traits := list(strsplit(traits, ",", fixed = TRUE))]
res$nr_snps_included <- nrow(betas)

res <- res %>% 
separate_rows(traits, sep = "\\| ") %>% 
mutate(traits = str_trim(traits))

res <- res[res$eQTL_gene != res$traits, ]

message("Splitting traits to separate rows...done!")

res_final <- rbind(res_final, res)
rm(res)
gc()
}
}

message("Adding phenotype annotations...")
res_final <- merge(res_final, pheno_manifest2[, c(1, 2), with = FALSE], by.x = "traits", by.y = "phenotype")
res_final <- res_final[, c(2:12, 1, ncol(res_final), 13:(ncol(res_final) - 1)), with = FALSE]
res_final <- res_final[order(type, lead_SNP_chr, lead_SNP_pos, traits)]
message("Adding phenotype annotations...done!")

message("Writing output...")
fwrite(res_final, paste0(gene, "_coloc.txt"), sep = "\t")
message("Writing output...done!")
