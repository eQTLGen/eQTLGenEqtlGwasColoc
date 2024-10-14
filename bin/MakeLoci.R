#!/usr/bin/env Rscript

library(data.table)
library(stringr)
library(arrow)
library(dplyr)
library(rtracklayer)
library(IGUtilityPackage)
library(argparse)

setDTthreads(1)
parser <- ArgumentParser(description = 'Find eQTL loci.')

parser$add_argument('--reference', metavar = 'file', type = 'character',
                    help = 'eQTLGen SNP reference file in parquet format.')
parser$add_argument('--sig_res', metavar = 'file', type = 'character',
                    help = 'Significant eQTL results in eQTLGen format. Can be gzipped.')
parser$add_argument('--eqtl_folder', metavar = 'file', type = 'character',
                    help = 'Parquet folder with all files.')
parser$add_argument('--gtf', metavar = 'file', type = 'character',
                    help = "ENSEMBL .gtf file, needs to be hg38.")
parser$add_argument('--lead_variant_win', type = 'numeric', default = 1000000,
                    help = 'Distance threshold for distance pruning.')
parser$add_argument('--cis_win', type = 'numeric', default = 1000000,
                    help = 'Distance threshold around cis eQTL lead variant. This window is used in coloc analysis.')
parser$add_argument('--trans_win', type = 'numeric', default = 5000000,
                    help = 'Distance threshold to declare variant trans.')
parser$add_argument('--p_thresh', type = 'numeric', default = 5e-8,
                    help = 'P-value threshold for significant effects.')
parser$add_argument('--i2_thresh', type = 'numeric', default = 100,
                    help = 'Heterogeneity threshold. Defaults to <=100%.')
parser$add_argument('--maxN_thresh', type = 'numeric', default = 0.8,
                    help = 'Per gene maximal sample size threshold. Defaults to 0.8 (SNPs with >=0.8*max(N))')
parser$add_argument('--minN_thresh', type = 'numeric', default = 0,
                    help = 'Minimal sample size threshold. Defaults to 0 (no filtering)')
parser$add_argument('--gene_filter', metavar = 'file', type = 'character',
                    help = "Filter with ENSEMBL IDs to include to the analysis.")
# parser$add_argument('--max_lead_distance', type = 'numeric', default = 250000,
#                     help = 'Maximum distance between primary cis and trans lead variants.')

args <- parser$parse_args()

# TEMP
# args <- list(reference = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2/projectfiles/data/derived_data/2023-01-28_MetaAnalysis/data/1000G-30x.parquet", 
# sig_res = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2/projectfiles/data/derived_data/2023-06-04_MetaAnalysis/data/Subset_P5e8_2023-06-04.txt.gz", 
# eqtl_folder = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenCisTransColoc/tests/input/sumstats/", 
# gtf = "/Users/urmovosa/Documents/projects/2019/eQTLGenPhase2_pipelines/eQTLGenCisTransColoc/tests/Homo_sapiens.GRCh38.108.gtf.gz", 
# lead_variant_win = 500000, cis_win = 500000, trans_win = 5000000, p_thresh = 5e-8, 
# i2_thresh = 40, maxN_thresh = 0.5, minN_thresh = 6000
# )

message("Reading in sig. results...")
sig <- fread(args$sig_res, key = "SNP")
sig <- sig[P < args$p_thresh & (i_squared <= args$i2_thresh | is.na(i_squared))]

message(nrow(sig))

# gene filter
gene_filter <- fread(args$gene_filter, header = TRUE)
colnames(gene_filter)[1] <- "gene"
if (nrow(gene_filter) > 0){

    sig <- sig[phenotype %in% gene_filter$gene]

}

message("Filter results to genes available in full files...")
eqtl_genes <- str_replace(list.files(args$eqtl), ".*phenotype=", "")
sig <- sig[phenotype %in% eqtl_genes]
message(paste(length(unique(sig$phenotype)), "genes in full results"))
if (length(unique(sig$phenotype)) < 1){stop("No genes in full results, terminating.")}
message("Filter results to genes available in full files...done!")

sig <- sig %>% 
    group_by(phenotype) %>% 
    filter(N >= args$maxN_thresh * max(N) & N >= args$minN_thresh) %>% 
    as.data.table()

message(paste(nrow(sig), "rows among significant results"))

message("Reading in sig. results...done!")

message("Reading in SNP list in the full files...")

gene <- list.files(args$eqtl_folder)[1]
snp_list  <- arrow::open_dataset(list.files(paste0(args$eqtl_folder, "/", gene), full.names = TRUE))

snp_list <- snp_list %>% select(variant) %>% collect() %>% as.data.table()
message("Reading in SNP list in the full files...done!")

message(nrow(snp_list))

message("Reading in reference...")
ref <- arrow::open_dataset(args$reference)

ref <- ref %>% 
    filter(ID %in% !!snp_list$variant) %>% 
    collect()

message(nrow(ref))

ref <- data.table(ref, key = "ID")
ref <- ref[, c(1, 5, 2, 3, 4), with = FALSE]
message("Reading in reference...done!")

sig <- merge(sig, ref[, c(1:3), with = FALSE], by.x = "SNP", by.y = "ID")
message(paste(nrow(sig), "rows among significant results, after merging with reference"))

message("Finding lead variants for each gene...")
LeadVariants <- sig %>% 
    group_by(phenotype) %>% 
    group_modify(~ IdentifyLeadSNPs(.x, 
    snp_id_col = "SNP", 
    snp_chr_col = "CHR", 
    snp_pos_col = "bp", 
    eff_all_col = "alt_all", 
    other_all_col = "ref_all",
    beta_col = "beta", 
    se_col = "se", 
    p_col = "P", 
    window = args$lead_variant_win))

message(paste(nrow(LeadVariants), "loci"))    
message("Finding lead variants for each gene...done!")

rm(sig)
gc()

message("Removing lead variants mapping to MHC region...")
message(paste(nrow(LeadVariants[LeadVariants$chr == 6 & LeadVariants$pos > 25000000 & LeadVariants$pos < 34000000, ]), "such variants removed."))
LeadVariants <- LeadVariants[!(LeadVariants$chr == 6 & LeadVariants$pos > 25000000 & LeadVariants$pos < 34000000), ]
message("Removing lead variants mapping to MHC region...done!")

# TODO: check why threre  are multiple cis variants for some genes. Probably due to smaller than 1Mb "clumping" window.
# Annotate cis/trans
message("Annotating lead variants to cis/trans...")
ensg <- readGFF(args$gtf)
ensg <- as.data.table(ensg)
ensg <- ensg[type == "gene"]
ensg$tss <- ensg$start
ensg[strand == "-"]$tss <- ensg[strand == "-"]$end
ensg <- ensg[, c(9, 1, 27), with = FALSE]

Lead2 <-  merge(LeadVariants, ensg, by.x = "phenotype", by.y = "gene_id")
Lead2$type <- "interim"
Lead2[Lead2$chr == Lead2$seqid | abs(Lead2$pos - Lead2$tss) < args$cis_win, ]$type <- "cis"
Lead2[Lead2$chr != Lead2$seqid | abs(Lead2$pos - Lead2$tss) > args$trans_win, ]$type <- "trans"
Lead2 <- Lead2[Lead2$type == "cis" | (Lead2$type == "trans" & Lead2$P < args$p_thresh), ]

Lead2 <- as.data.table(Lead2)

fwrite(Lead2, "eQtlLeadVariants.txt.gz", sep = "\t")

summary_per_gene <- Lead2 %>% 
group_by(phenotype) %>% 
summarise(Nr_eQTL_lead_variants = length(SNP),
Nr_cis_eQTL_lead_variants = length(SNP[type == "cis"]),
Nr_trans_eQTL_lead_variants = length(SNP[type == "trans"]))

summary_per_gene <- as.data.table(summary_per_gene[order(summary_per_gene$Nr_eQTL_lead_variants, decreasing = TRUE), ])


fwrite(summary_per_gene, "eQtlLeadSummary.txt.gz", sep = "\t")

message("Annotating lead variants cis/trans...done!")

ref <- data.table(SNP = ref$ID, chr = ref$CHR, start = ref$bp, 
end = ref$bp + 1)
ref$chr <- as.factor(ref$chr)
setkey(ref, chr, start, end)

if (length(Lead2$phenotype) <= 1000){
batches = 1
}else{
batches <- c(seq(from = 1, to = length(Lead2$phenotype), by = 1000), length(Lead2$phenotype))
}

res <- data.table(lead_SNP = NA, SNPs = NA)[-1]

rm(LeadVariants)
gc()

message(paste0("Iteratively finding overlapping variants in ", length(batches), " batches..."))

for (i in 1:length(batches)){

message(paste0("Overlapping batch ", i, "..."))
if(i < length(batches) & length(batches) > 1 & i != 1){
        temp_Lead2 <- Lead2[batches[i]:batches[i+1]]
        message("Multiple batches")
    } else if(length(batches) > 1 & i == 1){
        temp_Lead2 <- Lead2[batches[i]]
        message("Multiple batches")
    } else if(length(batches) > 1 & i == length(batches)){
        temp_Lead2 <- Lead2[batches[i]]
        message("Multiple batches")
    } else if(length(batches) == 1 & i == 1){
        temp_Lead2 <- Lead2
        message("Single batch")
    } else {
        message("Debug!")
        }

ref_temp <- ref[chr %in% temp_Lead2$chr]

temp_Lead2 <- data.table(
phenotype = temp_Lead2$phenotype,
lead_SNP = temp_Lead2$SNP, 
chr = temp_Lead2$chr, 
start = temp_Lead2$pos - args$cis_win,
end = temp_Lead2$pos + args$cis_win,
ea = temp_Lead2$ea, 
nea = temp_Lead2$nea, 
beta = temp_Lead2$beta,
se = temp_Lead2$se,
type = temp_Lead2$type
)

setkey(temp_Lead2, chr, start, end)

ref_temp$chr <- as.numeric(ref_temp$chr)


overlap_temp <- foverlaps(ref_temp, temp_Lead2, type = "within", nomatch = NULL)

overlap_temp <- overlap_temp[, .(SNPs = list(unique(SNP))), by = lead_SNP]
res <- rbind(res, overlap_temp)

message(paste0("Done"))

}

res <- unique(as.data.frame(res))
res <- as.data.table(res)
res <- merge(as.data.table(Lead2), res, by.x = "SNP", by.y = "lead_SNP")

if(nrow(res[duplicated(res$cis_gene), ]) == 0){

message("Saving results...")
split_res <- split(res, by = "phenotype")

write_gene_files <- function(phenotype, res) {
  # Construct the file name
  file_name <- paste0(phenotype, ".txt.gz")
  # Write the data.table to a file
  fwrite(res, file = file_name, sep = "\t")
}

invisible(lapply(names(split_res), function(phenotype) {
  write_gene_files(phenotype, split_res[[phenotype]])
}))

message("Saving results...done!")

} else {message("Debug!")}
