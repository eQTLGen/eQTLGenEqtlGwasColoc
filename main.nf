#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2


def helpmessage() {

log.info"""

EqtlGwasColocalisation v${workflow.manifest.version}"
===========================================================
Pipeline for running HyprColoc colocalisation analyses (https://www.nature.com/articles/s41467-020-20885-8) between cis/trans-eQTL loci from eQTLGen eQTL summary statistics and GWAS sumstats.

Usage:

nextflow run main.nf 
--eqtl_files \
--sig_eqtls \
--allele_info \
--gtf \
--OutputDir

Mandatory arguments:
--eqtl_files                eQTLGen parquet dataset.
--gwas_files                GWAS parquet dataset. Needs to be harmonised to follow same format as eQTL dataset.
--sig_eqtls                 Text file with significant eQTLGen eQTLs.
--allele_info               Parquet file with alleles and SNP positions for eQTL dataset.
--gwas_manifest             GWAS manifest file.
--gtf                       GTF file for gene annotation.

Optional arguments:
--OutputDir                 Output directory. Defaults to "results".
--leadvar_window            Window used in distance pruning to find locus lead variants. Defaults to 1000000.
--cis_window                cis-eQTL window. Defaults to 1000000.
--trans_window              trans-eQTL window. Defaults to 5000000.
--p_thresh                  P-value threshold for significant eQTL effects. Defaults to 5e-8.
--i2_thresh                 Heterogeneity threshold. Defaults to 40 (<40%).
--maxN_thresh               Per gene maximal sample size threshold to include variants. Defaults to 0.8 (SNPs with >=0.8*max(N))
--minN_thresh               Minimal sample size threshold. Defaults to 0 (no filtering).

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.OutputDir = 'results'
params.leadvar_window = 1000000
params.cis_window = 1000000
params.trans_window = 5000000
params.p_thresh = 5e-8
params.i2_thresh = 40
params.maxN_thresh = 0.8
params.minN_thresh = 0


//Show parameter values
log.info """=======================================================
eQTLGen cis-trans colocalisation pipeline v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Output directory']                         = params.OutputDir
summary['eQTL folder']                              = params.eqtl_files
summary['sig. eQTL results']                        = params.sig_eqtls
summary['GWAS folder']                              = params.gwas_files
summary['GWAS manifest file']                       = params.gwas_manifest
summary['Allele info file']                         = params.allele_info
summary['GTF file']                                 = params.gtf
summary['Pruning window']                           = params.leadvar_window
summary['cis-eQTL window']                          = params.cis_window
summary['trans-eQTL window']                        = params.trans_window
summary['P threshold']                              = params.p_thresh
summary['I2 threshold']                             = params.i2_thresh
summary['maxN threshold']                           = params.maxN_thresh
summary['minN threshold']                           = params.minN_thresh

// import modules
include { MAKELOCI; COLOC; MakeLoci; Coloc } from './modules/EqtlGwasColocalisation.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
// Get eQTL channel
empirical_ch = Channel.fromPath(params.eqtl_files, type: 'dir').ifEmpty { exit 1, "eQTL files not found!" }
gwas_ch = Channel.fromPath(params.gwas_files, type: 'dir').ifEmpty { exit 1, "GWAS files not found!" }
gwas_manifest_ch = Channel.fromPath(params.gwas_manifest, type: 'file').ifEmpty { exit 1, "GWAS manifest file not found!" }
sig_ch = Channel.fromPath(params.sig_eqtls).ifEmpty { exit 1, "Sig. eQTL files not found!" }
allele_ch = Channel.fromPath(params.allele_info).ifEmpty { exit 1, "eQTLGen reference not found!" }
gtf_ch = Channel.fromPath(params.gtf).ifEmpty { exit 1, "GTF not found!" }

input_ch = sig_ch.combine(empirical_ch).combine(gwas_ch).combine(gwas_manifest_ch).combine(allele_ch).combine(gtf_ch)

leadvar_window = Channel.value(params.leadvar_window)
cis_window = Channel.value(params.cis_window)
trans_window = Channel.value(params.trans_window)
p_thresh = Channel.value(params.p_thresh)
i2_thresh = Channel.value(params.i2_thresh)
maxN_thresh = Channel.value(params.maxN_thresh)
minN_thresh = Channel.value(params.minN_thresh)

input_ch = input_ch.combine(leadvar_window).combine(cis_window).combine(trans_window)
.combine(p_thresh).combine(i2_thresh).combine(maxN_thresh).combine(minN_thresh)


workflow {
        MAKELOCI(input_ch)
        
        coloc_input_ch = input_ch.combine(MAKELOCI.out.flatten())
        
        COLOC(coloc_input_ch)
        COLOC.out.collectFile(name: 'EqtlGwasColocResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
        }


workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
