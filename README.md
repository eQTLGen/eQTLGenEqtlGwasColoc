# Pipeline for running colocalisation analysis between eQTLGen eQTL and UKBB GWAS summary statistics

**TBA**

This pipeline runs eQTL summary statistics from eQTLGen phase 2 project and runs colocalisation analysis (HyprColoc) for every eQTL locus against corresponding loci in the set of UKBB GWAS summary statistic.

## Usage information

### Requirements for the system

- Have access to HPC with multiple cores.
- Have Bash >=3.2 installed.
- Have Java >=11 installed.
- Have Slurm scheduler managing the jobs in the HPC.
- HPC has Singularity installed and running.

### Setup of the pipeline

You can either clone it by using git (if available in HPC):

`git clone https://github.com/eQTLGen/eQTLGenEqtlGwasColoc.git`

Or just download this from the github download link and unzip.

### Required inputs

#### eQTLs

`--sig_eqtls` Gzipped text file with filtered SNP-gene combinations, showing some significance (e.g. P<5e-8). Output of [ExtractMetaAnalysisResults](https://github.com/eQTLGen/ExtractMetaAnalysisResults) pipeline. 

`--eqtl_files` Folder of genome-wide eQTL meta-analysis parquet files. Output of [MetaAnalysis](https://github.com/eQTLGen/MetaAnalysis) pipeline.

`--allele_info` eQTLGen SNP reference file in .parquet format, containing SNP ID, chr, pos, and alleles.

#### GWASs

`--gwas_files` Folder with UKBB GWAS summary statistics. These need to be in parquet format (hive structure) and harmonised to match eQTLGen eQTL summary statistics (IDs matched, genomic positions in hg38 and effect directions harmonised to match eQTLGen alternative allele).

`--gwas_manifest` [UKBB GWAS annotation file](https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/biomarkers.both_sexes.tsv.bgz), as provided by Neale Lab (https://www.nealelab.is/uk-biobank).

#### Annotation

`--gtf` ENSEMBL .gtf file for annotating gene symbols and cis/trans effects. In eQTLGen p2 project we use [v106 version](https://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz).

### Optional inputs

`--OutputDir` Optional output folder path. Defaults to `results` inside the pipeline directory.

### Additional settings

`--p_thresh` Additional P-value filter for significant eQTL effects. Defaults to 5e-8.

`--minN_thresh` Minimal sample size threshold for including genes and variants to the analysis. Defaults to 0 (no filter).

`--maxN_thresh` Variant filter for including eQTL per-locus variants into the analysis. Defaults to 0.5, meaning that variants having sample size <50% of the maximal sample size in this locus are excluded from the analysis. This is because eQTLGen data entials meta-analysis: some variants are tested in only part of the cohorts. This filter is meant to partly address the potential issues with emerging from different power for different variants. 

`--i2_thresh` Meta-analysis heterogeneity I2 threshold. Defaults 40 (<40%).

`--leadvar_window`  Distance for defining the lead variants and loci by distance clumping. Defaults to 1000000 (+/-1Mb from lead variant).

`--cis_window` Size of the cis-window. Defaults to 1000000 (+/-1Mb from lead variant).

`--trans_window` Size of the trans-window. Defaults to 5000000: gene TSS +/-5Mb from lead variant are declared to be trans effects.

### Running the pipeline

### Outputs

## Acknowledgements

This pipeline was written by Urmo Võsa and Robert Warmerdam.

HyprColoc is the work of Christopher N. Foley and colleagues:

[Foley, C. N., Staley, J. R., Breen, P. G., Sun, B. B., Kirk, P. D. W., Burgess, S., & Howson, J. M. M. (2021). A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. Nature Communications, 12(1), 1–18. https://doi.org/10.1038/s41467-020-20885-8](https://doi.org/10.1038/s41467-020-20885-8)

[Code repo](https://github.com/jrs95/hyprcoloc)
 
UK Biobank summary statistics are publicly available for research use on [Neale Lab website](https://www.nealelab.is/uk-biobank).
