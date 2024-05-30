#!/bin/bash nextflow

process MakeLoci {
    //scratch true

    publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "eQtlLead*"

    input:
        tuple path(sig_res), path(eqtls), path(gwas), path(gwas_manifest), path(ref), path(gtf), val(lead_variant_win), val(cis_win), val(trans_win), val(p_thresh), val(i2), val(maxN), val(minN)

    output:
       path("ENSG*.txt.gz")

    script:
        """
        MakeLoci.R \
        --sig_res ${sig_res} \
        --eqtl_folder ${eqtls} \
        --reference ${ref} \
        --gtf ${gtf} \
        --lead_variant_win ${lead_variant_win} \
        --cis_win ${cis_win} \
        --trans_win ${trans_win} \
        --p_thresh ${p_thresh} \
        --i2_thresh ${i2} \
        --maxN_thresh ${maxN} \
        --minN_thresh ${minN}
        """
}

process Coloc {
    input:
        tuple path(sig_res), path(eqtls), path(gwas), path(gwas_manifest), path(ref), path(gtf), val(lead_variant_win), val(cis_win), val(trans_win), val(p_thresh), val(i2), val(maxN), val(minN), path(loci)

    output:
        path("*_coloc.txt")

    script:
        """
        Hyprcoloc.R \
        --i2_thresh ${i2} \
        --maxN_thresh ${maxN} \
        --minN_thresh ${minN} \
        --loci ${loci} \
        --eqtl_folder ${eqtls} \
        --gwas_folder ${gwas} \
        --pheno_manifest ${gwas_manifest} \
        --gtf ${gtf}
        """
}

workflow MAKELOCI {
    take:
        data

    main:
        loci_ch = MakeLoci(data)
        
    emit:
        loci_ch

}

workflow COLOC {
    take:
        data

    main:
        coloc_output_ch = Coloc(data)

    emit:
        coloc_output_ch
}
