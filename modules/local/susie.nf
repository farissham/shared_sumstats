process SUSIE {
    tag "${meta.id}:chr${chr}:${start}-${end}"
    label 'process_low'

    // susieR is a conda-forge package (no bioconda biocontainer). Run via conda
    // (-profile conda) or an image that bundles r-susier + r-data.table +
    // r-optparse. TODO: pin/publish a SuSiE container and set it here.
    conda "conda-forge::r-susier=0.14.2 conda-forge::r-data.table conda-forge::r-optparse"
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'

    publishDir "${params.outdir}/susie", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hub), val(chr), val(start), val(end), path(ld), path(bim)

    output:
    tuple val(meta), path("*.credible_sets.csv"), emit: credible_sets
    tuple val(meta), path("*.align.tsv"),         emit: align
    tuple val(meta), path("*.pip.png"),           emit: pip_plot, optional: true
    path "versions.yml",                          emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.id}_${chr}_${start}_${end}"
    """
    susie.R \\
        --hub ${hub} \\
        --ld ${ld} \\
        --bim ${bim} \\
        --chr ${chr} --start ${start} --end ${end} \\
        --n ${meta.n} \\
        --l ${params.susie_l} \\
        --max_iter ${params.susie_max_iter} \\
        --out_cs ${prefix}.credible_sets.csv \\
        --out_plot ${prefix}.pip.png \\
        --out_align ${prefix}.align.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        susieR: \$(Rscript -e 'cat(tryCatch(as.character(packageVersion("susieR")), error=function(e) "NA"))')
    END_VERSIONS
    """
}
