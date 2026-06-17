process SBAYESRC_TIDY {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'

    publishDir "${params.outdir}/sbayesrc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(ma_file)
    path ld_ref_dir, stageAs: 'ld_ref'

    output:
    tuple val(meta), path("${meta.id}.tidy.ma"), emit: tidy
    path "${meta.id}.tidy.log",                  emit: log
    path "versions.yml",                         emit: versions

    script:
    """
    Rscript -e "SBayesRC::tidy(mafile='${ma_file}', LDdir='ld_ref', output='${meta.id}.tidy.ma', log2file=FALSE)" \\
        > ${meta.id}.tidy.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SBayesRC: \$(Rscript -e 'cat(as.character(packageVersion("SBayesRC")))')
    END_VERSIONS
    """
}
