process SBAYESRC_IMPUTE {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'

    publishDir "${params.outdir}/sbayesrc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(tidy_ma)
    path ld_ref_dir, stageAs: 'ld_ref'

    output:
    tuple val(meta), path("${meta.id}.imp.ma"), emit: imp
    path "${meta.id}.impute.log",               emit: log
    path "versions.yml",                        emit: versions

    script:
    """
    Rscript -e "SBayesRC::impute(mafile='${tidy_ma}', LDdir='ld_ref', output='${meta.id}.imp.ma', log2file=FALSE)" \\
        > ${meta.id}.impute.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SBayesRC: \$(Rscript -e 'cat(as.character(packageVersion("SBayesRC")))')
    END_VERSIONS
    """
}
