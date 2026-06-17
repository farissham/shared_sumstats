process HUB_TO_MA {
    tag "$meta.id"
    label 'process_single'
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'

    publishDir "${params.outdir}/sbayesrc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hub)

    output:
    tuple val(meta), path("${meta.id}.ma"), emit: ma
    path "versions.yml",                    emit: versions

    script:
    // Project the canonical hub (rsid chr pos ea oa eaf beta se p n z) into SBayesRC
    // COJO .ma columns (SNP A1 A2 freq b se p N). Pure column selection, no R needed.
    """
    zcat ${hub} \\
        | awk -F'\\t' -v OFS='\\t' 'NR==1 { print "SNP","A1","A2","freq","b","se","p","N"; next } { print \$1,\$4,\$5,\$6,\$7,\$8,\$9,\$10 }' \\
        > ${meta.id}.ma

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n1)
    END_VERSIONS
    """
}
