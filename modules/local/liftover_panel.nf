process LIFTOVER_PANEL {
    tag "$target_build"
    label 'process_low'
    container 'quay.io/biocontainers/ucsc-liftover:469--h9b8f530_0'

    publishDir "${params.outdir}/liftover", mode: params.publish_dir_mode

    input:
    tuple val(target_build), path(chain)
    path snp_info

    output:
    tuple val(target_build), path("${target_build}.rsid_map.tsv"), emit: rsid_map
    path "versions.yml",                                           emit: versions

    script:
    """
    panel_liftover.sh ${snp_info} ${chain} ${target_build}

    printf '"%s":\\n    ucsc-liftover: %s\\n' \\
        "${task.process}" \\
        "\$(liftOver 2>&1 | grep -oiE 'version [0-9]+' | grep -oE '[0-9]+' | head -1 || echo NA)" \\
        > versions.yml
    """
}
