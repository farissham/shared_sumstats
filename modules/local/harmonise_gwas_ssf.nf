process HARMONISE_GWAS_SSF {
    tag "$meta.id"
    label 'process_low'
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'   // TODO Phase 4: shared namespace + pinned tag

    publishDir "${params.outdir}/harmonise", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sumstats)
    path snp_info

    output:
    tuple val(meta), path("${meta.id}.harmonised.tsv.gz"), emit: harmonised
    path "versions.yml",                                   emit: versions

    script:
    """
    harmonise_gwas_ssf.R \\
        --input ${sumstats} \\
        --snp_info ${snp_info} \\
        --n ${meta.n} \\
        --output ${meta.id}.harmonised.tsv.gz

    printf '"%s":\\n    r-base: %s\\n    data.table: %s\\n' \\
        "${task.process}" \\
        "\$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')" \\
        "\$(Rscript -e 'cat(as.character(packageVersion("data.table")))')" \\
        > versions.yml
    """
}
