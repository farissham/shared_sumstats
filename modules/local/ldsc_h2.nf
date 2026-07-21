process LDSC_H2 {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/farissham/ldsc:dev'

    input:
    tuple val(meta), path(sumstats)
    path ld_dir, stageAs: 'ld_ref'

    output:
    tuple val(meta), path("${meta.id}.h2.tsv"), emit: summary
    path "${meta.id}.h2.log",                   emit: log
    path "versions.yml",                        emit: versions

    script:
    """
    ldsc_h2.R \\
        --sumstats ${sumstats} \\
        --ld_dir ld_ref \\
        --id ${meta.id} \\
        --out_summary ${meta.id}.h2.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')
    END_VERSIONS
    """
}
