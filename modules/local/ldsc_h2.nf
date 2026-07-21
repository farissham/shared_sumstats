process LDSC_H2 {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/farissham/ldsc:1.0.1'

    input:
    tuple val(meta), path(sumstats)
    path ld_dir, stageAs: 'ld_ref'

    output:
    tuple val(meta), path("${meta.id}.h2.tsv"), emit: summary
    path "${meta.id}.h2.log",                   emit: log
    path "versions.yml",                        emit: versions

    script:
    // pop_prev comes from the samplesheet (per-trait, same field coloc/other modules
    // could use later); samp_prev is run-level since it isn't in the hub or samplesheet.
    def samp_prev_arg = params.ldsc_samp_prev != null ? "--samp_prev ${params.ldsc_samp_prev}" : ''
    def pop_prev_arg  = meta.pop_prev != null ? "--pop_prev ${meta.pop_prev}" : ''
    """
    ldsc_h2.R \\
        --sumstats ${sumstats} \\
        --ld_dir ld_ref \\
        --id ${meta.id} \\
        ${samp_prev_arg} \\
        ${pop_prev_arg} \\
        --out_summary ${meta.id}.h2.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')
    END_VERSIONS
    """
}
