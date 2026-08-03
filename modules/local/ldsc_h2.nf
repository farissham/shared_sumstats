process LDSC_H2 {
    tag "$meta.id"
    label 'process_medium'
    container 'ghcr.io/farissham/ldsc:1.0.1'

    input:
    tuple val(meta), path(sumstats)
    path ld_dir, stageAs: 'ld_ref'

    output:
    tuple val(meta), path("${meta.id}.h2.tsv"), emit: summary
    // optional: note_exit() can return before ldsc.py ever runs (e.g. it produced
    // no log at all), in which case this file never gets created.
    path "${meta.id}.h2.log",                   emit: log, optional: true
    path "versions.yml",                        emit: versions

    script:
    // Both prevalences come from the samplesheet (per-trait: different cohort subsets
    // of the same disease can have very different case fractions in their own GWAS).
    // A blank samplesheet cell for a non-nullable "number"-typed optional schema
    // field surfaces here as [] (empty list), not null, depending on the
    // nf-schema samplesheet-to-channel conversion - `!= null` alone doesn't catch
    // that, and interpolating an empty list gives the literal string "[]" on the
    // CLI, which ldsc.py then fails to parse as a float. Require an actual Number
    // instead, regardless of what shape "absent" takes.
    def samp_prev_arg = (meta.samp_prev instanceof Number) ? "--samp_prev ${meta.samp_prev}" : ''
    def pop_prev_arg  = (meta.pop_prev  instanceof Number) ? "--pop_prev ${meta.pop_prev}"   : ''
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
