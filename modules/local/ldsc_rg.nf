process LDSC_RG {
    label 'process_medium'
    container 'ghcr.io/farissham/ldsc:dev'

    input:
    path sumstats   // list of [meta.id].munge.sumstats.gz — unique basenames, no staging collision
    val  ids
    path ld_dir, stageAs: 'ld_ref'

    output:
    path "cohort.rg_all_pairs.tsv", emit: summary
    path "cohort.rg.log",           emit: log
    path "versions.yml",            emit: versions

    script:
    def sumstats_arg = sumstats instanceof List ? sumstats.join(',') : sumstats
    def ids_arg       = ids.join(',')
    """
    ldsc_rg.R \\
        --sumstats ${sumstats_arg} \\
        --ids ${ids_arg} \\
        --ld_dir ld_ref \\
        --out_summary cohort.rg_all_pairs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')
    END_VERSIONS
    """
}
