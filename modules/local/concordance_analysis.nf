process CONCORDANCE_ANALYSIS {
    label 'process_medium'
    container 'ghcr.io/farissham/concordance:1.0.0'

    input:
    path hubs   // list of harmonised hub files, unique basenames
    val  ids

    output:
    path "concordance_all_pairs.tsv", emit: summary
    path "versions.yml",              emit: versions

    script:
    def hubs_arg = hubs instanceof List ? hubs.join(',') : hubs
    def ids_arg  = ids.join(',')
    """
    concordance_analysis.py \\
        --hubs ${hubs_arg} \\
        --ids ${ids_arg} \\
        --pval_thresh ${params.concordance_pval_thresh} \\
        --out_summary concordance_all_pairs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
