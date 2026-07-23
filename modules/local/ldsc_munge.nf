process LDSC_MUNGE {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2'

    input:
    tuple val(meta), path(hub)
    path snplist

    output:
    tuple val(meta), path("${meta.id}.munge.sumstats.gz"), emit: sumstats
    path "${meta.id}.munge.log",                           emit: log
    path "versions.yml",                                   emit: versions

    script:
    // Hub schema is fixed (rsid chr pos ea oa eaf beta se p n z), so columns are
    // hardcoded here rather than exposed as params — same spirit as HUB_TO_MA.
    """
    munge_sumstats.py \\
        --sumstats ${hub} \\
        --snp rsid --a1 ea --a2 oa --p p \\
        --signed-sumstats z,0 \\
        --N-col n \\
        --merge-alleles ${snplist} \\
        --out ${meta.id}.munge

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ldsc: \$(munge_sumstats.py --help 2>&1 | grep -m1 -oE 'v[0-9.]+' || echo "1.0.1")
    END_VERSIONS
    """
}
