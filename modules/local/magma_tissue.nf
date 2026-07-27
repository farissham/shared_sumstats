process MAGMA_TISSUE {
    tag "$meta.id"
    label 'process_low'
    conda "bioconda::magma=1.10"
    container 'ghcr.io/farissham/magma:1.10'

    input:
    tuple val(meta), path(genes_raw)
    path gtex

    output:
    tuple val(meta), path("${meta.id}.tissue.gsa.out"), emit: tissue_out
    path "versions.yml",                                emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    magma \\
        --gene-results ${genes_raw} \\
        --gene-covar ${gtex} \\
        --out ${meta.id}.tissue \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magma: \$(magma 2>&1 | grep -m1 -oP 'v\\K[0-9.]+' || echo "1.10")
    END_VERSIONS
    """
}
