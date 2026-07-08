process MAGMA_GENESET_ANALYSIS {
    tag "$meta.id"
    label 'process_low'
    conda "bioconda::magma=1.10"
    container 'quay.io/biocontainers/magma:1.10--h9f5acd7_0'

    input:
    tuple val(meta), path(genes_raw)
    path gene_sets

    output:
    tuple val(meta), path("${meta.id}.gsa.out"),                       emit: geneset_out
    tuple val(meta), path("${meta.id}.gsa.genes.out"), optional: true, emit: geneset_genes_out
    path "versions.yml",                                               emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    magma \\
        --gene-results ${genes_raw} \\
        --set-annot ${gene_sets} \\
        --out ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magma: \$(magma 2>&1 | grep -m1 -oP 'v\\K[0-9.]+' || echo "1.10")
    END_VERSIONS
    """
}
