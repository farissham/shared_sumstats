process MAGMA_GENE_ANALYSIS {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::magma=1.10"
    container 'ghcr.io/farissham/magma:1.10'

    publishDir "${params.outdir}/magma", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(hub), path(genes_annot)
    path bed
    path bim
    path fam

    output:
    tuple val(meta), path("${meta.id}.genes.raw"), emit: genes_raw
    tuple val(meta), path("${meta.id}.genes.out"), emit: genes_out
    path "versions.yml",                           emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = bed.baseName
    """
    zcat ${hub} \\
        | awk -F'\\t' -v snp="${params.magma_snp_col}" -v pval="${params.magma_p_col}" '
            NR==1 { for(i=1;i<=NF;i++) { if(\$i==snp) sc=i; if(\$i==pval) pc=i }
                    print "SNP P"; next }
            pc && \$pc != "NA" && \$pc != "" { print \$sc, \$pc }
        ' \\
        > ${meta.id}.pval.txt

    magma \\
        --bfile ${prefix} \\
        --gene-annot ${genes_annot} \\
        --pval ${meta.id}.pval.txt use='SNP,P' N=${meta.n} \\
        --out ${meta.id} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magma: \$(magma 2>&1 | grep -m1 -oP 'v\\K[0-9.]+' || echo "1.10")
    END_VERSIONS
    """
}
