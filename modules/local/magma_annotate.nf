process MAGMA_ANNOTATE {
    tag "magma_annotate"
    label 'process_low'
    conda "bioconda::magma=1.10"
    container 'ghcr.io/farissham/magma:1.10'

    publishDir "${params.outdir}/magma", mode: params.publish_dir_mode

    input:
    path bim      // PLINK .bim from reference panel (SNP positions)
    path gene_loc // NCBI gene location file (ENTREZ_ID CHR START STOP STRAND SYMBOL)

    output:
    path "panel.genes.annot", emit: genes_annot
    path "versions.yml",      emit: versions

    script:
    def annotate_arg = (params.magma_window as int) > 0 ? "--annotate window=${params.magma_window}" : "--annotate"
    """
    awk '{print \$2, \$1, \$4}' ${bim} > panel.snp.loc

    magma \\
        ${annotate_arg} \\
        --snp-loc panel.snp.loc \\
        --gene-loc ${gene_loc} \\
        --out panel

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magma: \$(magma 2>&1 | head -n2 | grep -oP 'v\\K[0-9]+\\.[0-9]+' | head -n1 || echo "1.10")
    END_VERSIONS
    """
}
