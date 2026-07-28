process COLOC {
    tag "${meta.id}:${gene_id}:chr${chr}:${start}-${end}"
    label 'process_low'

    // No reliable bioconda build for coloc>5.1.0.1 with optparse/data.table
    // bundled together (see docker/coloc/Dockerfile). Run via conda
    // (-profile conda) or the custom image below.
    conda "conda-forge::r-optparse conda-forge::r-data.table bioconda::r-coloc=5.1.0.1"
    container 'ghcr.io/farissham/coloc:5.2.3'

    input:
    tuple val(meta), val(gene_id), path(hub), val(chr), val(start), val(end), path(ref)

    output:
    tuple val(meta), val(gene_id), path("*.coloc_summary.tsv"), emit: summary
    tuple val(meta), val(gene_id), path("*.snp_pp.tsv"),        emit: snp_pp
    path "versions.yml",                                         emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = "${meta.id}_${gene_id}_${chr}_${start}_${end}"
    def n2_arg      = params.coloc_n2 != null ? "--n2 ${params.coloc_n2}" : ''
    def s1_arg      = params.coloc_s1 != null ? "--s1 ${params.coloc_s1}" : ''
    def s2_arg      = params.coloc_s2 != null ? "--s2 ${params.coloc_s2}" : ''
    def zscore_arg  = params.coloc_ref_zscore_col != null ? "--ref_zscore_col ${params.coloc_ref_zscore_col}" : ''
    """
    coloc.R \\
        --hub         ${hub} \\
        --ref         ${ref} \\
        --chr         ${chr} \\
        --start       ${start} \\
        --end         ${end} \\
        --gene_id     ${gene_id} \\
        --type1       ${params.coloc_trait1_type} \\
        --type2       ${params.coloc_ref_type} \\
        --ref_snp_col ${params.coloc_ref_snp_col} \\
        --ref_chr_col ${params.coloc_ref_chr_col} \\
        --ref_pos_col ${params.coloc_ref_pos_col} \\
        --ref_ea_col  ${params.coloc_ref_ea_col} \\
        --ref_oa_col  ${params.coloc_ref_oa_col} \\
        --ref_n_col   ${params.coloc_ref_n_col} \\
        --p1          ${params.coloc_p1} \\
        --p2          ${params.coloc_p2} \\
        --p12         ${params.coloc_p12} \\
        --min_snps    ${params.coloc_min_snps} \\
        ${n2_arg} \\
        ${s1_arg} \\
        ${s2_arg} \\
        ${zscore_arg} \\
        --out_summary ${prefix}.coloc_summary.tsv \\
        --out_snps    ${prefix}.snp_pp.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coloc: \$(Rscript -e 'cat(as.character(packageVersion("coloc")))')
    END_VERSIONS
    """
}
