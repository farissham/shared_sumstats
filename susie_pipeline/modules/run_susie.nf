process RUN_SUSIE {
    tag "chr${chr}:${start}-${end}"
    publishDir "${params.outdir}/susie", mode: 'copy'

    input:
    tuple val(chr), val(start), val(end), path(ld_matrix), path(snp_list)
    path sumstats
    path susie_script

    output:
    path "credible_sets_${chr}_${start}_${end}.csv"
    path "pip_plot_${chr}_${start}_${end}.png"

    script:
    """
    Rscript ${susie_script} \\
        --sumstats ${sumstats} \\
        --ld       ${ld_matrix} \\
        --snps     ${snp_list} \\
        --chr      ${chr} \\
        --start    ${start} \\
        --end      ${end} \\
        --L        ${params.L} \\
        --max_iter ${params.max_iter} \\
        --out_cs   credible_sets_${chr}_${start}_${end}.csv \\
        --out_plot pip_plot_${chr}_${start}_${end}.png
    """
}
