process COMPUTE_LD {
    tag "chr${chr}:${start}-${end}"
    publishDir "${params.outdir}/ld", mode: 'copy'

    input:
    tuple val(chr), val(start), val(end)
    path ref_bed
    path ref_bim
    path ref_fam

    output:
    tuple val(chr), val(start), val(end),
          path("locus_${chr}_${start}_${end}.ld"),
          path("locus_snps_${chr}_${start}_${end}.txt")

    script:
    def prefix = ref_bed.baseName
    """
    awk '\$1==${chr} && \$4>=${start} && \$4<=${end} {print \$2}' ${ref_bim} \
        > locus_snps_${chr}_${start}_${end}.txt

    plink \\
          --bfile ${prefix} \\
          --extract locus_snps_${chr}_${start}_${end}.txt \\
          --r square \\
          --out locus_${chr}_${start}_${end}
    """
}
