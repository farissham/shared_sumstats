process COMPUTE_LD {
    tag "chr${chr}:${start}-${end}"
    label 'process_low'
    container 'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1'

    publishDir "${params.outdir}/susie/ld", mode: params.publish_dir_mode

    input:
    tuple val(chr), val(start), val(end)
    path bed
    path bim
    path fam

    output:
    tuple val(chr), val(start), val(end),
          path("locus_${chr}_${start}_${end}.ld"),
          path("locus_${chr}_${start}_${end}.locus.bim"), emit: ld
    path "versions.yml", emit: versions

    script:
    def prefix = bed.baseName
    // The locus .bim is awk'd straight from the panel .bim, so its row order
    // matches plink --extract's (positional) order == the .ld matrix order.
    """
    awk -v c=${chr} -v s=${start} -v e=${end} '\$1==c && \$4>=s && \$4<=e' ${bim} \\
        > locus_${chr}_${start}_${end}.locus.bim
    cut -f2 locus_${chr}_${start}_${end}.locus.bim > snps_${chr}_${start}_${end}.txt

    plink \\
        --bfile ${prefix} \\
        --extract snps_${chr}_${start}_${end}.txt \\
        --r square \\
        --ld-window 99999 \\
        --ld-window-kb 99999 \\
        --out locus_${chr}_${start}_${end}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | head -n1 | sed 's/^PLINK //; s/ .*//')
    END_VERSIONS
    """
}
