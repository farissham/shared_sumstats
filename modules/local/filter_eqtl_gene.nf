process FILTER_EQTL_GENE {
    tag "${gene}"
    label 'process_low'

    conda "bioconda::r-coloc=5.1.0.1 conda-forge::r-optparse"
    container 'quay.io/biocontainers/r-coloc:5.1.0.1--r42h3121a25_1'

    input:
    val  gene
    path eqtl

    output:
    tuple val(gene), path("${gene}.eqtl.txt"), emit: filtered
    path "versions.yml",                        emit: versions

    script:
    def gene_col = params.coloc_ref_gene_col
    def zcat_cmd = eqtl.name.endsWith('.gz') ? 'zcat' : 'cat'
    if (!params.coloc_ref_parse_variant_id)
        // eQTLGen-style: exact match on gene column, all columns kept
        """
        ${zcat_cmd} ${eqtl} | awk -F'\\t' -v gene="${gene}" -v col="${gene_col}" '
            NR==1 { for(i=1;i<=NF;i++) if(\$i==col) gc=i; print; next }
            \$gc==gene' > ${gene}.eqtl.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(awk --version 2>&1 | head -1)
        END_VERSIONS
        """
    else
        // GTEx-style: prefix match on gene_id column + parse variant_id
        // variant_id: chr7_128830322_G_A_b38  →  chr=7 pos=.. oa=G ea=A
        """
        ${zcat_cmd} ${eqtl} | awk -F'\\t' -v OFS='\\t' -v gene="${gene}" -v col="${gene_col}" '
            NR==1 {
                for(i=1;i<=NF;i++) if(\$i==col) gc=i
                print "variant_id","chr","pos","ea","oa","beta","se","maf","gene_id"
                next
            }
            \$gc ~ gene {
                split(\$1, v, "_")
                gsub("chr","",v[1])
                print \$1, v[1]+0, v[2]+0, v[4], v[3], \$8, \$9, \$6, \$2
            }' > ${gene}.eqtl.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(awk --version 2>&1 | head -1)
        END_VERSIONS
        """
}
