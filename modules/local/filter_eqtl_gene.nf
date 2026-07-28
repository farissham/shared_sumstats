process FILTER_EQTL_GENE {
    tag "${gene}"
    label 'process_low'

    // No reliable bioconda build for coloc>5.1.0.1 with optparse/data.table
    // bundled together (see docker/coloc/Dockerfile). Run via conda
    // (-profile conda) or the custom image below.
    conda "conda-forge::r-optparse conda-forge::r-data.table bioconda::r-coloc=5.1.0.1"
    container 'ghcr.io/farissham/coloc:5.2.3'

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
        // GTEx-style: prefix match on gene_id column + parse variant_id.
        // variant_id: chr7_128830322_G_A_b38  →  chr=7 pos=.. oa=G ea=A
        // All source columns (variant_id/slope/slope_se/maf) are resolved by
        // header name, not hardcoded position - GTEx's real allpairs.txt.gz
        // layout is gene_id,variant_id,tss_distance,ma_samples,ma_count,maf,
        // pval_nominal,slope,slope_se (gene_id first, NOT variant_id first).
        """
        ${zcat_cmd} ${eqtl} | awk -F'\\t' -v OFS='\\t' -v gene="${gene}" -v col="${gene_col}" '
            NR==1 {
                for(i=1;i<=NF;i++) {
                    if(\$i==col)          gc=i
                    if(\$i=="variant_id") vc=i
                    if(\$i=="slope")      bc=i
                    if(\$i=="slope_se")   sc=i
                    if(\$i=="maf")        fc=i
                }
                print "variant_id","chr","pos","ea","oa","beta","se","maf","gene_id"
                next
            }
            \$gc ~ gene {
                split(\$vc, v, "_")
                gsub("chr","",v[1])
                print \$vc, v[1]+0, v[2]+0, v[4], v[3], \$bc, \$sc, \$fc, \$gc
            }' > ${gene}.eqtl.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(awk --version 2>&1 | head -1)
        END_VERSIONS
        """
}
