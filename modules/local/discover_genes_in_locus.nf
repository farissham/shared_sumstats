process DISCOVER_GENES_IN_LOCUS {
    tag "${meta.id}:chr${chr}:${start}-${end}"
    label 'process_low'

    conda "bioconda::r-coloc=5.2.3 conda-forge::r-optparse"
    container 'quay.io/biocontainers/r-coloc:5.2.3--r44h3121a25_0'

    input:
    tuple val(meta), val(chr), val(start), val(end), path(eqtl)

    output:
    tuple val(meta), val(chr), val(start), val(end),
          path("genes_${meta.id}_chr${chr}_${start}_${end}.txt"), emit: gene_list
    path "versions.yml",                                           emit: versions

    script:
    def gene_col = params.coloc_ref_gene_col
    def zcat_cmd = eqtl.name.endsWith('.gz') ? 'zcat' : 'cat'
    if (!params.coloc_ref_parse_variant_id)
        """
        ${zcat_cmd} ${eqtl} | awk -F'\\t' \\
            -v chr="${chr}" -v lstart="${start}" -v lend="${end}" \\
            -v chrcol="${params.coloc_ref_chr_col}" \\
            -v poscol="${params.coloc_ref_pos_col}" \\
            -v genecol="${gene_col}" '
            NR==1 {
                for(i=1;i<=NF;i++) {
                    if(\$i==chrcol) cc=i
                    if(\$i==poscol) pc=i
                    if(\$i==genecol) gc=i
                }
                next
            }
            \$cc+0==chr+0 && \$pc+0>=lstart+0 && \$pc+0<=lend+0 { genes[\$gc]=1 }
            END { for(g in genes) print g }
        ' | sort > genes_${meta.id}_chr${chr}_${start}_${end}.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(awk --version 2>&1 | head -1)
        END_VERSIONS
        """
    else
        """
        ${zcat_cmd} ${eqtl} | awk -F'\\t' \\
            -v chr="${chr}" -v lstart="${start}" -v lend="${end}" \\
            -v genecol="${gene_col}" '
            NR==1 {
                for(i=1;i<=NF;i++) {
                    if(\$i==genecol)      gc=i
                    if(\$i=="variant_id") vc=i
                }
                next
            }
            {
                split(\$vc, v, "_")
                gsub("chr","",v[1])
                c = v[1]+0
                p = v[2]+0
            }
            c==chr+0 && p>=lstart+0 && p<=lend+0 { genes[\$gc]=1 }
            END { for(g in genes) print g }
        ' | sort > genes_${meta.id}_chr${chr}_${start}_${end}.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(awk --version 2>&1 | head -1)
        END_VERSIONS
        """
}
