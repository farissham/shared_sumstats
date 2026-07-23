//
// COLOC_ANALYSIS: Bayesian colocalization of harmonised GWAS traits against a
// set of reference datasets (eQTL/pQTL) at user-defined loci.
//
// Four ref modes (checked in order):
//   1. coloc_ref_dir              : pre-filtered per-gene *.txt files (fastest)
//   2. coloc_ref + coloc_genes    : full eQTL + explicit gene list → filter in pipeline
//   3. coloc_ref + auto_genes=true: full eQTL → genes auto-discovered per locus;
//                                   each gene is tested only at its discovered locus
//                                   (not fanned across all loci)
//   4. coloc_ref (alone)          : single gene file (backward-compat)
//

include { DISCOVER_GENES_IN_LOCUS } from '../../../modules/local/discover_genes_in_locus'
include { FILTER_EQTL_GENE        } from '../../../modules/local/filter_eqtl_gene'
include { COLOC                   } from '../../../modules/local/coloc'

workflow COLOC_ANALYSIS {

    take:
    ch_harmonised   // channel: [ meta, harmonised.tsv.gz ]
    coloc_ref       // path or null: single/full eQTL file
    coloc_ref_dir   // path or null: directory of pre-filtered per-gene *.txt files
    coloc_genes     // path or null: explicit gene-list file (one gene per line)
    loci_csv        // path: loci CSV (columns: chr, start, end)

    main:
    ch_versions = Channel.empty()

    // Parse loci synchronously — Channel.value() never closes so late-arriving
    // gene channels always get the full loci list.
    def loci_rows = file(loci_csv, checkIfExists: true)
        .readLines()
        .drop(1)
        .findAll { it.trim() }
        .collect { line ->
            def f = line.split(',')
            [ f[0].trim() as int, f[1].trim() as int, f[2].trim() as int ]
        }

    // -------------------------------------------------------------------------
    // Build ch_coloc_in per ref mode
    // -------------------------------------------------------------------------
    def ch_coloc_in

    if (coloc_ref_dir) {
        // Mode 1: pre-filtered files already in a directory → fan gene × all loci
        def ch_refs = Channel.fromPath("${coloc_ref_dir}/*.txt", checkIfExists: true)
            .map { f -> [ f.baseName, f ] }

        ch_coloc_in = ch_harmonised
            .combine(ch_refs)
            .flatMap { meta, hub, gene_id, ref ->
                loci_rows.collect { locus ->
                    [ meta, gene_id, hub, locus[0], locus[1], locus[2], ref ]
                }
            }

    } else if (params.coloc_auto_genes) {
        // Mode 3: discover genes per locus, COLOC only at the discovered locus.
        // Preserves (gene, locus) pairing so genes are not fanned across all loci.
        def ch_eqtl   = Channel.fromPath(coloc_ref, checkIfExists: true).first()
        def eqtl_path = file(coloc_ref.toString(), checkIfExists: true)

        def ch_discovery_in = ch_harmonised
            .map { meta, hub -> meta }
            .flatMap { meta ->
                loci_rows.collect { locus -> [ meta, locus[0], locus[1], locus[2], eqtl_path ] }
            }

        DISCOVER_GENES_IN_LOCUS(ch_discovery_in)
        ch_versions = ch_versions.mix(DISCOVER_GENES_IN_LOCUS.out.versions.first())

        // Flatten to (meta, gene, chr, start, end) — keeps locus info per gene
        def ch_gene_locus = DISCOVER_GENES_IN_LOCUS.out.gene_list
            .flatMap { meta, chr, start, end, f ->
                f.readLines().collect { it.trim() }.findAll { it }
                    .collect { gene -> [ meta, gene, chr, start, end ] }
            }

        // FILTER runs once per unique gene name (same eQTL file for all traits/loci)
        def ch_unique_genes = ch_gene_locus
            .map  { meta, gene, chr, start, end -> gene }
            .unique()

        FILTER_EQTL_GENE(ch_unique_genes, ch_eqtl)
        ch_versions = ch_versions.mix(FILTER_EQTL_GENE.out.versions.first())

        // Join gene-locus pairs → filtered refs → harmonised hub
        ch_coloc_in = ch_gene_locus
            .map    { meta, gene, chr, start, end -> [ gene, meta, chr, start, end ] }
            .combine(FILTER_EQTL_GENE.out.filtered, by: [0])
            .map    { gene, meta, chr, start, end, ref -> [ meta, gene, chr, start, end, ref ] }
            .combine(ch_harmonised, by: [0])
            .map    { meta, gene, chr, start, end, ref, hub ->
                [ meta, gene, hub, chr, start, end, ref ]
            }

    } else if (coloc_genes) {
        // Mode 2: explicit gene list → filter full eQTL per gene, fan to all loci
        def ch_eqtl  = Channel.fromPath(coloc_ref, checkIfExists: true).first()
        def ch_genes = Channel.fromPath(coloc_genes, checkIfExists: true)
            .splitText()
            .map  { it.trim() }
            .filter { it && !it.startsWith('#') }

        FILTER_EQTL_GENE(ch_genes, ch_eqtl)
        ch_versions = ch_versions.mix(FILTER_EQTL_GENE.out.versions.first())

        ch_coloc_in = ch_harmonised
            .combine(FILTER_EQTL_GENE.out.filtered)
            .flatMap { meta, hub, gene_id, ref ->
                loci_rows.collect { locus ->
                    [ meta, gene_id, hub, locus[0], locus[1], locus[2], ref ]
                }
            }

    } else {
        // Mode 4: single pre-filtered gene file (backward-compatible) → fan to all loci
        def ch_refs = Channel.fromPath(coloc_ref, checkIfExists: true)
            .map { f -> [ f.baseName, f ] }

        ch_coloc_in = ch_harmonised
            .combine(ch_refs)
            .flatMap { meta, hub, gene_id, ref ->
                loci_rows.collect { locus ->
                    [ meta, gene_id, hub, locus[0], locus[1], locus[2], ref ]
                }
            }
    }

    COLOC(ch_coloc_in)
    ch_versions = ch_versions.mix(COLOC.out.versions)

    emit:
    summary  = COLOC.out.summary   // [ meta, gene_id, *.coloc_summary.tsv ]
    snp_pp   = COLOC.out.snp_pp   // [ meta, gene_id, *.snp_pp.tsv ]
    versions = ch_versions
}
