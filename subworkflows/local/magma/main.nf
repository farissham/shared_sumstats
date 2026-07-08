//
// MAGMA: genome-wide gene-based, gene-set, and tissue-specific enrichment analysis
// from a harmonised GWAS hub. SNP-to-gene annotation runs once per reference panel;
// gene analysis, gene-set analysis, and tissue analysis run per trait.
//

include { MAGMA_ANNOTATE         } from '../../../modules/local/magma_annotate'
include { MAGMA_GENE_ANALYSIS    } from '../../../modules/local/magma_gene_analysis'
include { MAGMA_GENESET_ANALYSIS } from '../../../modules/local/magma_geneset_analysis'
include { MAGMA_TISSUE           } from '../../../modules/local/magma_tissue'

workflow MAGMA {

    take:
    ch_harmonised   // channel: [ meta, harmonised.tsv.gz ]
    genotype        // value:   staged PLINK panel files (bed/bim/fam) as a list
    gene_loc        // path:    NCBI gene location file (e.g. NCBI38.gene.loc)
    gene_sets       // path:    gene-set .gmt file; pass [] to skip gene-set analysis
    gtex            // path:    GTEx expression file; pass [] to skip tissue analysis

    main:
    ch_versions = Channel.empty()

    // Resolve the panel triplet as value channels so they broadcast over traits.
    def ch_bed = genotype.map { fs -> fs.find { it.name.endsWith('.bed') } }.first()
    def ch_bim = genotype.map { fs -> fs.find { it.name.endsWith('.bim') } }.first()
    def ch_fam = genotype.map { fs -> fs.find { it.name.endsWith('.fam') } }.first()

    // SNP → gene annotation runs once (panel-only, independent of trait).
    MAGMA_ANNOTATE(ch_bim, gene_loc)
    ch_versions = ch_versions.mix(MAGMA_ANNOTATE.out.versions)

    // Fan every trait's hub across the single annotation → gene analysis per trait.
    def ch_gene_in = ch_harmonised
        .combine(MAGMA_ANNOTATE.out.genes_annot)
        .map { meta, hub, annot -> [ meta, hub, annot ] }

    MAGMA_GENE_ANALYSIS(ch_gene_in, ch_bed, ch_bim, ch_fam)
    ch_versions = ch_versions.mix(MAGMA_GENE_ANALYSIS.out.versions)

    // Gene-set analysis: skip if gene_sets is an empty list.
    ch_geneset_out       = Channel.empty()
    ch_geneset_genes_out = Channel.empty()
    if (!(gene_sets instanceof List && gene_sets.isEmpty())) {
        MAGMA_GENESET_ANALYSIS(MAGMA_GENE_ANALYSIS.out.genes_raw, gene_sets)
        ch_geneset_out       = MAGMA_GENESET_ANALYSIS.out.geneset_out
        ch_geneset_genes_out = MAGMA_GENESET_ANALYSIS.out.geneset_genes_out
        ch_versions          = ch_versions.mix(MAGMA_GENESET_ANALYSIS.out.versions)
    }

    // Tissue-specific analysis: skip if gtex is an empty list.
    ch_tissue_out = Channel.empty()
    if (!(gtex instanceof List && gtex.isEmpty())) {
        MAGMA_TISSUE(MAGMA_GENE_ANALYSIS.out.genes_raw, gtex)
        ch_tissue_out = MAGMA_TISSUE.out.tissue_out
        ch_versions   = ch_versions.mix(MAGMA_TISSUE.out.versions)
    }

    emit:
    genes_raw         = MAGMA_GENE_ANALYSIS.out.genes_raw   // [ meta, *.genes.raw ]
    genes_out         = MAGMA_GENE_ANALYSIS.out.genes_out   // [ meta, *.genes.out ]
    geneset_out       = ch_geneset_out                      // [ meta, *.gsa.out ]
    geneset_genes_out = ch_geneset_genes_out                // [ meta, *.gsa.genes.out ]
    tissue_out        = ch_tissue_out                       // [ meta, *.tissue.gsa.out ]
    versions          = ch_versions
}
