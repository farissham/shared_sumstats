//
// SUSIE_FINEMAP: per-locus SuSiE-RSS fine-mapping off the harmonised hub.
// LD is computed once per locus from the genotype panel (COMPUTE_LD) and reused
// across every trait; SUSIE re-aligns each trait's hub z onto the panel coding
// before fine-mapping. Ported from the susie_pipeline prototype.
//

include { COMPUTE_LD } from '../../../modules/local/compute_ld'
include { SUSIE      } from '../../../modules/local/susie'

workflow SUSIE_FINEMAP {

    take:
    ch_harmonised   // channel: [ meta, harmonised.tsv.gz ]
    genotype        // value:   staged PLINK panel files (bed/bim/fam)
    loci_csv        // path:    loci CSV (chr,start,end)

    main:
    ch_versions = Channel.empty()

    // loci -> [ chr, start, end ]
    def ch_loci = Channel.fromPath(loci_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> [ row.chr as int, row.start as int, row.end as int ] }

    // Resolve the panel triplet as value channels so they broadcast over loci.
    def ch_bed = genotype.map { fs -> fs.find { it.name.endsWith('.bed') } }.first()
    def ch_bim = genotype.map { fs -> fs.find { it.name.endsWith('.bim') } }.first()
    def ch_fam = genotype.map { fs -> fs.find { it.name.endsWith('.fam') } }.first()

    // LD per locus (panel-only; independent of trait).
    COMPUTE_LD(ch_loci, ch_bed, ch_bim, ch_fam)
    ch_versions = ch_versions.mix(COMPUTE_LD.out.versions)

    // Fan every trait's hub across every locus LD -> fine-map each pair.
    def ch_susie_in = ch_harmonised
        .combine(COMPUTE_LD.out.ld)
        .map { meta, hub, chr, start, end, ld, bim -> [ meta, hub, chr, start, end, ld, bim ] }

    SUSIE(ch_susie_in)
    ch_versions = ch_versions.mix(SUSIE.out.versions)

    emit:
    credible_sets = SUSIE.out.credible_sets   // [ meta, *.credible_sets.csv ]
    align         = SUSIE.out.align           // [ meta, *.align.tsv ]
    versions      = ch_versions
}
