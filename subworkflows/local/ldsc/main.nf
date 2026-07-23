//
// LDSC: LD Score Regression off the harmonised hub, via the real bulik/ldsc
// software. hub -> munge -> h2 (per-trait) + rg (cohort-level, all trait pairs).
//

include { LDSC_MUNGE } from '../../../modules/local/ldsc_munge'
include { LDSC_H2    } from '../../../modules/local/ldsc_h2'
include { LDSC_RG    } from '../../../modules/local/ldsc_rg'

workflow LDSC {

    take:
    ch_harmonised   // channel: [ meta, harmonised.tsv.gz ]
    ldsc_ld_dir     // value:   directory of eur_w_ld_chr LD scores
    ldsc_snplist    // value:   w_hm3.snplist for munge_sumstats --merge-alleles

    main:
    ch_versions = Channel.empty()

    LDSC_MUNGE(ch_harmonised, ldsc_snplist)
    ch_versions = ch_versions.mix(LDSC_MUNGE.out.versions)

    LDSC_H2(LDSC_MUNGE.out.sumstats, ldsc_ld_dir)
    ch_versions = ch_versions.mix(LDSC_H2.out.versions)

    // rg is cross-trait: collect every trait's munged sumstats and run once for
    // the whole cohort. LDSC_RG itself handles the <2-trait degenerate case.
    LDSC_RG(
        LDSC_MUNGE.out.sumstats.map { meta, f -> f }.collect(),
        LDSC_MUNGE.out.sumstats.map { meta, f -> meta.id }.collect(),
        ldsc_ld_dir
    )
    ch_versions = ch_versions.mix(LDSC_RG.out.versions)

    emit:
    h2       = LDSC_H2.out.summary   // [ meta, id.h2.tsv ]
    rg       = LDSC_RG.out.summary   // cohort.rg_all_pairs.tsv
    versions = ch_versions
}
