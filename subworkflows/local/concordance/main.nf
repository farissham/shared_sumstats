//
// CONCORDANCE: pairwise effect-direction concordance across every trait's
// harmonised hub, via belowlab/Concordance-Analysis (bin/concordance.py,
// vendored). Tests whether independent GWAS agree in effect direction at
// their significant loci more than expected under a permutation null -
// same "collect every trait, compare all pairs" shape as LDSC_RG, and
// gracefully no-ops with <2 traits rather than failing the run.
//

include { CONCORDANCE_ANALYSIS } from '../../../modules/local/concordance_analysis'

workflow CONCORDANCE {

    take:
    ch_harmonised   // channel: [ meta, harmonised.tsv.gz ]

    main:
    ch_versions = Channel.empty()

    CONCORDANCE_ANALYSIS(
        ch_harmonised.map { meta, f -> f }.collect(),
        ch_harmonised.map { meta, f -> meta.id }.collect()
    )
    ch_versions = ch_versions.mix(CONCORDANCE_ANALYSIS.out.versions)

    emit:
    summary  = CONCORDANCE_ANALYSIS.out.summary   // cohort-wide, not per-trait
    versions = ch_versions
}
