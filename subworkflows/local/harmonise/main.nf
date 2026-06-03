//
// HARMONISE: branch raw sumstats by input format, normalise each to the
// canonical harmonised-sumstats hub (<id>.harmonised.tsv.gz).
// Phase 1 implements the gwas-ssf path; other formats fail loudly until ported.
//

include { HARMONISE_GWAS_SSF } from '../../../modules/local/harmonise_gwas_ssf'

workflow HARMONISE {

    take:
    ch_samplesheet   // channel: [ meta, sumstats ]
    snp_info         // value:   path to LD reference snp.info

    main:
    ch_versions = Channel.empty()

    def ch_branched = ch_samplesheet.branch {
        gwas_ssf: it[0].format == 'gwas-ssf'
        other:    true
    }

    // Phase 1 supports gwas-ssf only. Surface anything else as a hard error
    // instead of silently dropping the trait.
    ch_branched.other.map { meta, sumstats ->
        error "HARMONISE: format '${meta.format}' is not yet ported (Phase 1 = gwas-ssf only) for id='${meta.id}'"
    }

    HARMONISE_GWAS_SSF(ch_branched.gwas_ssf, snp_info)
    ch_versions = ch_versions.mix(HARMONISE_GWAS_SSF.out.versions)

    emit:
    harmonised = HARMONISE_GWAS_SSF.out.harmonised   // [ meta, harmonised.tsv.gz ]
    versions   = ch_versions                          // path versions.yml
}
