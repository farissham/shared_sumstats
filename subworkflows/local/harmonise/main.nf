//
// HARMONISE: branch raw sumstats by input format, normalise each to the
// canonical harmonised-sumstats hub (<id>.harmonised.tsv.gz).
// Ported paths: gwas-ssf, gwama. Remaining formats fail loudly until ported.
//

include { HARMONISE_GWAS_SSF } from '../../../modules/local/harmonise_gwas_ssf'
include { HARMONISE_GWAMA    } from '../../../modules/local/harmonise_gwama'

workflow HARMONISE {

    take:
    ch_samplesheet   // channel: [ meta, sumstats ]
    snp_info         // value:   path to LD reference snp.info

    main:
    ch_versions = Channel.empty()

    def ch_branched = ch_samplesheet.branch {
        gwas_ssf: it[0].format == 'gwas-ssf'
        gwama:    it[0].format == 'gwama'
        other:    true
    }

    // Surface unported formats as a hard error instead of silently dropping the trait.
    ch_branched.other.map { meta, sumstats ->
        error "HARMONISE: format '${meta.format}' is not yet ported (supported: gwas-ssf, gwama) for id='${meta.id}'"
    }

    HARMONISE_GWAS_SSF(ch_branched.gwas_ssf, snp_info)
    HARMONISE_GWAMA(ch_branched.gwama, snp_info)

    ch_harmonised = HARMONISE_GWAS_SSF.out.harmonised.mix(HARMONISE_GWAMA.out.harmonised)
    ch_versions   = ch_versions
        .mix(HARMONISE_GWAS_SSF.out.versions)
        .mix(HARMONISE_GWAMA.out.versions)

    emit:
    harmonised = ch_harmonised   // [ meta, harmonised.tsv.gz ]
    versions   = ch_versions     // path versions.yml
}
