//
// HARMONISE: normalise raw sumstats of any supported format to the canonical
// harmonised-sumstats hub (<id>.harmonised.tsv.gz). Format-specific parsing,
// the build matcher, and allele orientation all live in bin/harmonise.R; this
// layer just validates the format and dispatches.
//

include { HARMONISE_SUMSTATS } from '../../../modules/local/harmonise'

workflow HARMONISE {

    take:
    ch_samplesheet   // channel: [ meta, sumstats ]
    snp_info         // value:   path to LD reference snp.info

    main:
    ch_versions = Channel.empty()

    def supported = ['gwas-ssf', 'gwama', 'hail', 'cvdkp']
    def ch_in = ch_samplesheet.map { meta, sumstats ->
        if (!(meta.format in supported))
            error "HARMONISE: unsupported format '${meta.format}' for id='${meta.id}' (supported: ${supported.join(', ')})"
        [ meta, sumstats ]
    }

    HARMONISE_SUMSTATS(ch_in, snp_info)
    ch_versions = ch_versions.mix(HARMONISE_SUMSTATS.out.versions)

    emit:
    harmonised = HARMONISE_SUMSTATS.out.harmonised   // [ meta, harmonised.tsv.gz ]
    versions   = ch_versions                          // path versions.yml
}
