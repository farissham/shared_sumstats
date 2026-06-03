/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { HARMONISE } from '../subworkflows/local/harmonise'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow SUMSTATS {

    take:
    ch_samplesheet // channel: [ meta, sumstats ]

    main:
    if (!params.ld_ref_dir) {
        error "params.ld_ref_dir must be set (select an SNP-set profile, e.g. -profile hm3, or pass --ld_ref_dir)"
    }
    def ch_snp_info = channel.fromPath("${params.ld_ref_dir}/snp.info", checkIfExists: true).first()

    HARMONISE(ch_samplesheet, ch_snp_info)

    emit:
    harmonised = HARMONISE.out.harmonised
    versions   = HARMONISE.out.versions
}
