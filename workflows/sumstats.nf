/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { HARMONISE         } from '../subworkflows/local/harmonise'
include { SBAYESRC          } from '../subworkflows/local/sbayesrc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow SUMSTATS {

    take:
    ch_samplesheet // channel: [ meta, sumstats ]

    main:
    ch_versions = Channel.empty()

    // Resolve the LD reference (Mode A 'prebuilt' | Mode B 'build'). snp_info gates
    // HARMONISE; ld_dir/annot feed SBayesRC; genotype feeds SuSiE.
    PREPARE_REFERENCE()
    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    HARMONISE(ch_samplesheet, PREPARE_REFERENCE.out.snp_info)
    ch_versions = ch_versions.mix(HARMONISE.out.versions)

    // Analysis modules off the harmonised hub, gated by params.modules.
    def mods = (params.modules ?: '').tokenize(',')*.trim()

    ch_sbayesrc = Channel.empty()
    if ('sbayesrc' in mods) {
        SBAYESRC(HARMONISE.out.harmonised, PREPARE_REFERENCE.out.ld_dir, PREPARE_REFERENCE.out.annot)
        ch_sbayesrc = SBAYESRC.out.results
        ch_versions = ch_versions.mix(SBAYESRC.out.versions)
    }

    emit:
    harmonised = HARMONISE.out.harmonised
    sbayesrc   = ch_sbayesrc
    versions   = ch_versions
}
