/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { HARMONISE         } from '../subworkflows/local/harmonise'
include { SBAYESRC          } from '../subworkflows/local/sbayesrc'
include { SUSIE_FINEMAP     } from '../subworkflows/local/susie'

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

    ch_susie = Channel.empty()
    if ('susie' in mods) {
        if (!params.genotype || !params.loci) {
            error "SUSIE: --genotype (PLINK prefix) and --loci (CSV: chr,start,end) are required for the 'susie' module"
        }
        SUSIE_FINEMAP(HARMONISE.out.harmonised, PREPARE_REFERENCE.out.genotype, file(params.loci, checkIfExists: true))
        ch_susie    = SUSIE_FINEMAP.out.credible_sets
        ch_versions = ch_versions.mix(SUSIE_FINEMAP.out.versions)
    }

    emit:
    harmonised = HARMONISE.out.harmonised
    sbayesrc   = ch_sbayesrc
    susie      = ch_susie
    versions   = ch_versions
}
