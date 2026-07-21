/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { HARMONISE         } from '../subworkflows/local/harmonise'
include { SBAYESRC          } from '../subworkflows/local/sbayesrc'
include { SUSIE_FINEMAP     } from '../subworkflows/local/susie'
include { MAGMA             } from '../subworkflows/local/magma'
include { LDSC              } from '../subworkflows/local/ldsc'

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

    ch_magma_genes      = Channel.empty()
    ch_magma_geneset    = Channel.empty()
    ch_magma_tissue     = Channel.empty()
    if ('magma' in mods) {
        if (!params.genotype || !params.magma_gene_loc) {
            error "MAGMA: --genotype (PLINK prefix) and --magma_gene_loc are required for the 'magma' module"
        }
        def gene_loc_ch  = file(params.magma_gene_loc,  checkIfExists: true)
        def gene_sets_ch = params.magma_gene_sets ? file(params.magma_gene_sets, checkIfExists: true) : []
        def gtex_ch      = params.magma_gtex      ? file(params.magma_gtex,      checkIfExists: true) : []
        MAGMA(HARMONISE.out.harmonised, PREPARE_REFERENCE.out.genotype, gene_loc_ch, gene_sets_ch, gtex_ch)
        ch_magma_genes   = MAGMA.out.genes_out
        ch_magma_geneset = MAGMA.out.geneset_out
        ch_magma_tissue  = MAGMA.out.tissue_out
        ch_versions      = ch_versions.mix(MAGMA.out.versions)
    }

    ch_ldsc_h2 = Channel.empty()
    ch_ldsc_rg = Channel.empty()
    if ('ldsc' in mods) {
        if (!params.ldsc_ld_dir || !params.ldsc_snplist) {
            error "LDSC: --ldsc_ld_dir and --ldsc_snplist are required for the 'ldsc' module"
        }
        LDSC(HARMONISE.out.harmonised,
             file(params.ldsc_ld_dir, checkIfExists: true),
             file(params.ldsc_snplist, checkIfExists: true))
        ch_ldsc_h2  = LDSC.out.h2
        ch_ldsc_rg  = LDSC.out.rg
        ch_versions = ch_versions.mix(LDSC.out.versions)
    }

    emit:
    harmonised     = HARMONISE.out.harmonised
    sbayesrc       = ch_sbayesrc
    susie          = ch_susie
    magma_genes    = ch_magma_genes
    magma_geneset  = ch_magma_geneset
    magma_tissue   = ch_magma_tissue
    ldsc_h2        = ch_ldsc_h2
    ldsc_rg        = ch_ldsc_rg
    versions       = ch_versions
}
