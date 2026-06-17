//
// PREPARE_REFERENCE: resolve the LD reference under the selected mode and emit a
// uniform set of channels for the rest of the pipeline, regardless of mode.
//   Mode A ('prebuilt'): a pre-staged ld_ref_dir (ld_ref_url download = Task 1.3).
//   Mode B ('build'):    build eigen-LD + snp.info from a genotype panel via
//                        SBayesRC::LDstep (Stage 3 - not yet implemented).
// snp.info gates HARMONISE; ld_dir feeds SBayesRC; genotype feeds SuSiE.
//

include { FETCH_LD_REF } from '../../../modules/local/fetch_ld_ref'

workflow PREPARE_REFERENCE {

    main:
    ch_versions = Channel.empty()

    def mode = params.ld_mode ?: 'prebuilt'

    def ch_ld_dir
    if (mode == 'prebuilt') {
        if (params.ld_ref_dir && file(params.ld_ref_dir).exists()) {
            ch_ld_dir = Channel.fromPath(params.ld_ref_dir, type: 'dir', checkIfExists: true).first()
        } else if (params.ld_ref_url) {
            FETCH_LD_REF(Channel.value(params.ld_ref_url))
            ch_ld_dir   = FETCH_LD_REF.out.ld_dir.first()
            ch_versions = ch_versions.mix(FETCH_LD_REF.out.versions)
        } else {
            error "PREPARE_REFERENCE: ld_mode='prebuilt' requires a pre-staged --ld_ref_dir (select an SNP-set profile, e.g. -profile hm3)"
        }
    } else if (mode == 'build') {
        error "PREPARE_REFERENCE: ld_mode='build' (Mode B) is not implemented yet (Stage 3)"
    } else {
        error "PREPARE_REFERENCE: unknown ld_mode='${mode}' (expected 'prebuilt' or 'build')"
    }

    // Derived, mode-agnostic outputs.
    def ch_snp_info = ch_ld_dir.map { dir -> file("${dir}/snp.info", checkIfExists: true) }
    def ch_annot    = params.annot_file ? Channel.fromPath(params.annot_file, checkIfExists: true).first() : Channel.value([])
    def ch_genotype = params.genotype  ? Channel.fromPath("${params.genotype}*", checkIfExists: true).collect() : Channel.value([])

    emit:
    ld_dir   = ch_ld_dir      // value: LD reference directory (eigen blocks + snp.info)
    snp_info = ch_snp_info    // value: path to snp.info (gates HARMONISE)
    annot    = ch_annot       // value: BaselineLD annotation file, or []
    genotype = ch_genotype    // value: staged PLINK panel files, or []
    versions = ch_versions    // path versions.yml
}
