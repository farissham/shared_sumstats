//
// SBAYESRC: genome-wide Bayesian SNP-effect re-estimation from the harmonised hub.
// hub -> .ma -> tidy -> impute -> sbayesrc. Consumes the eigen-LD reference
// (ld_dir) and an optional BaselineLD annotation file.
//

include { HUB_TO_MA       } from '../../../modules/local/hub_to_ma'
include { SBAYESRC_TIDY   } from '../../../modules/local/sbayesrc_tidy'
include { SBAYESRC_IMPUTE } from '../../../modules/local/sbayesrc_impute'
include { SBAYESRC_MAIN   } from '../../../modules/local/sbayesrc_main'

workflow SBAYESRC {

    take:
    ch_harmonised   // channel: [ meta, harmonised.tsv.gz ]
    ld_dir          // value:   LD reference directory (eigen blocks + snp.info)
    annot           // value:   BaselineLD annotation file, or [] for none

    main:
    ch_versions = Channel.empty()

    HUB_TO_MA(ch_harmonised)
    ch_versions = ch_versions.mix(HUB_TO_MA.out.versions)

    SBAYESRC_TIDY(HUB_TO_MA.out.ma, ld_dir)
    ch_versions = ch_versions.mix(SBAYESRC_TIDY.out.versions)

    SBAYESRC_IMPUTE(SBAYESRC_TIDY.out.tidy, ld_dir)
    ch_versions = ch_versions.mix(SBAYESRC_IMPUTE.out.versions)

    SBAYESRC_MAIN(SBAYESRC_IMPUTE.out.imp, ld_dir, annot)
    ch_versions = ch_versions.mix(SBAYESRC_MAIN.out.versions)

    emit:
    results  = SBAYESRC_MAIN.out.results   // [ meta, .txt, .par ]
    versions = ch_versions                 // path versions.yml
}
