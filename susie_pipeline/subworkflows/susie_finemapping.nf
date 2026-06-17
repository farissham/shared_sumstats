include { COMPUTE_LD } from '../modules/compute_ld'
include { RUN_SUSIE  } from '../modules/run_susie'

workflow SUSIE_FINEMAPPING {
    take:
    loci_ch       // channel: [chr, start, end]
    ref_bed       // path: reference .bed
    ref_bim       // path: reference .bim
    ref_fam       // path: reference .fam
    sumstats_ch   // path: summary statistics
    susie_script  // path: susie.R script

    main:
    ld_ch = COMPUTE_LD(loci_ch, ref_bed, ref_bim, ref_fam)
    RUN_SUSIE(ld_ch, sumstats_ch, susie_script)

    emit:
    credible_sets = RUN_SUSIE.out[0]
    pip_plots     = RUN_SUSIE.out[1]
}
