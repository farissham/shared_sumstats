nextflow.enable.dsl=2

include { SUSIE_FINEMAPPING } from './subworkflows/susie_finemapping'

params.sumstats    = null
params.ref         = null
params.loci        = "loci.csv"
params.outdir      = "results"
params.susie_script = "${projectDir}/susie.R"
params.L           = 10
params.max_iter    = 1000

workflow {
    loci_ch = Channel
        .fromPath(params.loci)
        .splitCsv(header: true)
        .map { row -> tuple(row.chr.toInteger(), row.start.toInteger(), row.end.toInteger()) }

    ref_bed = Channel.value(file("${params.ref}.bed"))
    ref_bim = Channel.value(file("${params.ref}.bim"))
    ref_fam = Channel.value(file("${params.ref}.fam"))

    sumstats_ch  = Channel.value(file(params.sumstats))
    susie_script = Channel.value(file(params.susie_script))

    SUSIE_FINEMAPPING(loci_ch, ref_bed, ref_bim, ref_fam, sumstats_ch, susie_script)
}
