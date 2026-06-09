process HARMONISE_SUMSTATS {
    tag "$meta.id"
    label 'process_low'
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'   // TODO Phase 4: shared namespace + pinned tag

    publishDir "${params.outdir}/harmonise", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sumstats)
    path snp_info

    output:
    tuple val(meta), path("${meta.id}.harmonised.tsv.gz"), emit: harmonised
    path "versions.yml",                                   emit: versions

    script:
    def rsid_map_arg = meta.rsid_map              ? "--rsid_map ${meta.rsid_map}" : ""
    def build_arg    = meta.build                 ? "--build ${meta.build}"       : ""
    def fwd_arg      = meta.assume_forward_strand ? "--assume_forward_strand"     : ""
    """
    harmonise.R \\
        --input ${sumstats} \\
        --format ${meta.format} \\
        --snp_info ${snp_info} \\
        --n ${meta.n} \\
        ${build_arg} \\
        --panel_build ${params.ld_ref_build} \\
        ${rsid_map_arg} \\
        ${fwd_arg} \\
        --output ${meta.id}.harmonised.tsv.gz

    printf '"%s":\\n    r-base: %s\\n    data.table: %s\\n' \\
        "${task.process}" \\
        "\$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')" \\
        "\$(Rscript -e 'cat(as.character(packageVersion("data.table")))')" \\
        > versions.yml
    """
}
