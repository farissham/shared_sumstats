process SBAYESRC_MAIN {
    tag "$meta.id"
    label 'process_high'
    // Serialise the MCMC across traits: each instance peaks ~5 GB; keep within
    // the cgroup memory limit by running one at a time.
    maxForks 1
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'

    publishDir "${params.outdir}/sbayesrc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(imp_ma)
    path ld_ref_dir, stageAs: 'ld_ref'
    path annot_file

    output:
    tuple val(meta), path("${meta.id}.txt"), path("${meta.id}.par"), emit: results
    path "${meta.id}.sbayesrc.log",            emit: log
    path "${meta.id}.AnnoPerSnpHsqEnrichment", emit: annot_enrich, optional: true
    path "${meta.id}.AnnoJointProb",           emit: annot_jprob,  optional: true
    path "versions.yml",                       emit: versions

    script:
    // annot is optional: when the staged file is [] (e.g. Mode B / non-EUR with no
    // matching BaselineLD annot), drop the annot= arg and run annotation-free.
    def annot_arg = annot_file ? "annot='${annot_file}', " : ""
    """
    Rscript -e "SBayesRC::sbayesrc(mafile='${imp_ma}', LDdir='ld_ref', outPrefix='${meta.id}', ${annot_arg}tuneStep=c(0.995,0.99,0.95,0.9,0.8,0.7,0.6,0.5), log2file=FALSE)" \\
        > ${meta.id}.sbayesrc.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SBayesRC: \$(Rscript -e 'cat(as.character(packageVersion("SBayesRC")))')
    END_VERSIONS
    """
}
