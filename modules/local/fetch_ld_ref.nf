process FETCH_LD_REF {
    tag "${url.tokenize('/')[-1]}"
    label 'process_single'
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'

    input:
    val url

    output:
    path "ld_ref",       emit: ld_dir
    path "versions.yml", emit: versions

    script:
    // Download + unzip via R (utils::unzip) so no extra wget/unzip container is
    // needed - the SBayesRC engine image is already pulled. Flatten a single
    // nested top-level dir so snp.info / block*.eigen.bin sit directly in ld_ref/.
    """
    Rscript -e "options(timeout=7200); download.file('${url}', 'ld_ref.zip', mode='wb', quiet=TRUE); unzip('ld_ref.zip', exdir='ld_ref_raw')"

    mkdir -p ld_ref
    if [ -f ld_ref_raw/snp.info ]; then
        mv ld_ref_raw/* ld_ref/
    else
        inner=\$(find ld_ref_raw -mindepth 1 -maxdepth 1 -type d | head -n1)
        if [ -n "\$inner" ]; then mv "\$inner"/* ld_ref/; else mv ld_ref_raw/* ld_ref/; fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')
    END_VERSIONS
    """
}
