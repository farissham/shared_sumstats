# quay.io/biocontainers/r-coloc only bundles exactly the coloc R package
# itself (whatever its recipe pins) - the module's `conda "bioconda::r-coloc=X
# conda-forge::r-optparse"` directive only takes effect under -profile conda;
# it does nothing for -profile singularity, which just pulls the pre-built
# container verbatim. Confirmed missing at runtime: "there is no package
# called 'optparse'" (bin/coloc.R also needs data.table, neither of which
# ship in the biocontainer). Same fix as susie/magma: build our own with
# everything bin/coloc.R actually needs bundled together.
#
# 5.2.3 is CRAN's real current version (confirmed via crandb.r-pkg.org) - the
# modules' original version number wasn't wrong, it just was never built as a
# biocontainer. Installing directly from CRAN restores the originally
# intended version instead of the 5.1.0.1 downgrade used as a stopgap.
FROM rocker/r-ver:4.4.0

RUN install2.r --error \
        data.table \
        optparse \
        coloc
