FROM rocker/r-ver:4.4.0

RUN install2.r --error \
        data.table \
        optparse \
        remotes && \
    Rscript -e 'remotes::install_version("susieR", version="0.14.2", repos="https://cran.r-project.org")'
