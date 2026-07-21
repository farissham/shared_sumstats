FROM quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2

RUN apt-get update && apt-get install -y --no-install-recommends r-base-core && \
    Rscript -e 'install.packages(c("optparse","data.table"), repos="https://cran.r-project.org")' && \
    rm -rf /var/lib/apt/lists/*
