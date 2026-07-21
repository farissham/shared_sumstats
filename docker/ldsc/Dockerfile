# bulik/ldsc (munge_sumstats.py, ldsc.py) + R (optparse, data.table) for LDSC_H2/LDSC_RG's log-parsing wrapper scripts.
FROM quay.io/biocontainers/ldsc:1.0.1--pyhdfd78af_2

RUN apt-get update && apt-get install -y --no-install-recommends r-base-core && \
    Rscript -e 'install.packages(c("optparse","data.table"), repos="https://cran.r-project.org")' && \
    rm -rf /var/lib/apt/lists/*
