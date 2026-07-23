# R (optparse, data.table) for LDSC_H2/LDSC_RG's log-parsing wrapper scripts,
# plus bulik/ldsc (munge_sumstats.py, ldsc.py) installed via conda into its own
# env (the bioconda ldsc recipe pins python<3, so it can't share the base env).
FROM rocker/r-ver:4.4.0

RUN install2.r --error \
        data.table \
        optparse

RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates bzip2 && \
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh && \
    /opt/conda/bin/conda create -y -n ldsc -c bioconda -c conda-forge ldsc=1.0.1 && \
    /opt/conda/bin/conda clean -afy && \
    rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/conda/envs/ldsc/bin:$PATH
