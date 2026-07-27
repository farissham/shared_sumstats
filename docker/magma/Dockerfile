# quay.io/biocontainers/magma:1.10--h9f5acd7_0 started returning "unauthorized"
# on pull (confirmed independent of CX3 - reproduced from outside the cluster
# too, and other biocontainers images pull fine, so this is specific to that
# repository/tag, not a general outage). Same fix as susie/ldsc: build our own
# from the same bioconda package instead of depending on that biocontainer.
FROM rocker/r-ver:4.4.0

RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates bzip2 && \
    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh && \
    /opt/conda/bin/conda create -y -n magma -c bioconda -c conda-forge magma=1.10 && \
    /opt/conda/bin/conda clean -afy && \
    rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/conda/envs/magma/bin:$PATH
