# quay.io/biocontainers/magma:1.10--h9f5acd7_0 started returning "unauthorized"
# on pull (confirmed independent of CX3 - reproduced from outside the cluster
# too, and other biocontainers images pull fine, so this is specific to that
# repository/tag, not a general outage). Turns out there is no bioconda
# package named "magma" for this tool at all (every bioconda/conda-forge
# "magma" hit is the unrelated linear-algebra library) - the recipe the
# biocontainer was built from appears to have been removed, likely due to
# that exact name collision. Fetch the real static binary directly from its
# official distributor (CTG lab, VU Amsterdam) instead.
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends wget unzip ca-certificates && \
    wget -q "https://vu.data.surfsara.nl/index.php/s/lxDgt2dNdNr6DYt/download" -O /tmp/magma.zip && \
    unzip -q /tmp/magma.zip -d /opt/magma && \
    chmod +x /opt/magma/magma && \
    rm /tmp/magma.zip && \
    rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/magma:$PATH
