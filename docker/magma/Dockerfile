# quay.io/biocontainers/magma:1.10--h9f5acd7_0 started returning "unauthorized"
# on pull (confirmed independent of CX3 - reproduced from outside the cluster
# too, and other biocontainers images pull fine, so this is specific to that
# repository/tag, not a general outage). Turns out there is no bioconda
# package named "magma" for this tool at all (every bioconda/conda-forge
# "magma" hit is the unrelated linear-algebra library) - the recipe the
# biocontainer was built from appears to have been removed, likely due to
# that exact name collision.
#
# Fetching the real static binary at build time (from CTG's official SURFsara
# distribution) failed repeatedly in GitHub Actions CI specifically (wget
# exit 4, network failure) even with retries + forced IPv4, while the
# identical download succeeds immediately from outside CI - most likely that
# host blocks/rate-limits known CI IP ranges. Sidestep the whole class of
# problem: the zip is small (~3.2MB) and vendored directly in this directory
# (magma_v1.10.zip, downloaded and verified locally: real application/zip,
# unzips to a flat magma/README/CHANGELOG/manual_v1.10.pdf layout, magma
# binary is a valid statically-linked x86-64 ELF executable), so the build
# never needs to reach that host at all.
FROM ubuntu:22.04

COPY magma_v1.10.zip /tmp/magma.zip

RUN apt-get update && apt-get install -y --no-install-recommends unzip ca-certificates && \
    unzip -q /tmp/magma.zip -d /opt/magma && \
    chmod +x /opt/magma/magma && \
    rm /tmp/magma.zip && \
    rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/magma:$PATH
