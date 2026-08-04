# Python3 + scipy for concordance_analysis.py / concordance.py
# (belowlab/Concordance-Analysis, vendored in bin/).
FROM python:3.11-slim

# procps provides `ps`, which Nextflow requires inside every container to
# collect task resource metrics - python:3.11-slim strips it out by default,
# unlike this pipeline's other (larger) base images, which happen to include
# it already. Without it, every task fails immediately (~10s) regardless of
# retries, since it's a structural container gap, not a resource issue.
RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir scipy
