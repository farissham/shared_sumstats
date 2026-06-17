#!/bin/bash
#PBS -N susie_pipeline
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=24:00:00
#PBS -q v1_small24
#PBS -V

module load Nextflow/25.10.2

export TOWER_ACCESS_TOKEN=$(cat /rds/general/user/<username>/home/.tower_token)

PIPELINE_DIR=/rds/general/project/<project>/live/<user>/susie_pipeline

cd ${PIPELINE_DIR}/scripts/

nextflow run ${PIPELINE_DIR}/main.nf \
    --sumstats /path/to/sumstats.tsv.gz \
    --ref      /path/to/g1000_eur \
    --loci     ${PIPELINE_DIR}/scripts/loci.csv \
    --outdir   /path/to/output \
    --L        10 \
    --max_iter 1000 \
    -c         ${PIPELINE_DIR}/scripts/hpc.config \
    -with-tower
