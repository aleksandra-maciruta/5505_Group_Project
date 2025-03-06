#!/bin/bash

# Download raw reads with fetchngs -- should create samplesheet??
nextflow run nf-core/fetchngs \
    --input data/ids.csv \
    --outdir ./results/fetchngs \
    --max_cpus 32 --max_memory 128.GB \
    --download_method sratools \
    --nf_core_pipeline rnaseq \
    -profile docker \
    -w work/fetchngs \
