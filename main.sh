#!/bin/bash

# Generate samplesheet
nextflow run nf-core/fetchngs \
--input data/ids.csv \
--outdir ./results/fetchngs \
--nf_core_pipeline rnaseq \
-profile docker \
-w work/fetchngs

