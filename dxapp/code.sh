#!/bin/bash
set -euo pipefail

echo "Starting Nextflow RD Genomic Pipeline by R Newman"

dx-download-all-inputs

# Use explicit DNAnexus-provided inputs
SAMPLESHEET="$samplesheet"
REF="$reference"

echo "Samplesheet: $SAMPLESHEET"
echo "Reference: $REF"

# -----------------------------
# Output to a specific folder
# -----------------------------
OUTDIR="./results/RN_results"

mkdir -p "$OUTDIR"

# -----------------------------
# Run Nextflow
# -----------------------------
nextflow run main.nf \
  --samplesheet "$SAMPLESHEET" \
  --genome_file "$REF" \
  --aligner "$aligner" \
  --variant_caller "$variant_caller" \
  --fastp "$fastp" \
  --fastqc "$fastqc" \
  --outdir "$OUTDIR" \
  -profile docker \
  -resume

# -----------------------------
# Upload outputs
# -----------------------------
dx-upload-all-outputs

echo "Done"