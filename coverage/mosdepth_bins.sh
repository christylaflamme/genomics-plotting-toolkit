#!/usr/bin/env bash
# Binned read depth across the whole genome with mosdepth -- input for
# plot_normalized_coverage.R. Whole genome (not just the region of interest)
# because the plot normalizes to the genome-wide autosomal median.
#
# Run on a compute node, not a login node (reads the entire BAM; ~10-30 min for 30x WGS).
#
# Usage:
#   bash mosdepth_bins.sh <bam|cram> <out_prefix> [bin_size=1000] [min_mapq=10] [threads=4] [reference.fa]
#
# Output: <out_prefix>.regions.bed.gz  (chrom, start, end, mean depth per bin)
# A CRAM needs the reference FASTA it was compressed against (6th argument).

set -euo pipefail

if [[ $# -lt 2 ]]; then
  sed -n '2,13p' "$0"; exit 1
fi

# Positional arguments, with defaults for the optional ones
bam=$1
out=$2
bin=${3:-1000}
mapq=${4:-10}
threads=${5:-4}
ref=${6:-}

mkdir -p "$(dirname "$out")"

# --no-per-base: skip the large per-base output (only bins are needed)
# --by:          bin size in bp; --mapq: minimum mapping quality
# --fasta:       only passed when a reference is given (needed for CRAM)
mosdepth --threads "$threads" --no-per-base --by "$bin" --mapq "$mapq" \
  ${ref:+--fasta "$ref"} "$out" "$bam"

echo "Wrote ${out}.regions.bed.gz (bin=${bin} bp, MAPQ>=${mapq})"
