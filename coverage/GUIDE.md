# Normalized read depth / inferred copy number

Requires `mosdepth` and R (`../envs/rplot.yml`).

```bash
# 1. binned depth over the whole genome
#    bash mosdepth_bins.sh <bam|cram> <out_prefix> [bin=1000] [mapq=10] [threads=4] [ref.fa]
bash mosdepth_bins.sh sample.bam out/sample 1000 10 4

# 2. plot a region
Rscript plot_normalized_coverage.R --bed out/sample.regions.bed.gz \
  --region chr1:10000000-12000000 --buffer 1000000 --smooth auto \
  --out out/sample_chr1_10000000-12000000
```

| Option | |
|---|---|
| `--buffer` | flank in bp on each side (default 0) |
| `--smooth` | running-median window in bins (default 25), or `auto` = `--smooth-frac` × event length |
| `--smooth-frac` | fraction used by `--smooth auto` (default 0.5) |
| `--title` | plot title |
| `--highlight` | shade `start-end` |
| `--mapq` | caption label only |

Outputs: `.png`, `.pdf`, `_bins.tsv`.

Normalized depth = bin depth / median depth of autosomal bins (depth > 0). Copy number = 2 × normalized depth.

Mosdepth: Pedersen & Quinlan 2018, *Bioinformatics*, doi:10.1093/bioinformatics/btx699.
