#!/usr/bin/env Rscript
# Normalized read-depth plot for one region from mosdepth binned depth
# (mosdepth_bins.sh). Depth is divided by the median of autosomal bins
# (depth > 0); the right axis shows 2 x ratio as inferred copy number.
#
# Usage:
#   Rscript plot_normalized_coverage.R --bed <prefix.regions.bed.gz> \
#     --region chr1:10000000-12000000 --out <out_prefix> \
#     [--buffer 0] [--smooth 25|auto] [--smooth-frac 0.5] [--mapq N] \
#     [--title "text"] [--highlight start-end]
#
# Outputs: <out_prefix>.pdf, .png, _bins.tsv

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
  library(dplyr)
  library(readr)
})

# ---- Command-line arguments ----
# Simple "--flag value" parsing; missing flags fall back to the default.
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else args[i + 1]
}

bed_file  <- get_arg("--bed")
region    <- get_arg("--region")
out       <- get_arg("--out")
smooth    <- get_arg("--smooth", "25")
smooth_frac <- as.numeric(get_arg("--smooth-frac", "0.5"))
mapq      <- get_arg("--mapq")
title     <- get_arg("--title", basename(out))
highlight <- get_arg("--highlight")
buffer    <- get_arg("--buffer", "0")

if (is.null(bed_file) || is.null(region) || is.null(out)) {
  stop("required: --bed <regions.bed.gz> --region chr:start-end --out <prefix>")
}

# ---- Region and plotting window ----
# "start-end" -> c(start, end); commas in numbers are allowed.
parse_range <- function(x) as.numeric(gsub(",", "", strsplit(x, "-")[[1]]))
m <- regmatches(region, regexec("^([^:]+):([0-9,]+-[0-9,]+)$", region))[[1]]
if (length(m) == 0) stop("--region must look like chr1:10000000-12000000")
region_chrom <- m[2]
region_range <- parse_range(m[3])

# Plotting window = event +/- buffer on each side (clipped at 0).
buffer_bp   <- as.numeric(buffer)
if (is.na(buffer_bp) || buffer_bp < 0) stop("--buffer must be a number of bp, e.g. 1000000")
plot_range  <- c(max(0, region_range[1] - buffer_bp), region_range[2] + buffer_bp)

# ---- Read mosdepth bins (chrom, start, end, mean depth) ----
bins <- read_tsv(bed_file, col_names = c("chrom", "start", "end", "depth"),
                 col_types = "cddd", progress = FALSE)

# Genome-wide baseline: autosomes only (sex chromosomes would skew it in XY
# samples); zero-depth bins are assembly gaps/unmappable, not real depth.
autosomal <- bins %>%
  filter(grepl("^(chr)?([1-9]|1[0-9]|2[0-2])$", chrom), depth > 0)
genome_median <- median(autosomal$depth)
# Bin size = most common bin width (the last bin of each chromosome is shorter).
bin_size <- as.numeric(names(which.max(table(bins$end - bins$start))))
message(sprintf("Autosomal median depth: %.2f x (%d bins of %g bp)",
                genome_median, nrow(autosomal), bin_size))

# ---- Bins in the plotting window, normalized to the baseline ----
# mid = bin centre (x position); ratio = normalized depth (1 = baseline).
reg <- bins %>%
  filter(chrom == region_chrom, end > plot_range[1], start < plot_range[2]) %>%
  arrange(start) %>%
  mutate(mid = (start + end) / 2, ratio = depth / genome_median)
if (nrow(reg) == 0) stop("no bins in ", region, " -- check chromosome naming (chr16 vs 16)")

# ---- Smoothing window (red line) ----
# Either a fixed number of bins, or 'auto' = smooth_frac x the event length in bins.
if (smooth == "auto") {
  if (is.na(smooth_frac) || smooth_frac <= 0) stop("--smooth-frac must be a positive number, e.g. 0.5")
  event_bins <- (region_range[2] - region_range[1]) / bin_size
  smooth_k <- max(3L, as.integer(round(event_bins * smooth_frac)))
  message(sprintf("--smooth auto: %.0f event bins x %g -> window %d bins",
                  event_bins, smooth_frac, smooth_k))
} else {
  smooth_k <- as.integer(smooth)
  if (is.na(smooth_k)) stop("--smooth must be a whole number of bins or 'auto'")
}
# runmed needs an odd window no larger than the number of bins
k <- min(smooth_k, nrow(reg))
if (k %% 2 == 0) k <- k - 1
reg$smooth <- if (k >= 3) as.numeric(runmed(reg$ratio, k, endrule = "median")) else reg$ratio

# Per-bin values behind the plot, for re-plotting or checking numbers.
write_tsv(reg, paste0(out, "_bins.tsv"))

# ---- Plot ----
# Caption records the bin size, MAPQ, baseline depth and smoothing window used.
caption <- sprintf(
  "Grey: mean depth per %g-bp bin (mosdepth%s) / autosomal median (%.1fx). Red: running median of %d bins%s.",
  bin_size, if (is.null(mapq)) "" else paste0(", MAPQ>=", mapq), genome_median, k,
  if (smooth == "auto") sprintf(" (auto: %g x event length)", smooth_frac) else "")

p <- ggplot(reg, aes(x = mid / 1e6, y = ratio))

# Optional shaded span (drawn first so points sit on top).
if (!is.null(highlight)) {
  hl <- parse_range(highlight)
  p <- p + annotate("rect", xmin = hl[1] / 1e6, xmax = hl[2] / 1e6,
                    ymin = -Inf, ymax = Inf, fill = "#1565C0", alpha = 0.08)
}

# Dotted lines at the event start and end.
p <- p + geom_vline(xintercept = region_range / 1e6, linetype = "dotted",
                    linewidth = 0.6, colour = "black")

# Dashed reference lines at CN 1, 2, 3, 4; grey points = bins; red line = smoothed.
# Right axis = 2 x normalized depth. The y-range comes from the data.
p <- p +
  geom_hline(yintercept = c(0.5, 1, 1.5, 2), linetype = "22",
             linewidth = 0.5, colour = "grey60") +
  geom_point(size = 0.5, colour = "grey55", alpha = 0.6) +
  geom_line(aes(y = smooth), colour = "#C62828", linewidth = 0.7) +
  scale_x_continuous(name = paste0(region_chrom, " position (Mb)")) +
  scale_y_continuous(name = "Normalized depth",
                     sec.axis = sec_axis(~ . * 2, name = "Inferred copy number")) +
  labs(title = title, caption = caption) +
  theme_prism(base_size = 13, base_fontface = "plain", base_line_size = 0.5) +
  theme(plot.caption = element_text(size = 9, colour = "grey30", hjust = 0))

# ---- Save ----
ggsave(paste0(out, ".pdf"), p, width = 11, height = 4.5)
ggsave(paste0(out, ".png"), p, width = 11, height = 4.5, dpi = 200)
message("Wrote ", out, ".pdf / .png / _bins.tsv")
