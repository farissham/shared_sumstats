#!/usr/bin/env Rscript
# Derive independent genome-wide-significant loci from raw GWAS sumstats
# (CHR POS REF ALT ALT_AF beta chi_sq_stat p_value N), for use as SuSiE/coloc
# --loci input. Filters to a p-value threshold, greedily picks lead SNPs in
# significance order, expands each to a +/-window around the lead, then merges
# any windows that still overlap (two independent leads closer than 2*window).

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: derive_loci.R <gwas_sum_stats.tsv> <out loci.csv> [pval_threshold] [window_bp]")
infile  <- args[1]
outfile <- args[2]
pthresh <- if (length(args) >= 3) as.numeric(args[3]) else 5e-8
window  <- if (length(args) >= 4) as.numeric(args[4]) else 500000

# Pre-filter via awk before it ever hits R's memory - this file is whole-genome.
cat(sprintf("Filtering %s to p_value < %g ...\n", infile, pthresh))
dt <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || $8 < %g' %s", pthresh, shQuote(infile)))
cat(sprintf("%d SNPs pass threshold.\n", nrow(dt)))
if (nrow(dt) == 0) stop("No SNPs pass this threshold - loosen pval_threshold and rerun.")

dt[, CHR := sub("^chr", "", CHR)]
setorder(dt, p_value)

# Greedy lead-SNP selection: take the most significant remaining SNP, drop
# every other SNP within +/-window of it, repeat.
remaining <- copy(dt)
leads <- data.table(chr = character(), lead_pos = numeric(), lead_p = numeric())
while (nrow(remaining) > 0) {
    lead <- remaining[1]
    leads <- rbind(leads, data.table(chr = lead$CHR, lead_pos = lead$POS, lead_p = lead$p_value))
    remaining <- remaining[!(CHR == lead$CHR & POS >= lead$POS - window & POS <= lead$POS + window)]
}
cat(sprintf("%d independent lead SNP(s) after clumping.\n", nrow(leads)))

loci <- leads[, .(chr = chr, start = pmax(1, lead_pos - window), end = lead_pos + window, lead_pos, lead_p)]
setorder(loci, chr, start)

# Merge any windows that still overlap (two leads 500kb-1Mb apart on the same chr).
merged <- data.table(chr = character(), start = numeric(), end = numeric())
i <- 1
while (i <= nrow(loci)) {
    cur_chr <- loci$chr[i]; cur_start <- loci$start[i]; cur_end <- loci$end[i]
    j <- i + 1
    while (j <= nrow(loci) && loci$chr[j] == cur_chr && loci$start[j] <= cur_end) {
        cur_end <- max(cur_end, loci$end[j])
        j <- j + 1
    }
    if (j > i + 1) cat(sprintf("Merged %d overlapping windows on chr%s into %d-%d\n", j - i, cur_chr, cur_start, cur_end))
    merged <- rbind(merged, data.table(chr = cur_chr, start = cur_start, end = cur_end))
    i <- j
}

fwrite(merged, outfile)
cat(sprintf("Wrote %d loci -> %s\n", nrow(merged), outfile))
print(merged)
