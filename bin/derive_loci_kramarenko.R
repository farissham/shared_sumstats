#!/usr/bin/env Rscript
# Derive independent genome-wide-significant loci from the Kramarenko CC_MTAG
# file (rsID CHR BP SNP A1 A2 EAFREQ BETA SE Neff Z _samples P), for use as
# SuSiE/coloc --loci input. Same greedy lead-SNP clumping logic as
# bin/derive_loci.R, adapted to this file's own column layout/positions
# (CHR=$2, BP=$3, P=$13) rather than the hail_step3-specific one that script
# hardcodes.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: derive_loci_kramarenko.R <CC_MTAG file> <out loci.csv> [pval_threshold] [window_bp]")
infile  <- args[1]
outfile <- args[2]
pthresh <- if (length(args) >= 3) as.numeric(args[3]) else 5e-8
window  <- if (length(args) >= 4) as.numeric(args[4]) else 500000

cat(sprintf("Filtering %s to P < %g ...\n", infile, pthresh))
dt <- fread(cmd = sprintf("awk -F'\\t' 'NR==1 || $13 < %g' %s", pthresh, shQuote(infile)))
cat(sprintf("%d SNPs pass threshold.\n", nrow(dt)))
if (nrow(dt) == 0) stop("No SNPs pass this threshold - loosen pval_threshold and rerun.")

setnames(dt, c("CHR", "BP", "P"), c("CHR", "POS", "p_value"), skip_absent = TRUE)
dt[, CHR := sub("^chr", "", as.character(CHR))]
setorder(dt, p_value)

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
merged <- loci[1]
if (nrow(loci) > 1) {
    for (i in 2:nrow(loci)) {
        last <- nrow(merged)
        if (loci$chr[i] == merged$chr[last] && loci$start[i] <= merged$end[last]) {
            merged$end[last] <- max(merged$end[last], loci$end[i])
            merged$lead_p[last] <- min(merged$lead_p[last], loci$lead_p[i])
        } else {
            merged <- rbind(merged, loci[i])
        }
    }
}
cat(sprintf("%d final locus/loci after merging overlaps.\n", nrow(merged)))

fwrite(merged[, .(chr, start, end)], outfile)
cat(sprintf("Written to %s\n", outfile))
