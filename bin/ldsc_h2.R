#!/usr/bin/env Rscript
# ldsc_h2.R
# Per-trait SNP-heritability via LD Score Regression, using the real bulik/ldsc
# software (ldsc.py --h2) rather than a from-scratch reimplementation, so results
# are directly comparable to published LDSC numbers (e.g. Jurgens et al. 2024).
#
# This script does not do the regression itself: it shells out to ldsc.py on the
# already-munged sumstats (see LDSC_MUNGE / munge_sumstats.py) and reshapes the
# text log it writes into a clean one-row TSV summary.
#
# ldsc.py failures (e.g. too few SNPs after merging with the LD reference) do not
# abort the run: a placeholder row is written and the process exits 0, same as
# coloc.R's note_exit() pattern for single-locus failures.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--sumstats",    type = "character", help = "Munged sumstats (munge_sumstats.py output, .sumstats.gz)"),
    make_option("--ld_dir",      type = "character", help = "Directory of LDSC LD scores (eur_w_ld_chr: per-chromosome .l2.ldscore.gz + .l2.M_5_50)"),
    make_option("--id",          type = "character", help = "Trait id, written to the output and used as the ldsc.py --out prefix"),
    make_option("--out_summary", type = "character", help = "Output TSV: one row of h2 + diagnostics")
)))

# Write a single placeholder row and exit 0.
note_exit <- function(msg) {
    cat("[ldsc_h2]", opt$id, "NOTE:", msg, "\n")
    fwrite(data.table(id = opt$id,
                       h2 = NA_real_, h2_se = NA_real_,
                       lambda_gc = NA_real_, mean_chi2 = NA_real_,
                       intercept = NA_real_, intercept_se = NA_real_,
                       ratio = NA_real_, ratio_se = NA_real_,
                       note = msg),
           opt$out_summary, sep = "\t")
    quit(status = 0)
}

# Pull a "value (se)" pair off a log line, e.g. "Total Observed scale h2: 0.1234 (0.0231)".
val_se <- function(log_lines, pattern) {
    line <- grep(pattern, log_lines, value = TRUE)
    if (length(line) == 0) return(c(NA_real_, NA_real_))
    m <- regmatches(line[1], regexec("(-?[0-9.eE+-]+)\\s*\\(([0-9.eE+-]+)\\)", line[1]))[[1]]
    if (length(m) < 3) return(c(NA_real_, NA_real_))
    as.numeric(m[2:3])
}

# Pull a single value off a log line with no SE, e.g. "Lambda GC: 1.0512".
val_only <- function(log_lines, pattern) {
    line <- grep(pattern, log_lines, value = TRUE)
    if (length(line) == 0) return(NA_real_)
    suppressWarnings(as.numeric(sub(".*:\\s*", "", line[1])))
}

## ---- run ldsc.py --h2 ---------------------------------------------------------
prefix   <- paste0(opt$id, ".h2")
ld_ref   <- paste0(opt$ld_dir, "/")
log_file <- paste0(prefix, ".log")

status <- system2("ldsc.py",
    c("--h2", opt$sumstats,
      "--ref-ld-chr", ld_ref,
      "--w-ld-chr",   ld_ref,
      "--out", prefix),
    stdout = TRUE, stderr = TRUE)

if (!file.exists(log_file)) {
    note_exit(paste("ldsc.py produced no log file; output:", paste(status, collapse = " | ")))
}

log_lines <- readLines(log_file)

if (length(grep("ERROR", log_lines)) > 0) {
    note_exit(paste("ldsc.py reported an error:", paste(grep("ERROR", log_lines, value = TRUE), collapse = "; ")))
}

## ---- parse the log -------------------------------------------------------------
h2        <- val_se(log_lines, "^Total Observed scale h2:")
lambda_gc <- val_only(log_lines, "^Lambda GC:")
mean_chi2 <- val_only(log_lines, "^Mean Chi\\^2:")
intercept <- val_se(log_lines, "^Intercept:")

# "Ratio" is undefined (and ldsc.py prints an explanatory sentence instead of a
# number) when mean chi^2 <= 1; treat that as NA with a note rather than a parse error.
ratio_line <- grep("^Ratio", log_lines, value = TRUE)
if (length(ratio_line) == 0) {
    ratio <- c(NA_real_, NA_real_)
    ratio_note <- ""
} else if (grepl("usually indicates", ratio_line[1])) {
    ratio <- c(NA_real_, NA_real_)
    ratio_note <- "ratio undefined: mean chi^2 <= 1"
} else {
    ratio <- val_se(log_lines, "^Ratio:")
    ratio_note <- ""
}

## ---- write output ---------------------------------------------------------------
fwrite(data.table(id = opt$id,
                   h2 = h2[1], h2_se = h2[2],
                   lambda_gc = lambda_gc, mean_chi2 = mean_chi2,
                   intercept = intercept[1], intercept_se = intercept[2],
                   ratio = ratio[1], ratio_se = ratio[2],
                   note = ratio_note),
       opt$out_summary, sep = "\t")

cat(sprintf("[ldsc_h2] %s  h2=%.4f (%.4f)  intercept=%.4f (%.4f)\n",
            opt$id, h2[1], h2[2], intercept[1], intercept[2]))
