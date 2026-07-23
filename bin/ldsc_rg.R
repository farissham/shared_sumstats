#!/usr/bin/env Rscript
# ldsc_rg.R
# Cohort-level pairwise genetic correlation (rg) across all traits' munged
# sumstats, via the real bulik/ldsc software (ldsc.py --rg). Cross-trait, so it
# is not a per-trait fragment like ldsc_h2.R's output — one row per trait pair,
# written once for the whole cohort.
#
# ldsc.py --rg computes every pairwise rg across all inputs in a single call and
# prints a "Summary of Genetic Correlation Results" table at the end of its log;
# this script parses that table straight into a clean TSV.
#
# Fewer than two traits (rg is undefined) or an ldsc.py failure do not abort the
# run: an empty summary is written and the process exits 0, same as coloc.R's
# note_exit() pattern for single-locus failures.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--sumstats",    type = "character", help = "Comma-separated munged sumstats files (munge_sumstats.py output, .sumstats.gz), in the same order as --ids"),
    make_option("--ids",         type = "character", help = "Comma-separated trait ids, same order as --sumstats"),
    make_option("--ld_dir",      type = "character", help = "Directory of LDSC LD scores (eur_w_ld_chr)"),
    make_option("--out_summary", type = "character", help = "Output TSV: one row per trait pair")
)))

empty_summary <- function() {
    data.table(id1 = character(), id2 = character(), rg = numeric(), rg_se = numeric(),
               z = numeric(), p = numeric(),
               h2_obs1 = numeric(), h2_obs1_se = numeric(),
               h2_int1 = numeric(), h2_int1_se = numeric(),
               gcov_int = numeric(), gcov_int_se = numeric())
}

note_exit <- function(msg) {
    cat("[ldsc_rg] NOTE:", msg, "\n")
    fwrite(empty_summary(), opt$out_summary, sep = "\t")
    quit(status = 0)
}

sumstats <- strsplit(opt$sumstats, ",")[[1]]
ids      <- strsplit(opt$ids, ",")[[1]]

if (length(sumstats) != length(ids)) {
    note_exit(sprintf("mismatched --sumstats (%d) and --ids (%d) counts", length(sumstats), length(ids)))
}
if (length(sumstats) < 2) {
    note_exit(sprintf("only %d trait(s) — genetic correlation needs at least 2", length(sumstats)))
}

## ---- run ldsc.py --rg -----------------------------------------------------------
prefix   <- "cohort.rg"
ld_ref   <- paste0(opt$ld_dir, "/")
log_file <- paste0(prefix, ".log")

status <- system2("ldsc.py",
    c("--rg", paste(sumstats, collapse = ","),
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

## ---- parse "Summary of Genetic Correlation Results" table -----------------------
hdr_i <- grep("^Summary of Genetic Correlation Results", log_lines)
if (length(hdr_i) == 0) {
    note_exit("no 'Summary of Genetic Correlation Results' table found in ldsc.py log")
}

tbl_block <- log_lines[(hdr_i[1] + 1):length(log_lines)]
# ldsc.py follows the table with a blank line and then "Analysis finished at ..."
# / "Total time elapsed: ..." — stop at the first blank line so those don't get
# fed into fread() as bogus table rows.
first_blank <- which(!nzchar(trimws(tbl_block)))[1]
if (!is.na(first_blank)) tbl_block <- tbl_block[seq_len(first_blank - 1)]

tbl <- fread(paste(tbl_block, collapse = "\n"))

# p1/p2 columns are the munged sumstats file paths ldsc.py was given; map them
# back to trait ids via the --sumstats/--ids ordering (basename match, since
# ldsc.py may echo relative or absolute paths).
path_to_id <- setNames(ids, basename(sumstats))
tbl[, id1 := path_to_id[basename(p1)]]
tbl[, id2 := path_to_id[basename(p2)]]

out <- tbl[, .(id1, id2, rg, rg_se = se, z, p,
               h2_obs1 = h2_obs, h2_obs1_se = h2_obs_se,
               h2_int1 = h2_int, h2_int1_se = h2_int_se,
               gcov_int, gcov_int_se)]

fwrite(out, opt$out_summary, sep = "\t")

cat(sprintf("[ldsc_rg] %d trait(s), %d pair(s) written to %s\n",
            length(ids), nrow(out), opt$out_summary))
