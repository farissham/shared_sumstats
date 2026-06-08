#!/usr/bin/env Rscript
# verify_harmonised.R
# Independent post-hoc checks on a harmonised-sumstats hub (<id>.harmonised.tsv.gz).
# Runs four verification layers and prints a PASS/FAIL report:
#   1. structural invariants  (columns, allele chars, value ranges, rsID)
#   2. orientation            (ea==A1 and oa==A2 against the LD reference snp.info)
#   3. internal consistency   (z == beta/se ; p ~= 2*pnorm(-|z|))
#   4. eaf sanity             (harmonised eaf correlates with snp.info A1Freq)
# Exit code is non-zero if any hard check FAILs.

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--harmonised", type = "character", help = "Harmonised hub .tsv(.gz)"),
    make_option("--snp_info",   type = "character", help = "LD reference snp.info"),
    make_option("--eaf_r_min",  type = "double", default = 0.7,
                help = "Min Pearson r between harmonised eaf and A1Freq (default 0.7)")
)))
stopifnot(file.exists(opt$harmonised), file.exists(opt$snp_info))

hub <- if (grepl("\\.gz$", opt$harmonised, ignore.case = TRUE))
    fread(cmd = sprintf("zcat %s", shQuote(opt$harmonised))) else fread(opt$harmonised)
snp <- fread(opt$snp_info, select = c("ID", "A1", "A2", "A1Freq"))
setnames(snp, "ID", "rsid")

fails <- character(0)
pass  <- function(tag, ok, msg) {
    cat(sprintf("[%s] %s  %s\n", if (ok) "PASS" else "FAIL", tag, msg))
    if (!ok) fails[[length(fails) + 1L]] <<- tag
}

cat(sprintf("Harmonised rows: %d\n\n", nrow(hub)))

## ---- 1. structural invariants ----
expect <- c("rsid","chr","pos","ea","oa","eaf","beta","se","p","n","z")
pass("1.header", identical(names(hub), expect),
     sprintf("columns = %s", paste(names(hub), collapse = ",")))
acgt <- function(x) all(x %in% c("A","C","G","T"))
pass("1.alleles", acgt(hub$ea) && acgt(hub$oa) && all(hub$ea != hub$oa),
     "ea/oa in {A,C,G,T} and ea != oa")
pass("1.ranges",
     all(hub$eaf >= 0 & hub$eaf <= 1, na.rm = TRUE) &&
     all(hub$p > 0 & hub$p <= 1, na.rm = TRUE) &&
     all(hub$se > 0, na.rm = TRUE) && all(hub$n > 0, na.rm = TRUE),
     "eaf in [0,1], p in (0,1], se>0, n>0")
ndup <- sum(duplicated(hub$rsid))
pass("1.rsid", all(nzchar(hub$rsid)) && !any(is.na(hub$rsid)) && ndup == 0,
     sprintf("non-missing, %d duplicate rsIDs", ndup))

## ---- 2. orientation vs snp.info ----
m <- merge(hub, snp, by = "rsid")
if (nrow(m) == 0) {
    pass("2.orientation", FALSE, "no overlap with snp.info to check")
} else {
    bad <- m[!(ea == A1 & oa == A2)]
    pass("2.orientation", nrow(bad) == 0,
         sprintf("%d/%d rows checked; %d with ea!=A1 or oa!=A2", nrow(m), nrow(hub), nrow(bad)))
}

## ---- 3. internal consistency ----
z_recalc <- hub$beta / hub$se
pass("3.z=beta/se", all(abs(z_recalc - hub$z) < 1e-6 * (1 + abs(hub$z)), na.rm = TRUE),
     sprintf("max |z - beta/se| = %.2e", max(abs(z_recalc - hub$z), na.rm = TRUE)))
p_recalc <- 2 * pnorm(-abs(hub$z))
ok_p <- is.finite(p_recalc) & hub$p > 0 & p_recalc > 0
agree <- mean(abs(log10(hub$p[ok_p]) - log10(p_recalc[ok_p])) < 0.5)
pass("3.p~2pnorm(-|z|)", agree > 0.95,
     sprintf("%.1f%% of SNPs agree within 0.5 on log10(p)", 100 * agree))

## ---- 4. eaf vs A1Freq ----
if (nrow(m) >= 3) {
    r <- suppressWarnings(cor(m$eaf, m$A1Freq, use = "complete.obs"))
    pass("4.eaf~A1Freq", is.finite(r) && r >= opt$eaf_r_min,
         sprintf("Pearson r = %.3f (min %.2f)", r, opt$eaf_r_min))
} else {
    pass("4.eaf~A1Freq", FALSE, "too few overlapping SNPs")
}

cat("\n")
if (length(fails) == 0) {
    cat("RESULT: ALL CHECKS PASSED\n"); quit(status = 0)
} else {
    cat(sprintf("RESULT: %d CHECK(S) FAILED: %s\n", length(fails), paste(fails, collapse = ", ")))
    quit(status = 1)
}
