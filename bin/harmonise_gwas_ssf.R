#!/usr/bin/env Rscript
# harmonise_gwas_ssf.R
# Convert a GWAS-SSF / GWAS Catalog harmonised .tsv(.gz) into the canonical
# harmonised-sumstats hub:
#   rsid chr pos ea oa eaf beta se p n z   (tab-separated, gzipped)
#
# Alleles/coords are aligned to the LD reference snp.info (joined by rsID):
#   - ea/oa  = snp.info A1/A2 (one common orientation across all traits)
#   - beta/eaf flipped to match ea = A1
#   - chr/pos = snp.info Chrom/PhysPos (carry the panel's build; recorded in meta.build)

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--input",    type = "character", help = "GWAS-SSF .tsv(.gz) sumstats"),
    make_option("--snp_info", type = "character", help = "LD reference snp.info"),
    make_option("--n",        type = "integer",   help = "Sample size fallback if file has no N"),
    make_option("--output",   type = "character", help = "Output harmonised .tsv.gz path")
)))

stopifnot(file.exists(opt$input), file.exists(opt$snp_info),
          !is.null(opt$n), !is.null(opt$output))

dt <- if (grepl("\\.gz$", opt$input, ignore.case = TRUE)) {
    fread(cmd = sprintf("zcat %s", shQuote(opt$input)))
} else fread(opt$input)

# GWAS-SSF permits either `rsid` or `variant_id` for the rsID column.
if (!"rsid" %in% names(dt)) {
    alt <- intersect(c("variant_id", "rs_id", "SNP", "snp"), names(dt))
    if (length(alt)) setnames(dt, alt[1], "rsid")
}

needed <- c("rsid", "effect_allele", "other_allele",
            "effect_allele_frequency", "beta", "standard_error", "p_value")
missing <- setdiff(needed, names(dt))
if (length(missing)) stop("Sumstats missing required columns: ",
                          paste(missing, collapse = ", "))

# ---- per-SNP N: effective N for case/control, else N_total, else fixed --n ----
n_case_col  <- grep("^n_case$",  names(dt), ignore.case = TRUE, value = TRUE)[1]
n_total_col <- grep("^n_total$", names(dt), ignore.case = TRUE, value = TRUE)[1]
n_col       <- grep("^n$",       names(dt), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(n_case_col) && !is.na(n_total_col)) {
    nc <- as.numeric(dt[[n_case_col]]); nt <- as.numeric(dt[[n_total_col]])
    n_eff <- 4 * nc * (nt - nc) / nt
    n_eff[!is.finite(n_eff) | n_eff <= 0] <- NA_real_
    dt[, n_pipeline := n_eff]
} else if (!is.na(n_total_col)) {
    dt[, n_pipeline := as.numeric(dt[[n_total_col]])]
} else if (!is.na(n_col)) {
    dt[, n_pipeline := as.numeric(dt[[n_col]])]
} else {
    dt[, n_pipeline := opt$n]
}

# ---- join to LD ref snp.info on rsID -> attach chr/pos + ref alleles ----
snpinfo <- fread(opt$snp_info, select = c("Chrom", "ID", "PhysPos", "A1", "A2"))
setnames(snpinfo, c("Chrom", "ID", "PhysPos"), c("chr", "rsid", "pos"))
setkey(snpinfo, rsid); setkey(dt, rsid)
m <- snpinfo[dt, nomatch = NULL]

# ---- orient to ea = ref A1 (flip beta + eaf when input effect allele is A2) ----
m[, flip := fcase(
    effect_allele == A1 & other_allele == A2, FALSE,
    effect_allele == A2 & other_allele == A1, TRUE
)]
m <- m[!is.na(flip)]                       # drop ambiguous / non-matching alleles
m[, beta_out := fifelse(flip, -beta, beta)]
m[, eaf_out  := fifelse(flip, 1 - effect_allele_frequency, effect_allele_frequency)]

out <- m[, .(rsid = rsid,
             chr  = chr,
             pos  = pos,
             ea   = A1,
             oa   = A2,
             eaf  = eaf_out,
             beta = beta_out,
             se   = standard_error,
             p    = p_value,
             n    = n_pipeline)]
out[, z := beta / se]

out <- out[!is.na(rsid) & nzchar(rsid) & rsid != "."]
out <- out[!is.na(beta) & !is.na(se) & se > 0 & !is.na(n) & n > 0]
setorder(out, chr, pos)

fwrite(out, opt$output, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[harmonise:gwas-ssf] wrote %d SNPs to %s\n", nrow(out), opt$output))
