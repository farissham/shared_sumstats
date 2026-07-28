#!/usr/bin/env Rscript
# coloc.R
# Per-locus Bayesian colocalization (coloc.abf) between two GWAS/QTL datasets.
#
# Column names in each input file are configurable via --hub_*_col and --ref_*_col
# options, defaulting to hub format (rsid chr pos ea oa beta se n eaf).
#
# For the DCM Jurgens sumstats (trait 1) use:
#   --hub_snp_col rsID --hub_chr_col CHR --hub_pos_col POS \
#   --hub_ea_col EA --hub_oa_col NEA --hub_beta_col BETA \
#   --hub_se_col SE --hub_n_col N --hub_eaf_col EAFREQ
#
# Allele harmonisation aligns trait 2 onto trait 1's coding before coloc.abf().
# Palindromic SNPs (A/T, C/G) are dropped because strand is unknowable.
# Single-locus failures write placeholder outputs and exit 0.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(coloc)
})

opt <- parse_args(OptionParser(option_list = list(

    # --- input files ---
    make_option("--hub",  type = "character", help = "Trait 1 sumstats .tsv(.gz)"),
    make_option("--ref",  type = "character", help = "Trait 2 sumstats .tsv(.gz)"),

    # --- locus ---
    make_option("--chr",   type = "integer", help = "Chromosome"),
    make_option("--start", type = "integer", help = "Locus start (bp)"),
    make_option("--end",   type = "integer", help = "Locus end (bp)"),

    # --- trait 1 column names ---
    make_option("--hub_snp_col",  type = "character", default = "rsid",  help = "SNP ID column in hub"),
    make_option("--hub_chr_col",  type = "character", default = "chr",   help = "Chromosome column in hub"),
    make_option("--hub_pos_col",  type = "character", default = "pos",   help = "Position column in hub"),
    make_option("--hub_ea_col",   type = "character", default = "ea",    help = "Effect allele column in hub"),
    make_option("--hub_oa_col",   type = "character", default = "oa",    help = "Other allele column in hub"),
    make_option("--hub_beta_col", type = "character", default = "beta",  help = "Beta column in hub"),
    make_option("--hub_se_col",   type = "character", default = "se",    help = "SE column in hub"),
    make_option("--hub_n_col",    type = "character", default = "n",     help = "Sample size column in hub"),
    make_option("--hub_eaf_col",  type = "character", default = "eaf",   help = "EAF column in hub (optional)"),

    # --- trait 2 column names ---
    make_option("--ref_snp_col",  type = "character", default = "rsid",  help = "SNP ID column in ref"),
    make_option("--ref_chr_col",  type = "character", default = "chr",   help = "Chromosome column in ref"),
    make_option("--ref_pos_col",  type = "character", default = "pos",   help = "Position column in ref"),
    make_option("--ref_ea_col",   type = "character", default = "ea",    help = "Effect allele column in ref"),
    make_option("--ref_oa_col",   type = "character", default = "oa",    help = "Other allele column in ref"),
    make_option("--ref_beta_col",   type = "character", default = "beta",  help = "Beta column in ref"),
    make_option("--ref_se_col",     type = "character", default = "se",    help = "SE column in ref"),
    make_option("--ref_n_col",      type = "character", default = "n",     help = "Sample size column in ref"),
    make_option("--ref_eaf_col",    type = "character", default = "eaf",   help = "EAF column in ref (optional)"),
    make_option("--ref_zscore_col", type = "character", default = NULL,
                help = "Z-score column in ref (eQTLGen). When set, approx beta=z/sqrt(N) and se=1/sqrt(N) — no separate beta/se needed."),

    # --- dataset types and case proportions ---
    make_option("--type1", type = "character", default = "quant",      help = "Trait 1 type: quant or cc"),
    make_option("--type2", type = "character", default = "quant",      help = "Trait 2 type: quant or cc"),
    make_option("--s1",    type = "double",    default = NA_real_,     help = "Proportion of cases in trait 1 (cc only)"),
    make_option("--s2",    type = "double",    default = NA_real_,     help = "Proportion of cases in trait 2 (cc only)"),
    make_option("--n1",    type = "integer",   default = NA_integer_,  help = "Sample size fallback for trait 1"),
    make_option("--n2",    type = "integer",   default = NA_integer_,  help = "Sample size fallback for trait 2"),

    # --- coloc priors ---
    make_option("--p1",       type = "double",  default = 1e-4, help = "Prior: P(SNP assoc with trait 1)"),
    make_option("--p2",       type = "double",  default = 1e-4, help = "Prior: P(SNP assoc with trait 2)"),
    make_option("--p12",      type = "double",  default = 1e-5, help = "Prior: P(SNP assoc with both)"),
    make_option("--min_snps", type = "integer", default = 50,   help = "Min overlapping SNPs per locus"),

    # --- outputs ---
    make_option("--out_summary", type = "character", help = "Locus-level coloc summary TSV"),
    make_option("--out_snps",    type = "character", help = "Per-SNP SNP.PP.H4 TSV"),
    make_option("--gene_id",     type = "character", default = NULL,
                help = "Gene identifier written to summary (e.g. ref file basename)")
)))

comp <- function(x) chartr("ACGTacgt", "TGCAtgca", x)

read_sumstats <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        fread(cmd = sprintf("gzip -dc %s", shQuote(path)))
    } else {
        fread(path)
    }
}

# Rename columns to internal standard names after reading.
standardise_cols <- function(dt, snp, chr, pos, ea, oa, beta, se, n, eaf) {
    mapping <- c(rsid = snp, chr = chr, pos = pos, ea = ea, oa = oa,
                 beta = beta, se = se, n = n, eaf = eaf)
    for (std in names(mapping)) {
        orig <- mapping[[std]]
        if (orig %in% names(dt) && orig != std) setnames(dt, orig, std)
    }
    dt
}

# Write empty placeholder outputs and exit 0.
note_exit <- function(msg) {
    cat("[coloc] NOTE:", msg, "\n")
    fwrite(data.table(gene_id = if (!is.null(opt$gene_id)) opt$gene_id else "",
                      chr = opt$chr, start = opt$start, end = opt$end,
                      nsnps = 0L,
                      PP.H0.abf = NA_real_, PP.H1.abf = NA_real_,
                      PP.H2.abf = NA_real_, PP.H3.abf = NA_real_, PP.H4.abf = NA_real_,
                      note = msg),
           opt$out_summary, sep = "\t")
    fwrite(data.table(snp = character(), SNP.PP.H4 = numeric()), opt$out_snps, sep = "\t")
    quit(status = 0)
}

## ---- read and standardise column names --------------------------------------
hub <- read_sumstats(opt$hub)
ref <- read_sumstats(opt$ref)

hub <- standardise_cols(hub,
    opt$hub_snp_col, opt$hub_chr_col, opt$hub_pos_col,
    opt$hub_ea_col,  opt$hub_oa_col,  opt$hub_beta_col,
    opt$hub_se_col,  opt$hub_n_col,   opt$hub_eaf_col)

ref <- standardise_cols(ref,
    opt$ref_snp_col, opt$ref_chr_col, opt$ref_pos_col,
    opt$ref_ea_col,  opt$ref_oa_col,  opt$ref_beta_col,
    opt$ref_se_col,  opt$ref_n_col,   opt$ref_eaf_col)

# Z-score mode for ref (e.g. eQTLGen): approximate beta = z/sqrt(N), se = 1/sqrt(N).
# Sign of beta is preserved so allele harmonisation can still flip direction.
if (!is.null(opt$ref_zscore_col)) {
    if (!opt$ref_zscore_col %in% names(ref))
        stop("ref missing z-score column: ", opt$ref_zscore_col)
    if (!"n" %in% names(ref) && is.na(opt$n2))
        stop("ref z-score mode requires n column (--ref_n_col) or --n2 fallback")
    ref[, zscore_tmp := as.numeric(get(opt$ref_zscore_col))]
    ref[, n_tmp      := if ("n" %in% names(ref)) as.numeric(n) else as.numeric(opt$n2)]
    ref[, se   := 1 / sqrt(n_tmp)]
    ref[, beta := zscore_tmp * se]
    ref[, c("zscore_tmp", "n_tmp") := NULL]
    cat("[coloc] ref z-score mode: beta and se approximated from", opt$ref_zscore_col, "\n")
}

for (col in c("rsid", "chr", "pos", "ea", "oa", "beta", "se")) {
    if (!col %in% names(hub)) stop("hub missing column after rename: ", col)
    if (!col %in% names(ref)) stop("ref missing column after rename: ", col)
}

## ---- subset to locus ---------------------------------------------------------
hub <- hub[chr == opt$chr & pos >= opt$start & pos <= opt$end]
ref <- ref[chr == opt$chr & pos >= opt$start & pos <= opt$end]

if (nrow(hub) == 0) note_exit("no hub SNPs in this locus")
if (nrow(ref) == 0) note_exit("no ref SNPs in this locus")

## ---- deduplicate by rsid (keep highest |z| per SNP) -------------------------
dedup_rsid <- function(dt, label) {
    if (!anyDuplicated(dt$rsid)) return(dt)
    dt[, absz_tmp := abs(as.numeric(beta) / as.numeric(se))]
    dt <- dt[order(-absz_tmp)][!duplicated(rsid)]
    dt[, absz_tmp := NULL]
    cat(sprintf("[coloc] %s: removed duplicate rsIDs, %d unique SNPs remain\n", label, nrow(dt)))
    dt
}
hub <- dedup_rsid(hub, "hub")
ref <- dedup_rsid(ref, "ref")

hub[, c("ea", "oa") := .(toupper(ea), toupper(oa))]
ref[, c("ea", "oa") := .(toupper(ea), toupper(oa))]

## ---- join on rsid (primary) or chr:pos (fallback) ---------------------------
m <- merge(hub, ref, by = "rsid", suffixes = c(".1", ".2"))

if (nrow(m) < opt$min_snps) {
    hub[, chrpos := paste0(chr, ":", pos)]
    ref[, chrpos := paste0(chr, ":", pos)]
    m2 <- merge(hub, ref, by = "chrpos", suffixes = c(".1", ".2"))
    if (nrow(m2) > nrow(m)) m <- m2
}

if (nrow(m) < opt$min_snps) {
    note_exit(sprintf("only %d overlapping SNPs after join (< min_snps %d)", nrow(m), opt$min_snps))
}

## ---- allele harmonisation ---------------------------------------------------
m[, action := fcase(
    ea.2 == ea.1 & oa.2 == oa.1,             "keep",
    ea.2 == oa.1 & oa.2 == ea.1,             "flip",
    comp(ea.2) == ea.1 & comp(oa.2) == oa.1, "keep",
    comp(ea.2) == oa.1 & comp(oa.2) == ea.1, "flip",
    default = "drop")]
m[ea.1 == comp(oa.1), action := "drop"]  # palindrome
m <- m[action != "drop"]
m[, beta2_aln := fifelse(action == "flip", -as.numeric(beta.2), as.numeric(beta.2))]

cat(sprintf("[coloc] chr%d:%d-%d  keep=%d  flip=%d\n",
            opt$chr, opt$start, opt$end,
            sum(m$action == "keep"), sum(m$action == "flip")))

if (nrow(m) < opt$min_snps) {
    note_exit(sprintf("only %d SNPs after allele harmonisation (< min_snps %d)", nrow(m), opt$min_snps))
}

## ---- sample sizes ------------------------------------------------------------
# When only one dataset has an n column, data.table merge leaves it unsuffixed.
# When both have n, they become n.1 and n.2. Handle both cases.
n1 <- suppressWarnings({
    col <- if ("n.1" %in% names(m)) "n.1" else if ("n" %in% names(m)) "n" else NA_character_
    if (!is.na(col)) max(as.numeric(m[[col]]), na.rm = TRUE) else NA_real_
})
n2 <- suppressWarnings({
    col <- if ("n.2" %in% names(m)) "n.2" else NA_character_
    if (!is.na(col)) max(as.numeric(m[[col]]), na.rm = TRUE) else NA_real_
})
if (!is.finite(n1) || n1 <= 0) {
    if (!is.na(opt$n1) && opt$n1 > 0) n1 <- opt$n1 else
        stop("trait 1: no valid N — check --hub_n_col or pass --n1")
}
if (!is.finite(n2) || n2 <= 0) {
    if (!is.na(opt$n2) && opt$n2 > 0) n2 <- opt$n2 else
        stop("trait 2: no valid N — check --ref_n_col or pass --n2")
}

## ---- coloc dataset objects --------------------------------------------------
snp_ids <- if ("rsid" %in% names(m)) m$rsid else m$chrpos

## eaf.1/eaf.2 only exist if BOTH hub and ref have an "eaf" column (merge()
## only suffixes colliding names). If only one side has it - e.g. eQTLGen's
## ref never has an eaf-equivalent column at all - it survives the merge
## unsuffixed as plain "eaf", which "eaf.1"/"eaf.2" would never match. Same
## ambiguity already handled correctly for n.1/n above; apply the same
## fallback here instead of silently losing MAF (and therefore coloc.abf()'s
## only way to estimate sdY for type="quant") whenever no collision occurs.
## Resolved against the pre-merge hub/ref tables (still in scope), not
## guessed from the merged name alone, since an unsuffixed "eaf" in m could
## belong to either side depending on which one actually had the column.
eaf1_col <- if ("eaf.1" %in% names(m)) "eaf.1" else if ("eaf" %in% names(hub)) "eaf" else NA_character_
eaf2_col <- if ("eaf.2" %in% names(m)) "eaf.2" else if ("eaf" %in% names(ref)) "eaf" else NA_character_

ds1 <- list(snp     = snp_ids,
            beta    = as.numeric(m$beta.1),
            varbeta = as.numeric(m$se.1)^2,
            type    = opt$type1,
            N       = n1)
if (!is.na(eaf1_col) && !all(is.na(m[[eaf1_col]])))
    ds1$MAF <- pmin(as.numeric(m[[eaf1_col]]), 1 - as.numeric(m[[eaf1_col]]))
if (opt$type1 == "cc" && !is.na(opt$s1)) ds1$s <- opt$s1

ds2 <- list(snp     = snp_ids,
            beta    = m$beta2_aln,
            varbeta = as.numeric(m$se.2)^2,
            type    = opt$type2,
            N       = n2)
if (!is.na(eaf2_col) && !all(is.na(m[[eaf2_col]])))
    ds2$MAF <- pmin(as.numeric(m[[eaf2_col]]), 1 - as.numeric(m[[eaf2_col]]))
if (opt$type2 == "cc" && !is.na(opt$s2)) ds2$s <- opt$s2
# eQTLGen uses inverse-normal transformed expression → sdY = 1 exactly.
# coloc.abf() requires sdY (or MAF+N) for type="quant"; set it when z-score
# mode is active and no MAF is available.
if (!is.null(opt$ref_zscore_col) && opt$type2 == "quant" && is.null(ds2$MAF))
    ds2$sdY <- 1

## ---- filter rows with invalid beta/se ---------------------------------------
ok <- is.finite(ds1$beta) & is.finite(ds1$varbeta) & ds1$varbeta > 0 &
      is.finite(ds2$beta) & is.finite(ds2$varbeta) & ds2$varbeta > 0
if (!any(ok)) note_exit("no SNPs with finite beta/varbeta after filtering")
for (nm in c("snp", "beta", "varbeta", "MAF")) {
    if (nm %in% names(ds1)) ds1[[nm]] <- ds1[[nm]][ok]
    if (nm %in% names(ds2)) ds2[[nm]] <- ds2[[nm]][ok]
}
if (sum(ok) < opt$min_snps)
    note_exit(sprintf("only %d SNPs with finite beta/varbeta (< min_snps %d)", sum(ok), opt$min_snps))

## ---- run coloc.abf ----------------------------------------------------------
result <- tryCatch(
    coloc.abf(ds1, ds2, p1 = opt$p1, p2 = opt$p2, p12 = opt$p12),
    error = function(e) note_exit(paste("coloc.abf error:", conditionMessage(e)))
)

## ---- write outputs ----------------------------------------------------------
smry <- data.table(
    gene_id   = if (!is.null(opt$gene_id)) opt$gene_id else "",
    chr       = opt$chr,
    start     = opt$start,
    end       = opt$end,
    nsnps     = as.numeric(result$summary[["nsnps"]]),
    PP.H0.abf = result$summary[["PP.H0.abf"]],
    PP.H1.abf = result$summary[["PP.H1.abf"]],
    PP.H2.abf = result$summary[["PP.H2.abf"]],
    PP.H3.abf = result$summary[["PP.H3.abf"]],
    PP.H4.abf = result$summary[["PP.H4.abf"]],
    note      = ""
)
fwrite(smry, opt$out_summary, sep = "\t")

fwrite(as.data.table(result$results[, c("snp", "SNP.PP.H4")]), opt$out_snps, sep = "\t")

cat(sprintf("[coloc] chr%d:%d-%d  nsnps=%d  PP.H4.abf=%.4f\n",
            opt$chr, opt$start, opt$end,
            result$summary[["nsnps"]],
            result$summary[["PP.H4.abf"]]))
