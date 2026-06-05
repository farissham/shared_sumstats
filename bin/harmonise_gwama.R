#!/usr/bin/env Rscript
# harmonise_gwama.R
# Convert a GWAMA meta-analysis output (_meta.out / _rsids.tsv) into the
# canonical harmonised-sumstats hub:
#   rsid chr pos ea oa eaf beta se p n z   (tab-separated, gzipped)
#
# Alleles/coords are aligned to the LD reference snp.info (joined by rsID):
#   - ea/oa  = snp.info A1/A2 (one common orientation across all traits)
#   - beta/eaf flipped to match ea = A1
#   - chr/pos = snp.info Chrom/PhysPos (carry the panel's build; recorded in meta.build)
#
# rsID resolution: GWAMA's `rs_number` is a chr:pos string whose coordinates do
# not reliably align to the LD reference (build/patch differences), so we route
# through rsID. The input must either already carry an `rsid` column (e.g.
# *_rsids.tsv) or supply --rsid_map (chr,pos,rsid) to look one up by chr:pos.
# N is injected from --n (GWAMA per-SNP n_samples is unreliable).

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--input",    type = "character", help = "GWAMA _meta.out or _rsids.tsv"),
    make_option("--snp_info", type = "character", help = "LD reference snp.info"),
    make_option("--rsid_map", type = "character", default = NA,
                help = "Optional chr,pos,rsid map; required if input lacks an rsid column."),
    make_option("--n",        type = "integer",   help = "Sample size to inject"),
    make_option("--output",   type = "character", help = "Output harmonised .tsv.gz path")
)))

stopifnot(file.exists(opt$input), file.exists(opt$snp_info),
          !is.null(opt$n), !is.null(opt$output))

dt <- if (grepl("\\.gz$", opt$input, ignore.case = TRUE)) {
    fread(cmd = sprintf("zcat %s", shQuote(opt$input)))
} else fread(opt$input)

needed <- c("reference_allele", "other_allele", "eaf", "beta", "se", "p-value")
missing <- setdiff(needed, names(dt))
if (length(missing)) stop("GWAMA sumstats missing required columns: ",
                          paste(missing, collapse = ", "))

# ---- attach rsid: use existing column, else look up via rsid_map on chr:pos ----
if (!"rsid" %in% names(dt)) {
    if (is.na(opt$rsid_map) || !file.exists(opt$rsid_map))
        stop("Input has no 'rsid' column and --rsid_map was not provided or does not exist.")
    if ("rs_number" %in% names(dt) && !all(c("CHR", "POS") %in% names(dt)))
        dt[, c("CHR", "POS") := tstrsplit(rs_number, ":", fixed = TRUE, type.convert = TRUE)]
    rmap <- fread(opt$rsid_map, select = c("CHR", "POS", "rsid"))
    setkey(rmap, CHR, POS); setkey(dt, CHR, POS)
    dt <- rmap[dt, nomatch = NULL]
}

# ---- join to LD ref snp.info on rsID -> attach chr/pos + ref alleles + A1Freq ----
snpinfo <- fread(opt$snp_info, select = c("Chrom", "ID", "PhysPos", "A1", "A2", "A1Freq"))
setnames(snpinfo, c("Chrom", "ID", "PhysPos"), c("chr", "rsid", "pos"))
setkey(snpinfo, rsid); setkey(dt, rsid)
m <- snpinfo[dt, nomatch = NULL]

# ---- orient to ea = ref A1 (flip beta + eaf when reference_allele is A2) ----
# GWAMA's `eaf` is the frequency of reference_allele.
m[, flip := fcase(
    reference_allele == A1 & other_allele == A2, FALSE,
    reference_allele == A2 & other_allele == A1, TRUE
)]
m <- m[!is.na(flip)]                       # drop ambiguous / non-matching alleles
m[, beta_out := fifelse(flip, -beta, beta)]

# eaf: GWAMA frequently writes -9 (or omits eaf) when input studies lack it, so a
# raw flip would emit -9/10 into the hub. Use the study eaf only when it is a valid
# frequency in [0,1]; otherwise fall back to the LD reference A1Freq, which is the
# frequency of A1 == ea (no flip needed). Mirrors preprocess_gwama.R's freq = A1Freq.
m[, eaf_num := suppressWarnings(as.numeric(eaf))]
m[, eaf_out := fifelse(!is.na(eaf_num) & eaf_num >= 0 & eaf_num <= 1,
                       fifelse(flip, 1 - eaf_num, eaf_num),
                       A1Freq)]

out <- m[, .(rsid = rsid,
             chr  = chr,
             pos  = pos,
             ea   = A1,
             oa   = A2,
             eaf  = eaf_out,
             beta = beta_out,
             se   = se,
             p    = `p-value`,
             n    = opt$n)]
out[, z := beta / se]

out <- out[!is.na(rsid) & nzchar(rsid) & rsid != "."]
out <- out[!is.na(beta) & !is.na(se) & se > 0 & !is.na(n) & n > 0]
setorder(out, chr, pos)

fwrite(out, opt$output, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[harmonise:gwama] wrote %d SNPs to %s\n", nrow(out), opt$output))
