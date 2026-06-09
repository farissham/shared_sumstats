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
# rsID resolution + build matching: GWAMA's `rs_number` may hold actual rsIDs or a
# chr:pos string. The rsID join to the panel is build-agnostic, but the chr:pos ->
# rsID step is build-sensitive. The built-in build matcher resolves it by build:
#   - input build == panel build  -> derive the chr:pos->rsID map from snp.info (free)
#   - input build != panel build  -> use the supplied build-matched --rsid_map
# A low match rate is flagged as a likely wrong declared build.
# N is injected from --n (GWAMA per-SNP n_samples is unreliable).

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--input",       type = "character", help = "GWAMA _meta.out or _rsids.tsv"),
    make_option("--snp_info",    type = "character", help = "LD reference snp.info"),
    make_option("--rsid_map",    type = "character", default = NA,
                help = "Optional chr,pos,rsid map (must match the input's build); used for chr:pos rs_number on a non-panel build."),
    make_option("--build",       type = "character", default = NA,
                help = "Declared genome build of the input chr:pos coordinates (e.g. GRCh37/GRCh38)."),
    make_option("--panel_build", type = "character", default = "GRCh37",
                help = "Genome build of the LD reference snp.info (default GRCh37)."),
    make_option("--min_match",   type = "double",    default = 0.01,
                help = "Minimum chr:pos match rate before flagging a build mismatch (default 0.01)."),
    make_option("--palindrome_maf", type = "double", default = 0.42,
                help = "Resolve palindromic (A/T, C/G) SNPs by frequency only when MAF <= this (default 0.42); else drop as unresolvable."),
    make_option("--assume_forward_strand", action = "store_true", default = FALSE,
                help = "Keep unresolved palindromic SNPs assuming forward strand instead of dropping them."),
    make_option("--n",           type = "integer",   help = "Sample size to inject"),
    make_option("--output",      type = "character", help = "Output harmonised .tsv.gz path")
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

# read the LD reference up front: it serves both the same-build chr:pos->rsID map
# (keyed by position) and the final rsID join (keyed by rsid).
snpinfo <- fread(opt$snp_info, select = c("Chrom", "ID", "PhysPos", "A1", "A2", "A1Freq"))
setnames(snpinfo, c("Chrom", "ID", "PhysPos"), c("chr", "rsid", "pos"))

# ---- attach rsid (with build matching) ----
# GWAMA's marker column `rs_number` may hold actual rsIDs OR a chr:pos string,
# depending on the input studies. Resolve in priority order:
#   1. an explicit `rsid` column
#   2. rs_number that already looks like rsIDs (rs#######) -> use directly (build-agnostic)
#   3. otherwise treat rs_number as chr:pos and translate to rsID (build-sensitive)
if (!"rsid" %in% names(dt)) {
    rsnum_is_rsid <- "rs_number" %in% names(dt) &&
        mean(grepl("^rs[0-9]+$", head(dt$rs_number, 1000L))) > 0.5
    if (rsnum_is_rsid) {
        setnames(dt, "rs_number", "rsid")
    } else {
        if (!("rs_number" %in% names(dt)))
            stop("Input lacks both 'rsid' and 'rs_number'; cannot resolve rsID.")
        if (!all(c("CHR", "POS") %in% names(dt)))
            dt[, c("CHR", "POS") := tstrsplit(rs_number, ":", fixed = TRUE, type.convert = TRUE)]

        # --- build matcher: pick the chr:pos -> rsID map by genome build ---
        build <- if (is.na(opt$build) || !nzchar(opt$build)) opt$panel_build else opt$build
        if (!is.na(opt$rsid_map) && file.exists(opt$rsid_map)) {
            cat(sprintf("[build-match] using external rsid map: %s\n", opt$rsid_map))
            rmap <- fread(opt$rsid_map, select = c("CHR", "POS", "rsid"))
        } else if (identical(build, opt$panel_build)) {
            cat(sprintf("[build-match] input build '%s' == panel build; deriving chr:pos->rsID map from snp.info\n", build))
            rmap <- snpinfo[, .(CHR = chr, POS = pos, rsid = rsid)]
        } else {
            stop(sprintf(paste0("rs_number is chr:pos on build '%s' but the panel is '%s' and no --rsid_map was given. ",
                                "Supply a build-matched chr,pos,rsid map (or liftover the input to %s)."),
                         build, opt$panel_build, opt$panel_build))
        }

        n_in <- nrow(dt)
        setkey(rmap, CHR, POS); setkey(dt, CHR, POS)
        dt <- rmap[dt, nomatch = NULL]
        match_rate <- if (n_in > 0) nrow(dt) / n_in else 0
        cat(sprintf("[build-match] %d/%d positions matched (%.2f%%)\n", nrow(dt), n_in, 100 * match_rate))
        if (match_rate < opt$min_match)
            stop(sprintf(paste0("Only %.2f%% of chr:pos matched the map - the declared build '%s' is likely wrong, ",
                                "or the map is on a different build."), 100 * match_rate, build))
    }
}

# ---- join to LD ref snp.info on rsID -> attach chr/pos + ref alleles + A1Freq ----
setkey(snpinfo, rsid); setkey(dt, rsid)
m <- snpinfo[dt, nomatch = NULL]

# ---- orient alleles to the panel (ea = A1) with strand-flip + palindrome handling ----
# `reference_allele` is the effect allele; `eaf` is its frequency. We must express
# every variant as ea=A1, oa=A2 (panel orientation), flipping beta/eaf when needed.
#   * exact match (same strand): ra/oa equal A1/A2 in some order.
#   * opposite strand: the complement of ra/oa equals A1/A2 (recovers reverse-strand
#     inputs that an exact-only match would silently drop).
#   * palindromic SNPs (A/T, C/G): strand is unknowable from alleles, so orientation
#     is inferred from allele frequency (study eaf vs reference A1Freq) and dropped
#     when the MAF is too close to 0.5 to call (unless --assume_forward_strand).
comp <- function(x) chartr("ACGTacgt", "TGCAtgca", x)
m[, ra := reference_allele]
m[, oa := other_allele]
m[, eaf_num   := suppressWarnings(as.numeric(eaf))]
m[, eaf_valid := !is.na(eaf_num) & eaf_num >= 0 & eaf_num <= 1]
m[, biallelic := nchar(ra) == 1L & nchar(oa) == 1L & nchar(A1) == 1L & nchar(A2) == 1L]
m[, cra := fifelse(biallelic, comp(ra), NA_character_)]
m[, coa := fifelse(biallelic, comp(oa), NA_character_)]
m[, palindromic := biallelic & (oa == cra)]          # {A,T} or {C,G}
m[, set_ok := (ra == A1 & oa == A2) | (ra == A2 & oa == A1)]

m[, flip := NA]
m[, aclass := NA_character_]

# palindromic SNPs first: resolve strand by frequency only (never trust the letters)
pal_maf <- opt$palindrome_maf
m[palindromic & set_ok & eaf_valid &
    pmin(eaf_num, 1 - eaf_num) <= pal_maf & pmin(A1Freq, 1 - A1Freq) <= pal_maf,
    `:=`(flip = abs(eaf_num - (1 - A1Freq)) < abs(eaf_num - A1Freq), aclass = "palindrome")]
if (isTRUE(opt$assume_forward_strand))   # opt-in: keep unresolved palindromes as-is
    m[palindromic & set_ok & is.na(flip),
        `:=`(flip = (ra == A2), aclass = "palindrome_assumed")]

# non-palindromic: exact (same-strand) match, then opposite-strand recovery
m[is.na(flip) & !palindromic & ra == A1 & oa == A2, `:=`(flip = FALSE, aclass = "same")]
m[is.na(flip) & !palindromic & ra == A2 & oa == A1, `:=`(flip = TRUE,  aclass = "same")]
m[is.na(flip) & biallelic & !palindromic & cra == A1 & coa == A2, `:=`(flip = FALSE, aclass = "strand")]
m[is.na(flip) & biallelic & !palindromic & cra == A2 & coa == A1, `:=`(flip = TRUE,  aclass = "strand")]

n_pre <- nrow(m)
m <- m[!is.na(flip)]                       # drop mismatches + unresolved palindromes
cat(sprintf("[orient] kept %d/%d  (same=%d strand=%d palindrome=%d), dropped %d\n",
            nrow(m), n_pre, sum(m$aclass %in% c("same")), sum(m$aclass == "strand"),
            sum(m$aclass %in% c("palindrome", "palindrome_assumed")), n_pre - nrow(m)))

m[, beta_out := fifelse(flip, -beta, beta)]
# eaf: prefer the (oriented) study eaf; fall back to A1Freq (= freq of ea==A1) when the
# study eaf is missing/sentinel (e.g. GWAMA writes -9). Mirrors preprocess_gwama.R.
m[, eaf_out := fifelse(eaf_valid, fifelse(flip, 1 - eaf_num, eaf_num), A1Freq)]

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
