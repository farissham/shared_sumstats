#!/usr/bin/env Rscript
# harmonise.R — unified GWAS summary-statistics harmoniser.
#
# Pipeline (format-agnostic after parsing):
#   parse <format>  ->  resolve rsID (build matcher)  ->  panel join + orient  ->  hub
#
# Each per-format parser normalises its layout to a common intermediate:
#   marker (rsID or chr:pos), ea_in (effect allele), oa_in, eaf_in, beta, se, p, n_row
# then the shared steps below are identical for every format. The canonical hub is:
#   rsid chr pos ea oa eaf beta se p n z   (tab-separated, gzipped)

suppressPackageStartupMessages({ library(data.table); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
    make_option("--input",       type = "character"),
    make_option("--format",      type = "character", help = "gwas-ssf | gwama | hail | hail_step3 | cvdkp"),
    make_option("--snp_info",    type = "character", help = "LD reference snp.info"),
    make_option("--rsid_map",    type = "character", default = NA,
                help = "Optional chr,pos,rsid map (must match the input build) for chr:pos markers."),
    make_option("--chain",       type = "character", default = NA,
                help = "Optional UCSC .chain (input build -> panel build); lifts chr:pos markers before resolution."),
    make_option("--build",       type = "character", default = NA,
                help = "Declared genome build of chr:pos markers (e.g. GRCh37/GRCh38)."),
    make_option("--panel_build", type = "character", default = "GRCh37",
                help = "Genome build of snp.info (default GRCh37)."),
    make_option("--min_match",   type = "double", default = 0.01,
                help = "Min chr:pos match rate before flagging a build mismatch."),
    make_option("--build_check_min", type = "integer", default = 20L,
                help = "Only apply the build-mismatch check when at least this many chr:pos markers are present."),
    make_option("--palindrome_maf", type = "double", default = 0.42,
                help = "Resolve palindromic SNPs by frequency only when MAF <= this; else drop."),
    make_option("--assume_forward_strand", action = "store_true", default = FALSE,
                help = "Keep unresolved palindromes on the forward-strand assumption instead of dropping."),
    make_option("--n",           type = "integer",   help = "Fallback sample size."),
    make_option("--output",      type = "character", help = "Output harmonised .tsv.gz path.")
)))
stopifnot(file.exists(opt$input), file.exists(opt$snp_info), !is.null(opt$n), !is.null(opt$output))

read_any <- function(path) if (grepl("\\.gz$", path, ignore.case = TRUE))
    fread(cmd = sprintf("zcat %s", shQuote(path))) else fread(path)
need <- function(d, cols) { miss <- setdiff(cols, names(d))
    if (length(miss)) stop(sprintf("[%s] missing columns: %s", opt$format, paste(miss, collapse = ", "))) }

# Minimal pure-R UCSC chain liftover (no external tools). Lifts a table with
# integer CHR/POS (1-based) from the chain's source (target) assembly to its
# destination (query) assembly. Handles +/- query strand; positions in gaps or
# on chromosomes/blocks the chain doesn't cover are dropped.
liftover_chain <- function(pos_dt, chain_path) {
    stopifnot(file.exists(chain_path))
    ln <- readLines(chain_path); ln <- ln[nzchar(trimws(ln))]
    n <- length(ln)
    chrom <- character(n); bstart <- numeric(n); bend <- numeric(n)
    qlo <- numeric(n); qstr <- character(n); qsz <- numeric(n); qnm <- character(n)
    j <- 0L; ct <- 0; cq <- 0; tN <- qN <- qStrand <- NA_character_; qSize <- 0
    for (k in seq_len(n)) {
        f <- strsplit(trimws(ln[k]), "[[:space:]]+")[[1]]
        if (f[1] == "chain") {
            tN <- sub("^chr", "", f[3]); ct <- as.numeric(f[6])
            qN <- sub("^chr", "", f[8]); qSize <- as.numeric(f[9]); qStrand <- f[10]; cq <- as.numeric(f[11])
        } else {
            p <- as.numeric(f); size <- p[1]
            j <- j + 1L
            chrom[j] <- tN; bstart[j] <- ct; bend[j] <- ct + size - 1
            qlo[j] <- cq; qstr[j] <- qStrand; qsz[j] <- qSize; qnm[j] <- qN
            ct <- ct + size + (if (length(p) >= 2) p[2] else 0)
            cq <- cq + size + (if (length(p) >= 3) p[3] else 0)
        }
    }
    if (j == 0L) return(pos_dt[0L])
    blocks <- data.table(chrom = chrom[1:j], bstart = bstart[1:j], bend = bend[1:j],
                         qlo = qlo[1:j], qstr = qstr[1:j], qsz = qsz[1:j], qnm = qnm[1:j])
    setkey(blocks, chrom, bstart, bend)
    q <- data.table(chrom = as.character(pos_dt$CHR),
                    pstart = as.numeric(pos_dt$POS) - 1, .idx = seq_len(nrow(pos_dt)))
    q[, pend := pstart]
    ov <- foverlaps(q, blocks, by.x = c("chrom", "pstart", "pend"),
                    by.y = c("chrom", "bstart", "bend"), type = "within", nomatch = NULL)
    ov[, q0 := fifelse(qstr == "+", qlo + (pstart - bstart), qsz - 1 - (qlo + (pstart - bstart)))]
    res <- pos_dt[ov$.idx]
    res[, CHR := suppressWarnings(as.integer(ov$qnm))]
    res[, POS := as.integer(ov$q0 + 1)]
    res[!is.na(CHR) & !is.na(POS)]
}

# ----------------------------- per-format parsers -----------------------------
# Each returns: marker, ea_in (effect allele), oa_in, eaf_in, beta, se, p, n_row

parse_gwas_ssf <- function(d) {
    if (!"rsid" %in% names(d)) {
        alt <- intersect(c("variant_id", "rs_id", "SNP", "snp"), names(d))
        if (length(alt)) setnames(d, alt[1], "rsid")
    }
    need(d, c("rsid", "effect_allele", "other_allele", "effect_allele_frequency",
              "beta", "standard_error", "p_value"))
    ncase <- grep("^n_case$",  names(d), ignore.case = TRUE, value = TRUE)[1]
    ntot  <- grep("^n_total$", names(d), ignore.case = TRUE, value = TRUE)[1]
    ncol  <- grep("^n$",       names(d), ignore.case = TRUE, value = TRUE)[1]
    n_row <- if (!is.na(ncase) && !is.na(ntot)) {
        nc <- as.numeric(d[[ncase]]); nt <- as.numeric(d[[ntot]])
        ne <- 4 * nc * (nt - nc) / nt; ne[!is.finite(ne) | ne <= 0] <- NA_real_; ne
    } else if (!is.na(ntot)) as.numeric(d[[ntot]])
      else if (!is.na(ncol)) as.numeric(d[[ncol]])
      else NA_real_
    data.table(marker = as.character(d$rsid), ea_in = d$effect_allele, oa_in = d$other_allele,
               eaf_in = suppressWarnings(as.numeric(d$effect_allele_frequency)),
               beta = as.numeric(d$beta), se = as.numeric(d$standard_error),
               p = as.numeric(d$p_value), n_row = n_row)
}

parse_gwama <- function(d) {
    need(d, c("reference_allele", "other_allele", "eaf", "beta", "se", "p-value"))
    marker <- if ("rsid" %in% names(d)) as.character(d$rsid) else as.character(d$rs_number)
    data.table(marker = marker, ea_in = d$reference_allele, oa_in = d$other_allele,
               eaf_in = suppressWarnings(as.numeric(d$eaf)),
               beta = as.numeric(d$beta), se = as.numeric(d$se),
               p = as.numeric(d[["p-value"]]), n_row = NA_real_)
}

parse_hail <- function(d) {
    need(d, c("locus", "alleles", "beta", "standard_error", "p_value", "AF"))
    loc <- sub("^chr", "", d$locus)                       # "chr1:822354" -> "1:822354"
    ac  <- gsub('[]["]', "", d$alleles)                   # '["C","T"]' -> "C,T"
    ref <- tstrsplit(ac, ",", fixed = TRUE)[[1]]
    alt <- tstrsplit(ac, ",", fixed = TRUE)[[2]]
    data.table(marker = loc, ea_in = alt, oa_in = ref,    # effect = ALT; beta/AF are per ALT
               eaf_in = suppressWarnings(as.numeric(d$AF)),
               beta = as.numeric(d$beta), se = as.numeric(d$standard_error),
               p = as.numeric(d$p_value), n_row = NA_real_)
}

parse_hail_step3 <- function(d) {
    # hail_gwas step3 export: separate CHR/POS/REF/ALT (not combined locus/alleles),
    # and a chi_sq_stat instead of standard_error - se is recovered via the Wald
    # relationship chi_sq = (beta/se)^2, so se = |beta| / sqrt(chi_sq_stat).
    need(d, c("CHR", "POS", "REF", "ALT", "ALT_AF", "beta", "chi_sq_stat", "p_value"))
    chr   <- sub("^chr", "", as.character(d$CHR))
    se    <- abs(as.numeric(d$beta)) / sqrt(as.numeric(d$chi_sq_stat))
    n_row <- if ("N" %in% names(d)) as.numeric(d$N) else NA_real_
    data.table(marker = paste0(chr, ":", d$POS), ea_in = d$ALT, oa_in = d$REF,
               eaf_in = suppressWarnings(as.numeric(d$ALT_AF)),
               beta = as.numeric(d$beta), se = se,
               p = as.numeric(d$p_value), n_row = n_row)
}

parse_cvdkp <- function(d) {
    fc <- function(p, req = TRUE) { m <- grep(p, names(d), ignore.case = TRUE, value = TRUE, perl = TRUE)
        if (!length(m)) { if (req) stop("[cvdkp] missing column matching: ", p); return(NA_character_) }; m[1] }
    c_rsid <- fc("^rsid$"); c_a1 <- fc("^allele1$"); c_a2 <- fc("^allele2$"); c_f <- fc("^freq1$")
    c_b <- fc("^effect$"); c_se <- fc("^stderr$"); c_lp <- fc("^log\\(p\\)$|^log\\.p$|^logp$")
    c_nc <- fc("^n_events?$|^n_case$", FALSE); c_nt <- fc("^n_total$", FALSE)
    n_row <- if (!is.na(c_nc) && !is.na(c_nt)) {
        nc <- as.numeric(d[[c_nc]]); nt <- as.numeric(d[[c_nt]])
        ne <- 4 * nc * (nt - nc) / nt; ne[!is.finite(ne) | ne <= 0] <- NA_real_; ne
    } else NA_real_
    data.table(marker = as.character(d[[c_rsid]]), ea_in = d[[c_a1]], oa_in = d[[c_a2]],
               eaf_in = suppressWarnings(as.numeric(d[[c_f]])),
               beta = as.numeric(d[[c_b]]), se = as.numeric(d[[c_se]]),
               p = 10 ^ as.numeric(d[[c_lp]]), n_row = n_row)
}

dt <- read_any(opt$input)
s  <- switch(opt$format,
    "gwas-ssf" = parse_gwas_ssf(dt),   "gwama" = parse_gwama(dt),
    "hail"     = parse_hail(dt),       "cvdkp" = parse_cvdkp(dt),
    "hail_step3" = parse_hail_step3(dt),
    stop("Unknown --format: ", opt$format))

# read panel up front: used by both the build matcher and the final rsID join
snpinfo <- fread(opt$snp_info, select = c("Chrom", "ID", "PhysPos", "A1", "A2", "A1Freq"))
setnames(snpinfo, c("Chrom", "ID", "PhysPos"), c("chr", "rsid", "pos"))

# --------------------- resolve rsID (build matcher), per row ---------------------
# rsID-formatted markers are used directly (build-agnostic). chr:pos markers are
# translated to rsID via a position map chosen by build: snp.info itself when the
# input is on the panel build, else an external --rsid_map.
s[, marker := as.character(marker)]
s[, is_rsid := grepl("^rs[0-9]+$", marker)]
s[, is_pos  := !is_rsid & grepl("^[0-9A-Za-z]+:[0-9]+$", marker)]
keep_cols <- c("rsid", "ea_in", "oa_in", "eaf_in", "beta", "se", "p", "n_row")
take <- function(d) if (nrow(d)) d[, ..keep_cols] else NULL

by_rsid <- s[is_rsid == TRUE]
if (nrow(by_rsid)) by_rsid[, rsid := marker]

by_pos <- s[is_pos == TRUE]
if (nrow(by_pos)) {
    by_pos[, c("CHR", "POS") := tstrsplit(marker, ":", fixed = TRUE, type.convert = TRUE)]
    build <- if (is.na(opt$build) || !nzchar(opt$build)) opt$panel_build else opt$build
    # liftover: when chr:pos is on a non-panel build and a chain is supplied, rewrite
    # CHR/POS to the panel build, collapsing this to the same-build case below.
    if (!identical(build, opt$panel_build) && !is.na(opt$chain) && nzchar(opt$chain) && file.exists(opt$chain)) {
        n0 <- nrow(by_pos)
        by_pos <- liftover_chain(by_pos, opt$chain)
        cat(sprintf("[liftover] %s -> %s via %s: %d/%d positions lifted\n",
                    build, opt$panel_build, basename(opt$chain), nrow(by_pos), n0))
        build <- opt$panel_build
    }
    if (!is.na(opt$rsid_map) && file.exists(opt$rsid_map)) {
        cat(sprintf("[build-match] external rsid map: %s\n", opt$rsid_map))
        rmap <- fread(opt$rsid_map, select = c("CHR", "POS", "rsid"))
    } else if (identical(build, opt$panel_build)) {
        cat(sprintf("[build-match] input build '%s' == panel build; map derived from snp.info\n", build))
        rmap <- snpinfo[, .(CHR = chr, POS = pos, rsid = rsid)]
    } else {
        stop(sprintf(paste0("chr:pos markers on build '%s' but the panel is '%s' and no --rsid_map. ",
                            "Supply a build-matched map or liftover to %s."),
                     build, opt$panel_build, opt$panel_build))
    }
    n_pos <- nrow(by_pos)
    # CHR type can disagree between sides (e.g. a genome-wide input with sex
    # chromosomes -> tstrsplit keeps its CHR character, vs an autosome-only
    # reference panel -> fread auto-detects integer) - data.table refuses to
    # join mismatched types outright. Force both to character so a real X/Y
    # chromosome value never gets silently coerced to NA on either side.
    rmap[, CHR := as.character(CHR)]; by_pos[, CHR := as.character(CHR)]
    setkey(rmap, CHR, POS); setkey(by_pos, CHR, POS)
    by_pos <- rmap[by_pos, nomatch = NULL]
    rate <- if (n_pos > 0) nrow(by_pos) / n_pos else 0
    cat(sprintf("[build-match] %d/%d chr:pos matched (%.2f%%)\n", nrow(by_pos), n_pos, 100 * rate))
    if (n_pos >= opt$build_check_min && rate < opt$min_match)
        stop(sprintf("Only %.2f%% of chr:pos matched the map - the declared build '%s' is likely wrong.",
                     100 * rate, build))
}

s <- rbind(take(by_rsid), take(by_pos))
if (is.null(s) || !nrow(s)) stop("No variants could be resolved to rsID.")

# --------------------- panel join + allele orientation ---------------------
setkey(snpinfo, rsid); setkey(s, rsid)
m <- snpinfo[s, nomatch = NULL]

comp <- function(x) chartr("ACGTacgt", "TGCAtgca", x)
m[, eaf_valid := !is.na(eaf_in) & eaf_in >= 0 & eaf_in <= 1]
m[, biallelic := nchar(ea_in) == 1L & nchar(oa_in) == 1L & nchar(A1) == 1L & nchar(A2) == 1L]
m[, cra := fifelse(biallelic, comp(ea_in), NA_character_)]
m[, coa := fifelse(biallelic, comp(oa_in), NA_character_)]
m[, palindromic := biallelic & (oa_in == cra)]
m[, set_ok := (ea_in == A1 & oa_in == A2) | (ea_in == A2 & oa_in == A1)]
m[, flip := NA]; m[, aclass := NA_character_]

pmaf <- opt$palindrome_maf
m[palindromic & set_ok & eaf_valid &
    pmin(eaf_in, 1 - eaf_in) <= pmaf & pmin(A1Freq, 1 - A1Freq) <= pmaf,
    `:=`(flip = abs(eaf_in - (1 - A1Freq)) < abs(eaf_in - A1Freq), aclass = "palindrome")]
if (isTRUE(opt$assume_forward_strand))
    m[palindromic & set_ok & is.na(flip), `:=`(flip = (ea_in == A2), aclass = "palindrome_assumed")]
m[is.na(flip) & !palindromic & ea_in == A1 & oa_in == A2, `:=`(flip = FALSE, aclass = "same")]
m[is.na(flip) & !palindromic & ea_in == A2 & oa_in == A1, `:=`(flip = TRUE,  aclass = "same")]
m[is.na(flip) & biallelic & !palindromic & cra == A1 & coa == A2, `:=`(flip = FALSE, aclass = "strand")]
m[is.na(flip) & biallelic & !palindromic & cra == A2 & coa == A1, `:=`(flip = TRUE,  aclass = "strand")]

n_pre <- nrow(m); m <- m[!is.na(flip)]
cat(sprintf("[orient] kept %d/%d (same=%d strand=%d palindrome=%d) dropped %d\n",
            nrow(m), n_pre, sum(m$aclass == "same"), sum(m$aclass == "strand"),
            sum(m$aclass %in% c("palindrome", "palindrome_assumed")), n_pre - nrow(m)))

m[, beta_out := fifelse(flip, -beta, beta)]
m[, eaf_out  := fifelse(eaf_valid, fifelse(flip, 1 - eaf_in, eaf_in), A1Freq)]
m[, n_out    := fifelse(!is.na(n_row) & n_row > 0, as.numeric(n_row), as.numeric(opt$n))]

out <- m[, .(rsid, chr, pos, ea = A1, oa = A2, eaf = eaf_out,
             beta = beta_out, se = se, p = p, n = n_out)]
out[, z := beta / se]
out <- out[!is.na(rsid) & nzchar(rsid) & rsid != "."]
out <- out[!is.na(beta) & !is.na(se) & se > 0 & !is.na(n) & n > 0]
setorder(out, chr, pos)

fwrite(out, opt$output, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[harmonise:%s] wrote %d SNPs to %s\n", opt$format, nrow(out), opt$output))
