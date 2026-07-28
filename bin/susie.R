#!/usr/bin/env Rscript
# susie.R
# Per-locus SuSiE-RSS fine-mapping off the harmonised hub.
#
# Ported from the susie_pipeline prototype and retargeted to the canonical hub
# (rsid chr pos ea oa eaf beta se p n z). The key addition over the prototype is
# a robust join + allele re-alignment of z onto the LD panel's coding:
#
#   * join hub -> panel by rsID (primary), falling back to chr:pos when the
#     panel .bim is not rsID-keyed (PLINK .bim IDs may be rsID or chr:pos[:alleles]);
#   * verify + ALIGN alleles against the panel .bim A1/A2: flip z where the hub
#     effect allele is the panel A2 (incl. reverse strand), drop on allele
#     mismatch, and drop palindromic (A/T, C/G) SNPs whose strand is unknowable;
#   * susie_rss then sees z and R on one consistent allele basis.
#
# --align_only writes the per-SNP alignment audit and exits 0 without susieR, so
# the join/alignment (the correctness-critical part) is testable without the
# fine-mapping container.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--hub",        type = "character", help = "Harmonised hub .tsv(.gz)"),
    make_option("--ld",         type = "character", help = "LD matrix (.ld, square, no header)"),
    make_option("--bim",        type = "character", help = "Locus .bim in LD-matrix order (chr id cm pos A1 A2)"),
    make_option("--chr",        type = "integer",   help = "Chromosome"),
    make_option("--start",      type = "integer",   help = "Locus start (bp)"),
    make_option("--end",        type = "integer",   help = "Locus end (bp)"),
    make_option("--n",          type = "integer",   default = NA, help = "Sample-size fallback if hub n missing"),
    make_option("--l",          type = "integer",   default = 10,   help = "Max causal variants (SuSiE L)"),
    make_option("--max_iter",   type = "integer",   default = 1000, help = "Max SuSiE IBSS iterations"),
    make_option("--out_cs",     type = "character", help = "Output credible-sets CSV"),
    make_option("--out_plot",   type = "character", help = "Output PIP plot PNG"),
    make_option("--out_align",  type = "character", help = "Output per-SNP alignment audit TSV"),
    make_option("--align_only", action = "store_true", default = FALSE,
                help = "Write the alignment audit and exit (no susieR, no fine-mapping)")
)))

comp <- function(x) chartr("ACGTacgt", "TGCAtgca", x)

## ---- read hub, subset to the locus ------------------------------------------
# Decompress .gz via a shell pipe so we don't depend on data.table's optional
# R.utils backend (not present in all images).
hub <- if (grepl("\\.gz$", opt$hub, ignore.case = TRUE)) {
    fread(cmd = sprintf("gzip -dc %s", shQuote(opt$hub)))
} else fread(opt$hub)
need <- c("rsid", "chr", "pos", "ea", "oa", "beta", "se", "p", "n", "z")
miss <- setdiff(need, names(hub))
if (length(miss)) stop("hub missing columns: ", paste(miss, collapse = ", "))
hub <- hub[chr == opt$chr & pos >= opt$start & pos <= opt$end]
hub[, c("ea", "oa") := .(toupper(ea), toupper(oa))]

## ---- read panel .bim (LD-matrix order) --------------------------------------
bim <- fread(opt$bim, header = FALSE,
             col.names = c("bchr", "bid", "bcm", "bpos", "A1", "A2"))
bim[, c("A1", "A2") := .(toupper(A1), toupper(A2))]
bim[, idx := .I]                      # row order == LD matrix row/col order

## ---- JOIN: rsID primary, chr:pos fallback -----------------------------------
# primary: hub.rsid == panel .bim ID
m_rsid <- merge(hub, bim[, .(bid, A1, A2, idx)],
                by.x = "rsid", by.y = "bid")
if (nrow(m_rsid)) m_rsid[, match_type := "rsid"]

# fallback: hub rows whose rsID matched no panel ID -> match on chr:pos
un <- hub[!rsid %in% bim$bid]
m_pos <- merge(un, bim[, .(bchr, bpos, A1, A2, idx)],
               by.x = c("chr", "pos"), by.y = c("bchr", "bpos"))
if (nrow(m_pos)) m_pos[, match_type := "pos"]

m <- rbind(m_rsid, m_pos, fill = TRUE)

## ---- allele verification + re-alignment of z onto the panel coding ----------
if (nrow(m)) {
    m[, action := fcase(
        ea == A1 & oa == A2,                 "keep",   # same coding
        ea == A2 & oa == A1,                 "flip",   # swapped -> flip z
        comp(ea) == A1 & comp(oa) == A2,     "keep",   # reverse strand
        comp(ea) == A2 & comp(oa) == A1,     "flip",   # reverse strand, swapped
        default = "drop")]                             # alleles don't match panel
    m[ea == comp(oa), action := "drop"]                # palindrome: strand unknowable
    m[, z_aln := fifelse(action == "flip", -as.numeric(z), as.numeric(z))]
    setorder(m, idx)
} else {
    m[, c("action", "z_aln", "match_type") := .(character(), numeric(), character())]
}

## ---- alignment audit (always written; the testable artifact) ----------------
audit <- if (nrow(m)) {
    m[, .(rsid, chr, pos, ea, oa, A1, A2, match_type, action, z, z_aln)]
} else {
    data.table(rsid = character(), chr = integer(), pos = integer(),
               ea = character(), oa = character(), A1 = character(), A2 = character(),
               match_type = character(), action = character(),
               z = numeric(), z_aln = numeric())
}
fwrite(audit, opt$out_align, sep = "\t")
cat(sprintf("[susie] locus chr%d:%d-%d  matched=%d  keep/flip/drop=%d/%d/%d\n",
            opt$chr, opt$start, opt$end, nrow(m),
            sum(audit$action == "keep"), sum(audit$action == "flip"),
            sum(audit$action == "drop")))

if (isTRUE(opt$align_only)) {
    # Emit a placeholder credible-sets file so the module's outputs exist; this
    # lets the wiring be integration-tested without the susieR container.
    fwrite(data.table(note = "align_only: alignment audit written, fine-mapping skipped"), opt$out_cs)
    quit(status = 0)
}

## ---- per-locus resilience: never let one locus kill the run -----------------
note <- function(msg) {
    cat("[susie] ", msg, " -> writing note, exiting 0\n", sep = "")
    fwrite(data.table(note = msg), opt$out_cs)
    quit(status = 0)
}

kept <- m[action != "drop"]
if (nrow(kept) < 1) note("no usable SNPs after allele alignment")

## ---- LD matrix, aligned to kept SNPs ----------------------------------------
R <- as.matrix(fread(opt$ld, header = FALSE))
R <- R[kept$idx, kept$idx, drop = FALSE]
R <- (R + t(R)) / 2
diag(R) <- 1.0
R <- R * 0.99 + diag(nrow(R)) * 0.01      # mild conditioning (as in the prototype)

n_eff <- suppressWarnings(max(kept$n, na.rm = TRUE))
if (!is.finite(n_eff) || n_eff <= 0) n_eff <- opt$n

## ---- fine-map (honours --max_iter; resilient to solver errors) --------------
suppressPackageStartupMessages(library(susieR))
result <- tryCatch(
    susie_rss(z = kept$z_aln, R = R, n = n_eff,
              L = opt$l, max_iter = opt$max_iter, tol = 1e-3),
    error = function(e) { note(paste("susie_rss failed:", conditionMessage(e))) }
)

cat("[susie] converged:", isTRUE(result$converged),
    " credible sets:", length(result$sets$cs), "\n")

cs <- summary(result)$cs
annotate <- function(vars_str, col) {
    v <- as.integer(strsplit(as.character(vars_str), ",")[[1]])
    paste(col[v], collapse = ",")
}
if (!is.null(cs)) {
    cs$rsid <- sapply(cs$variable, annotate, col = kept$rsid)
    cs$pos  <- sapply(cs$variable, annotate, col = kept$pos)
    cs$z    <- sapply(cs$variable, annotate, col = kept$z_aln)
    cs$p    <- sapply(cs$variable, annotate, col = kept$p)
    cs$pip  <- sapply(cs$variable, annotate, col = result$pip)
    fwrite(cs, opt$out_cs)
    cat("[susie] credible sets ->", opt$out_cs, "\n")
} else {
    fwrite(data.table(note = "no credible sets found"), opt$out_cs)
}

png(opt$out_plot, width = 1200, height = 600, res = 150)
susie_plot(result, y = "PIP", main = sprintf("chr%d:%d-%d", opt$chr, opt$start, opt$end))
dev.off()
cat("[susie] PIP plot ->", opt$out_plot, "\n")
