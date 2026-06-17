library(optparse)
library(data.table)
library(susieR)

option_list <- list(
    make_option("--sumstats",  type="character", help="Path to summary statistics file"),
    make_option("--ld",        type="character", help="Path to LD matrix (.ld)"),
    make_option("--snps",      type="character", help="Path to SNP list file"),
    make_option("--chr",       type="integer",   help="Chromosome"),
    make_option("--start",     type="integer",   help="Locus start position"),
    make_option("--end",       type="integer",   help="Locus end position"),
    make_option("--L",         type="integer",   default=10,   help="Max causal variants [default: 10]"),
    make_option("--max_iter",  type="integer",   default=1000, help="Max SuSiE iterations [default: 1000]"),
    make_option("--out_cs",    type="character", help="Output path for credible sets CSV"),
    make_option("--out_plot",  type="character", help="Output path for PIP plot PNG")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Loading sumstats...\n")
sumstats <- fread(opt$sumstats)
cat("Sumstats columns:", paste(names(sumstats), collapse=", "), "\n")
sumstats[, z_score := as.numeric(A1_beta) / as.numeric(se)]

cat("Subsetting to locus chr", opt$chr, ":", opt$start, "-", opt$end, "\n")
locus_data <- sumstats[chr == opt$chr & pos_b37 >= opt$start & pos_b37 <= opt$end]
cat("Locus SNPs in sumstats:", nrow(locus_data), "\n")

cat("Loading LD matrix and SNP list...\n")
snp_list <- fread(opt$snps, header=FALSE)$V1
R_matrix <- as.matrix(fread(opt$ld, header=FALSE))
cat("LD matrix dimensions:", nrow(R_matrix), "x", ncol(R_matrix), "\n")
cat("SNP list length:", length(snp_list), "\n")

cat("Matching SNPs...\n")
# Order by LD matrix so z-scores and R rows/cols are aligned
common_snps <- intersect(snp_list, locus_data$rsID)
cat("Common SNPs:", length(common_snps), "\n")

if (length(common_snps) == 0) {
    cat("ERROR: No common SNPs found. Check rsID format in sumstats vs LD reference.\n")
    cat("Example sumstats rsIDs:", head(locus_data$rsID, 5), "\n")
    cat("Example LD SNP IDs:", head(snp_list, 5), "\n")
    fwrite(data.table(note="No common SNPs found"), opt$out_cs)
    quit(status=1)
}

snp_idx    <- which(snp_list %in% common_snps)
R_matrix   <- R_matrix[snp_idx, snp_idx]
locus_data <- locus_data[match(common_snps, rsID)]  # align to LD matrix order

cat("Symmetrizing and regularizing LD matrix...\n")
R_matrix <- (R_matrix + t(R_matrix)) / 2
diag(R_matrix) <- 1.0

# Shrink off-diagonals slightly to improve conditioning
R_matrix <- R_matrix * 0.99 + diag(nrow(R_matrix)) * 0.01

cat("Running SuSiE RSS...\n")
result <- susie_rss(
    z        = locus_data$z_score,
    R        = R_matrix,
    n        = max(locus_data$N_total),
    L        = opt$L,
    max_iter = 10000,
    tol      = 1e-3
)

cat("Converged:", result$converged, "\n")
cat("Number of credible sets:", length(result$sets$cs), "\n")

cs <- summary(result)$cs

annotate_vars <- function(vars_str, col) {
    vars <- as.integer(strsplit(as.character(vars_str), ",")[[1]])
    paste(col[vars], collapse=",")
}

if (!is.null(cs)) {
    cs$rsID    <- sapply(cs$variable, annotate_vars, col=locus_data$rsID)
    cs$pos     <- sapply(cs$variable, annotate_vars, col=locus_data$pos_b37)
    cs$z_score <- sapply(cs$variable, annotate_vars, col=locus_data$z_score)
    cs$pval    <- sapply(cs$variable, annotate_vars, col=locus_data$pval)
    fwrite(cs, opt$out_cs)
    cat("Credible sets written to:", opt$out_cs, "\n")
} else {
    fwrite(data.table(note="No credible sets found"), opt$out_cs)
    cat("No credible sets found.\n")
}

png(opt$out_plot, width=1200, height=600, res=150)
susie_plot(result, y="PIP", main=paste0("Chr", opt$chr, ":", opt$start, "-", opt$end))
dev.off()
cat("PIP plot written to:", opt$out_plot, "\n")
