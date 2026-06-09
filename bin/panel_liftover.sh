#!/usr/bin/env bash
# panel_liftover.sh SNP_INFO CHAIN TARGET_BUILD
# Lift the LD reference panel (snp.info, on the panel build) to TARGET_BUILD with
# UCSC liftOver, emitting a build-matched chr,pos,rsid map that the harmoniser's
# build matcher consumes via --rsid_map. Operating on the panel (clean SNVs) keeps
# liftOver on its most accurate footing and is reusable across all inputs of that build.
set -euo pipefail
snp_info="$1"; chain="$2"; build="$3"

# snp.info -> BED (UCSC needs a chr prefix; BED start is 0-based). name = rsID.
#   snp.info cols: Chrom ID Index GenPos PhysPos A1 A2 A1Freq N Block
awk 'NR>1 { print "chr"$1"\t"($5-1)"\t"$5"\t"$2 }' "$snp_info" > panel.bed

liftOver panel.bed "$chain" lifted.bed unmapped.bed

# lifted BED -> chr,pos,rsid (strip chr prefix; 1-based pos = BED end col)
{ printf 'CHR\tPOS\trsid\n'
  awk '{ sub(/^chr/,"",$1); print $1"\t"$3"\t"$4 }' lifted.bed
} > "${build}.rsid_map.tsv"

n_in=$(($(wc -l < panel.bed)))
n_out=$(($(wc -l < lifted.bed)))
echo "[liftover-panel] ${build}: lifted ${n_out}/${n_in} panel SNPs ($(wc -l < unmapped.bed) unmapped lines)"
