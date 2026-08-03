#!/bin/bash
# Lift a PLINK panel's CHR/POS columns GRCh37->GRCh38, using the liftOver
# binary + hg19ToHg38.over.chain.gz. Unlike lift_snpinfo.sh (which preserves
# row order via a chr=0 sentinel, because ldm.info's eigen blocks index
# snp.info by exact row position), a PLINK panel has no such external index -
# .bed/.bim/.fam are always kept in sync by PLINK itself, so unmapped SNPs
# can simply be excluded outright via `plink2 --exclude` rather than kept as
# sentinel rows. We still use the same sentinel-then-exclude pattern
# internally (chr=0 never appears in real data), just as an implementation
# detail to keep row-index bookkeeping simple, not because anything external
# needs the row count preserved.
#
# Usage: ./lift_plink_panel.sh <bfile prefix> <liftOver wrapper> <chain file> <output prefix>
set -euo pipefail

BFILE="${1:?input PLINK bfile prefix (no .bed/.bim/.fam)}"
LIFTOVER="${2:?path to liftOver wrapper/binary}"
CHAIN="${3:?path to hg19ToHg38.over.chain.gz}"
OUT="${4:?output PLINK bfile prefix}"

# mktemp's default (system /tmp, PBS-redirected to /var/tmp/pbs.<job>/...) is
# NOT under the Singularity wrapper's --bind /rds/general:/rds/general, so
# liftOver running inside the container can't see it - same bug already hit
# and fixed once for lift_snpinfo.sh. Put the workdir under ephemeral instead.
WORKDIR=$(mktemp -d -p "$(dirname "$OUT")")
trap 'rm -rf "$WORKDIR"' EXIT

echo "[lift_plink_panel] input: ${BFILE}.bim ($(wc -l < "${BFILE}.bim") SNPs)"

# BED-style intervals keyed by original row number (not rsID) - avoids any
# rsID-collision ambiguity when mapping lifted positions back.
awk 'BEGIN{OFS="\t"} { print "chr"$1, $4-1, $4, NR }' "${BFILE}.bim" > "$WORKDIR/snps.hg19.bed"

"$LIFTOVER" "$WORKDIR/snps.hg19.bed" "$CHAIN" "$WORKDIR/snps.hg38.bed" "$WORKDIR/snps.unmapped.bed"

echo "[lift_plink_panel] $(wc -l < "$WORKDIR/snps.hg38.bed") / $(wc -l < "$WORKDIR/snps.hg19.bed") lifted"

# Map: row number -> new chr/pos (strip "chr" prefix; BED end == 1-based pos).
awk 'BEGIN{OFS="\t"} { sub("^chr","",$1); print $4, $1, $3 }' "$WORKDIR/snps.hg38.bed" \
    | sort -k1,1n > "$WORKDIR/lifted.map"

# Rebuild .bim row-for-row (preserves .bed/.fam sync): lifted rows get their
# new chr/pos, unlifted rows get the chr=0 sentinel and go on the exclude list.
# liftOver can legitimately map a GRCh37 position onto a GRCh38 alt-contig/
# patch scaffold (e.g. "1_KI270766v1_alt") that has no GRCh37 equivalent -
# none of this pipeline's other reference data (LD reference, hub) uses
# anything but standard autosomes 1-22, so treat a non-1-22 chromosome code
# the same as a failed lift rather than passing it through to PLINK2 (which
# rejects it outright anyway without --allow-extra-chr).
awk 'BEGIN{OFS="\t"; while((getline l < "'"$WORKDIR"'/lifted.map") > 0) {
        split(l, a, "\t"); newchr[a[1]]=a[2]; newpos[a[1]]=a[3]
    } }
    { ok = (NR in newchr) && (newchr[NR] ~ /^([1-9]|1[0-9]|2[0-2])$/)
      if (ok) { $1=newchr[NR]; $4=newpos[NR]; print }
      else { print $2 > "'"$WORKDIR"'/exclude.txt"; $1=0; $4=0; print } }' \
    "${BFILE}.bim" > "$WORKDIR/lifted.bim"

n_excl=$(wc -l < "$WORKDIR/exclude.txt" 2>/dev/null || echo 0)
echo "[lift_plink_panel] excluding $n_excl unmapped SNPs"

cp "${BFILE}.bed" "$WORKDIR/lifted.bed"
cp "${BFILE}.fam" "$WORKDIR/lifted.fam"

# A per-row liftover keeps the original GRCh37 row order while updating
# chr/pos in place, which can leave a handful of variants out of chromosome-
# contiguous order if liftOver reassigns them across a chromosome boundary.
# PLINK2 requires contiguous chromosome grouping for --make-bed; its own
# suggested remedy is --make-pgen + --sort-vars first, then convert to .bed.
plink2 --bed "$WORKDIR/lifted.bed" --bim "$WORKDIR/lifted.bim" --fam "$WORKDIR/lifted.fam" \
    --exclude "$WORKDIR/exclude.txt" \
    --sort-vars \
    --make-pgen --out "$WORKDIR/sorted"

plink2 --pfile "$WORKDIR/sorted" \
    --make-bed --out "$OUT"

echo "[lift_plink_panel] done: ${OUT}.bim ($(wc -l < "${OUT}.bim") SNPs retained)"
