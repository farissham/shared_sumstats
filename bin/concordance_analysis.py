#!/usr/bin/env python3
# concordance_analysis.py
# Cohort-level pairwise effect-direction concordance across all traits'
# harmonised hubs, via belowlab/Concordance-Analysis's concordance.py
# (bin/concordance.py, vendored unmodified). Mirrors ldsc_rg.R's interface
# and <2-trait graceful-exit pattern (note_exit(), same as coloc.R/ldsc_h2.R),
# since this is the same "collect every trait, compare all pairs" shape as
# LDSC_RG.
#
# concordance.py takes exactly 2 files per call and expects 9
# whitespace-separated columns by position (MARKERNAME CHR POS EA NEA EAF
# EFFECT STDERR PVAL) - see its own header for why POS's genome build doesn't
# matter here (it's only used for within-file distance grouping, never
# compared across files). This script reformats each trait's canonical hub
# (rsid chr pos ea oa eaf beta se p n z) into that layout once, then invokes
# concordance.py once per pair and parses its stdout into one combined TSV.

import argparse
import csv
import gzip
import itertools
import os
import subprocess
import sys

FIELDNAMES = ["id1", "id2", "pval_thresh", "total", "observed_concordance",
              "expected_concordance", "pvalue", "note"]


def note_exit(out_summary, msg):
    print(f"[concordance_analysis] NOTE: {msg}", flush=True)
    with open(out_summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader()
    sys.exit(0)


def prep_hub(hub_path, out_path):
    """Reformat the canonical hub to concordance.py's expected 9 columns."""
    opener = gzip.open if hub_path.endswith(".gz") else open
    need = ["rsid", "chr", "pos", "ea", "oa", "eaf", "beta", "se", "p"]
    with opener(hub_path, "rt") as fin, open(out_path, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        missing = [c for c in need if c not in idx]
        if missing:
            raise SystemExit(f"hub {hub_path} missing columns: {missing}")
        fout.write("MARKERNAME\tCHR\tPOS\tEA\tNEA\tEAF\tEFFECT\tSTDERR\tPVAL\n")
        for line in fin:
            f = line.rstrip("\n").split("\t")
            fout.write("\t".join(f[idx[c]] for c in need) + "\n")


def grab_value(lines, prefix):
    for line in lines:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return ""


def run_pair(script_dir, file1, file2, pval_thresh, workdir):
    concordance_py = os.path.join(script_dir, "concordance.py")
    proc = subprocess.run(
        ["python3", concordance_py, file1, file2, str(pval_thresh)],
        cwd=workdir, capture_output=True, text=True,
    )
    lines = proc.stdout.splitlines()
    return {
        "total": grab_value(lines, "total"),
        "observed_concordance": grab_value(lines, "concordance"),
        "expected_concordance": grab_value(lines, "expected concordance rate based on simluation"),
        "pvalue": grab_value(lines, "pvalue"),
        "returncode": proc.returncode,
        "stderr": proc.stderr,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hubs", required=True, help="comma-separated harmonised hub paths")
    ap.add_argument("--ids", required=True, help="comma-separated trait ids, same order as --hubs")
    ap.add_argument("--pval_thresh", type=float, default=0.005)
    ap.add_argument("--out_summary", required=True)
    args = ap.parse_args()

    hubs = args.hubs.split(",")
    ids = args.ids.split(",")
    out_summary = args.out_summary

    if len(hubs) != len(ids):
        note_exit(out_summary, f"mismatched --hubs ({len(hubs)}) and --ids ({len(ids)}) counts")
    if len(hubs) < 2:
        note_exit(out_summary, f"only {len(hubs)} trait(s) - concordance needs at least 2")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    workdir = os.getcwd()

    prepped = {}
    for id_, hub in zip(ids, hubs):
        prepped_path = os.path.join(workdir, f"{id_}.concordance_input.tsv")
        prep_hub(hub, prepped_path)
        prepped[id_] = prepped_path

    rows = []
    for id1, id2 in itertools.combinations(ids, 2):
        res = run_pair(script_dir, prepped[id1], prepped[id2], args.pval_thresh, workdir)
        note = "" if res["returncode"] == 0 else f"concordance.py exited {res['returncode']}: {res['stderr'][-500:]}"
        rows.append({
            "id1": id1, "id2": id2, "pval_thresh": args.pval_thresh,
            "total": res["total"], "observed_concordance": res["observed_concordance"],
            "expected_concordance": res["expected_concordance"], "pvalue": res["pvalue"],
            "note": note,
        })

    with open(out_summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"[concordance_analysis] {len(ids)} trait(s), {len(rows)} pair(s) written to {out_summary}", flush=True)


if __name__ == "__main__":
    main()
