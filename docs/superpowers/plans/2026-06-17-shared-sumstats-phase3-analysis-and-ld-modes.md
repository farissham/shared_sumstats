# shared_sumstats Phase 3 — Analysis Layer (SBayesRC + SuSiE) on a Dual-Mode LD Reference

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. TDD per the existing repo convention (nf-test) wherever the unit is testable without heavy container/data.

**Goal:** Add the analysis layer to `shared_sumstats` — **SBayesRC** (genome-wide Bayesian SNP effects / PGS weights) and **SuSiE** (per-locus fine-mapping → credible sets + PIP) — both consuming the canonical harmonised hub (`<id>.harmonised.tsv.gz`). At the same time, introduce a **`PREPARE_REFERENCE`** stage that resolves the LD reference under **two modes**:

- **Mode A (`prebuilt`, default):** a pre-staged `ld_ref_dir` (or a download of `ld_ref_url`).
- **Mode B (`build`):** construct the eigen-LD **and** `snp.info` from a genotype panel via `SBayesRC::LDstep1-4`.

**Architecture:** `PREPARE_REFERENCE` runs first and emits a *uniform* set of channels regardless of mode — `snp_info`, `ld_dir`, `genotype`, `annot`. `snp_info` **gates HARMONISE** (HARMONISE orients the hub to it; in Mode B it is a build output, so Nextflow auto-orders the build before HARMONISE). The two analysis subworkflows are gated by `params.modules` (`'sbayesrc'`, `'susie'`). SBayesRC consumes `ld_dir` (+ `annot`); SuSiE builds a per-locus `R` from `genotype`.

```
 PREPARE_REFERENCE  ──► snp_info ─► HARMONISE ─► hub ─┬─► SBAYESRC (ld_dir, annot) ─► PGS weights
   (Mode A | Mode B)     ld_dir ───────────────────── │
                         genotype ────────────────────┴─► SUSIE (genotype, loci) ─► credible sets + PIP
                         annot ─────────────────────────► (SBAYESRC MAIN)
```

**Tech Stack:** Nextflow DSL2 (nf-core conventions, `nf-schema`), R (`SBayesRC`, `susieR`, `data.table`, `optparse`), PLINK, nf-test, Docker/Singularity.

**Scope boundary:** This plan is the project's **Phase 3** (SBayesRC), expanded to also land SuSiE and the dual-mode LD reference. Internally **staged**:
- **Stage 1** — `PREPARE_REFERENCE` (Mode A) + SBayesRC port. *(This is "Phase 1" in the scoping discussion — the keystone + the cheap, high-value consumer.)*
- **Stage 2** — SuSiE subworkflow (port + allele re-alignment + per-locus resilience).
- **Stage 3** — Mode B LD build (`SBayesRC::LDstep`).

Out of scope (later project phases): LDSC (Phase 4), report assembler (Phase 5), contribution scaffolding + CI (Phase 6).

**Path constants:**
- `SRC` = `/Volumes/bernard_ssd/extracurricular/bioinformatics/georgelab/nf-core-sumstats` (source to port from)
- `DEST` = `/Volumes/bernard_ssd/extracurricular/bioinformatics/georgelab/draft_sumstats` (this repo; folder is `draft_sumstats`, remote is `bernooi/shared_sumstats`)
- nf-test binary = `/Volumes/bernard_ssd/extracurricular/bioinformatics/georgelab/masters/nf-test`; run with `export JAVA_HOME=/opt/homebrew/opt/openjdk`.

---

## Stage 1 — `PREPARE_REFERENCE` (Mode A) + SBayesRC port

> **Status (2026-06-17): ✅ complete.** Params + schema added; `PREPARE_REFERENCE` (pre-staged + `FETCH_LD_REF` download) wired; `HUB_TO_MA` + `SBAYESRC_TIDY/IMPUTE/MAIN` + `SBAYESRC` subworkflow ported with `versions.yml`; `sumstats.nf` rewired with `params.modules` gate. 15 nf-tests green (incl. new `prepare_reference` + `hub_to_ma`). **Remaining (manual, not unit-testable):** a real `-profile hm3 --modules sbayesrc` run to exercise the ported SBayesRC modules against the full LD reference.

### Task 1.1: Add params + schema entries

**Files:** `nextflow.config` ✏️, `nextflow_schema.json` ✏️

- [ ] Add to `params {}` (near `ld_ref_build`):
  ```groovy
  ld_mode        = 'prebuilt'   // 'prebuilt' (Mode A) | 'build' (Mode B)
  genotype       = null         // PLINK prefix: Mode-B build input AND SuSiE's R panel
  blocks_ref     = null         // LDstep refblock (Block Chrom StartBP EndBP); GRCh37 default
  loci           = null         // SuSiE loci CSV (chr,start,end)
  susie_l        = 10
  susie_max_iter = 1000
  ```
- [ ] Mirror them in `nextflow_schema.json` (reference/analysis groups).
- [ ] **Verify:** `nextflow config -profile test 2>&1 | head -30` parses with no unknown-param error.

### Task 1.2: `PREPARE_REFERENCE` subworkflow (Mode A pre-staged) + rewire `workflows/sumstats.nf` (TDD)

**Files:** `subworkflows/local/prepare_reference/main.nf` 🆕, `workflows/sumstats.nf` ✏️, `tests/subworkflows/local/prepare_reference.nf.test` 🆕

- [ ] Write the failing nf-test: prebuilt mode pointed at `assets/test` emits `snp_info` whose basename is `snp.info` and a non-empty `ld_dir`. (No process invoked → runs without Docker.)
- [ ] Implement `PREPARE_REFERENCE` (param-driven, no `take`):
  - `mode == 'prebuilt'`: `ld_dir = Channel.fromPath(params.ld_ref_dir, type:'dir', checkIfExists:true).first()`; if no dir but `ld_ref_url` set → error pointing at Task 1.3 (`FETCH_LD_REF`).
  - `mode == 'build'`: `error "not implemented (Stage 3)"`.
  - Derive `snp_info = ld_dir.map { file("${it}/snp.info", checkIfExists:true) }`.
  - `annot = params.annot_file ? Channel.fromPath(...).first() : Channel.value([])`.
  - `genotype = params.genotype ? Channel.fromPath("${params.genotype}*").collect() : Channel.value([])`.
  - emit `ld_dir, snp_info, annot, genotype, versions`.
- [ ] Rewire `SUMSTATS`: replace the inline `ch_snp_info` (sumstats.nf:19-22) with `PREPARE_REFERENCE()`; pass `PREPARE_REFERENCE.out.snp_info` to `HARMONISE`; hold `ld_dir`/`annot`/`genotype` for the analysis gate (Task 1.6). Keep emitting `harmonised`/`versions`.
- [ ] **Verify:** new test passes; existing `tests/subworkflows/local/harmonise.nf.test` and `tests/pipeline/main.nf.test` still pass (the pipeline test uses `ld_ref_dir=assets/test`, which Mode-A pre-staged resolves identically).

### Task 1.3: `FETCH_LD_REF` module (Mode A download)

**Files:** `modules/local/fetch_ld_ref.nf` 🆕

- [ ] Process: input `val url`; downloads + unzips to a dir; output `path "<dir>", emit: ld_dir` + `versions.yml`. Wire the `ld_ref_url` branch of `PREPARE_REFERENCE` to it.
- [ ] Container: pin a small wget/unzip image (TODO: confirm tag). **Not** covered by nf-test (network); document a manual smoke against `ld_ref_url`.

### Task 1.4: `HUB_TO_MA` module (TDD)

**Files:** `modules/local/hub_to_ma.nf` 🆕, `tests/modules/local/hub_to_ma.nf.test` 🆕, fixture `assets/test/harmonised_demo.tsv.gz` 🆕 (or generate from the existing harmonise output)

- [ ] Pure column projection (no R needed): `zcat hub | awk` → `.ma` (`SNP A1 A2 freq b se p N` ← hub `rsid ea oa eaf beta se p n`). Container `ghcr.io/bernooi/gctb-sbayesrc:dev`.
- [ ] Test: header is `SNP\tA1\tA2\tfreq\tb\tse\tp\tN`; row count = hub data rows; spot-check one row maps columns correctly.

### Task 1.5: Port `SBAYESRC_TIDY` / `SBAYESRC_IMPUTE` / `SBAYESRC_MAIN`

**Files:** `modules/local/sbayesrc_tidy.nf` 🆕, `modules/local/sbayesrc_impute.nf` 🆕, `modules/local/sbayesrc_main.nf` 🆕 (port from `$SRC/modules/local/`)

- [ ] Near-verbatim port (same container, same `SBayesRC::tidy/impute/sbayesrc` calls). Inputs: `tuple(meta, ma)` + `path ld_ref_dir, stageAs:'ld_ref'` (+ MAIN: `path annot_file`). MAIN keeps `maxForks 1`, `label process_high`, hardcoded `tuneStep`, optional `AnnoPerSnpHsqEnrichment`/`AnnoJointProb` outputs.
- [ ] **Add `versions.yml`** to each (source omitted it): record `SBayesRC` package version.
- [ ] MAIN: make `annot` optional — when the staged file is `[]`, drop the `annot=` arg (supports Mode B / non-EUR with no matching annot).
- [ ] *(No nf-test fixtures — real runs need the full 591-block LD ref + container; tested via a real `-profile hm3` run, matching source reality. Document this limitation.)*

### Task 1.6: `SBAYESRC` subworkflow + wire into `SUMSTATS`

**Files:** `subworkflows/local/sbayesrc/main.nf` 🆕, `workflows/sumstats.nf` ✏️

- [ ] `SBAYESRC` (`take: ch_harmonised, ld_dir, annot`): `HUB_TO_MA → SBAYESRC_TIDY → SBAYESRC_IMPUTE → SBAYESRC_MAIN`; mix all `versions`; `emit: results, versions`.
- [ ] In `SUMSTATS`: `def mods = (params.modules ?: '').tokenize(',')*.trim()`; `if ('sbayesrc' in mods) SBAYESRC(HARMONISE.out.harmonised, PREPARE_REFERENCE.out.ld_dir, PREPARE_REFERENCE.out.annot)`; mix its versions into the emitted `versions`.
- [ ] **Verify:** `nextflow config` parses; `-profile test` (modules empty) still harmonise-only; existing tests green.

### Task 1.7: Stage 1 tests green
- [ ] `nf-test test tests/subworkflows/local/prepare_reference.nf.test tests/modules/local/hub_to_ma.nf.test tests/subworkflows/local/harmonise.nf.test tests/pipeline/main.nf.test` all pass.
- [ ] Commit each task.

**Stage 1 Definition of Done:** `PREPARE_REFERENCE` is the single seam feeding `snp_info` to HARMONISE; `params.modules` gates `sbayesrc`; `HUB_TO_MA` + the three SBayesRC modules exist with `versions.yml`; all feasible nf-tests pass; no maintainer-local paths in committed configs.

---

## Stage 2 — SuSiE fine-mapping (port + the part-1 allele re-alignment)

> **Status (2026-06-22): ✅ ported, retargeted, join built.** `bin/susie.R` reads the hub and does the rsID-primary + chr:pos-fallback join with allele re-alignment (keep/flip/drop) — verified for real in the engine container (rs1 keep, rs2 flip→−z, rs3 via pos-fallback, rs4 palindrome-drop). `COMPUTE_LD` (plink biocontainer) + `SUSIE` (susieR via conda) modules, the `SUSIE_FINEMAP` subworkflow, and the `--modules susie` gate in `sumstats.nf` are wired and tested — **17 nf-tests green** incl. `compute_ld` and a `susie` subworkflow integration test (`--align_only`). **Remaining:** (a) the real `susie_rss` execution needs a susieR env (conda-forge `r-susier`) — wired but not yet run end-to-end (analogous to SBayesRC's real-run gap); (b) README SuSiE section; (c) decide whether to retire Faris's standalone `susie_pipeline/` (collaborator's committed dir — confirm before removing).

**Stage goal:** an opt-in `susie` module that fine-maps the harmonised hub per locus → credible sets + PIP, adding the per-SNP allele re-alignment to the LD panel that the prototype lacks. Unlike SBayesRC, all of this is unit-testable on tiny fixtures.

**Per-trait × per-locus data flow:**
```
 hub --subset locus--> z,ea,oa,eaf (snp.info coding)
 genotype panel --COMPUTE_LD (plink --r)--> R matrix + locus .bim (panel coding)
        -> bin/susie.R: intersect by rsid -> RE-ALIGN z to .bim -> susie_rss -> credible sets + PIP
```

### Task 2.0: Containers (prerequisite)

- [ ] **PLINK** for `COMPUTE_LD`: pin a biocontainer, e.g. `quay.io/biocontainers/plink:1.90b6.21--h779adbc_1`.
- [ ] **R + susieR** for `bin/susie.R`: needs `susieR` + `data.table` + `optparse`. Either (a) build a Wave/mulled biocontainer of the three, or (b) add `susieR` to `ghcr.io/bernooi/gctb-sbayesrc:dev` and re-push. Pin the chosen image; every SuSiE R step uses it.
- [ ] **Verify:** `docker run --rm <image> Rscript -e 'library(susieR); library(data.table); library(optparse)'` exits 0.

### Task 2.1: Fixtures (tiny panel + loci + hand-made LD/bim/hub for the R-only test)

**Files:** `assets/test/susie/loci.csv` 🆕, `assets/test/susie/panel.{bed,bim,fam}` 🆕, `assets/test/susie/locus.ld` 🆕, `assets/test/susie/locus.bim` 🆕, `assets/test/susie/hub_locus.tsv.gz` 🆕

- [ ] **Step 1: `loci.csv`** — `chr,start,end` covering the chr1 fixture SNPs, e.g. `1,1000,5000`.
- [ ] **Step 2: tiny PLINK panel** — hand-write a small text VCF (~60 synthetic samples; SNPs `rs1..rs4` on chr1 at the snp.info positions, with some correlation), then in the PLINK container: `plink --vcf panel.vcf --make-bed --out assets/test/susie/panel`. **Deliberately set the `.bim` A1/A2 so `rs2` is swapped vs the hub and `rs4` is a palindrome** (so COMPUTE_LD + alignment have something to exercise).
- [ ] **Step 3: hand-made unit fixtures for `bin/susie.R`** (no PLINK needed): a 4×4 whitespace `locus.ld`, a matching `locus.bim` (`chr rsid cm pos A1 A2` for `rs1..rs4`, with `rs2` swapped and `rs4` palindromic), and `hub_locus.tsv.gz` (hub columns for `rs1..rs4`, `ea/oa` matching the *unswapped* orientation so the test produces one keep, one flip, one palindrome-drop).
- [ ] **Verify:** `wc -l assets/test/susie/panel.bim` = 4; `gunzip -c .../hub_locus.tsv.gz | head`.
- [ ] **Commit.**

### Task 2.2: `bin/susie.R` — hub remap + re-alignment + resilience (TDD)

**Files:** `bin/susie.R` 🆕, `tests/modules/local/susie.nf.test` 🆕 (added in Task 2.4) — for now drive a **script-level** check.

**CLI:** `--hub --ld --bim --chr --start --end --n --l --max_iter --out_cs --out_plot --out_align`

- [ ] **Step 1: write the alignment audit-driven test first** (script-level, runnable in the susieR container): run `susie.R` on the Task 2.1 hand-made fixtures; assert `out_align` shows `rs1=keep`, `rs2=flip`, `rs4=drop`, and that `out_cs` exists.
- [ ] **Step 2: implement `bin/susie.R`.** Logic:
  1. read hub (gz) → subset to the locus by `chr` & `pos`;
  2. read `.ld` (matrix) and `.bim` (panel SNP order + `A1`/`A2`);
  3. intersect hub ∩ panel by `rsid`, order rows/cols to `.bim`;
  4. **re-align z onto the panel coding** (the crux):
     ```r
     comp <- function(x) chartr("ACGTacgt", "TGCAtgca", x)
     m[, action := fcase(
         ea == A1 & oa == A2,             "keep",
         ea == A2 & oa == A1,             "flip",
         comp(ea) == A1 & comp(oa) == A2, "keep",   # reverse strand
         comp(ea) == A2 & comp(oa) == A1, "flip",   # reverse strand, swapped
         default = "drop")]
     m[ea == comp(oa), action := "drop"]            # palindrome: strand unknowable -> drop (v1)
     m <- m[action != "drop"]
     m[, z_aln := fifelse(action == "flip", -z, z)]
     ```
     Write the per-SNP audit (`rsid, ea, oa, A1, A2, action, z, z_aln`) to `--out_align`.
  5. `susie_rss(z = m$z_aln, R = R[idx, idx], n = max(m$n), L = opt$l, max_iter = opt$max_iter)`;
  6. write credible sets (`--out_cs`) + PIP plot (`--out_plot`);
  7. **per-locus resilience:** if 0 common SNPs or `susie_rss` errors, write a one-line note to `--out_cs` and **`quit(status = 0)`** (never kill the run).
- [ ] **Step 3: verify** the script-level test passes (flip negates `z`, palindrome dropped, CS produced). **Commit.**

> v1 **drops** palindromes (conservative — a wrong sign is worse than a drop). Resolving them by frequency (`plink --freq` panel MAF vs hub `eaf`) is a later refinement; note it in the README.

### Task 2.3: `COMPUTE_LD` module (TDD)

**Files:** `modules/local/compute_ld.nf` 🆕, `tests/modules/local/compute_ld.nf.test` 🆕

- [ ] **Process** (PLINK container), input `tuple(val(locus), path(bed), path(bim), path(fam))` where `locus=[chr,start,end]`; script:
  ```bash
  awk -v c=${chr} -v s=${start} -v e=${end} '$1==c && $4>=s && $4<=e {print $2}' ${bim} > snps.txt
  plink --bfile ${bed.baseName} --extract snps.txt --r square --out locus_${chr}_${start}_${end}
  plink --bfile ${bed.baseName} --extract snps.txt --make-just-bim --out locus_${chr}_${start}_${end}
  ```
  output `tuple(val(locus), path("*.ld"), path("locus_*.bim")), emit: ld` + `versions.yml`. (`--make-just-bim` gives the panel alleles in matrix order for alignment.)
- [ ] **Test:** run on the panel fixture for `1,1000,5000`; assert the `.ld` has N rows = #SNPs in window and a `.bim` is emitted. **Commit.**

### Task 2.4: `SUSIE` module (TDD)

**Files:** `modules/local/susie.nf` 🆕, `tests/modules/local/susie.nf.test` 🆕

- [ ] Wrap `bin/susie.R` (susieR container). Input `tuple(val(meta), path(hub), val(locus), path(ld), path(bim))`; pass `--l ${params.susie_l} --max_iter ${params.susie_max_iter} --n ${meta.n}`. Output `tuple(val(meta), path("*credible_sets*.csv")), emit: credible_sets` + plot + align + `versions.yml`.
- [ ] **Test:** feed `[meta, hub_locus.tsv.gz, [1,1000,5000], locus.ld, locus.bim]`; assert credible-sets CSV and that the align audit marks the flip + drop. **Commit.**

### Task 2.5: `susie` subworkflow + wire into `SUMSTATS` (TDD)

**Files:** `subworkflows/local/susie/main.nf` 🆕, `workflows/sumstats.nf` ✏️, `tests/subworkflows/local/susie.nf.test` 🆕

- [ ] **Subworkflow** (`take: ch_harmonised, genotype, loci_csv`):
  ```groovy
  ch_loci = Channel.fromPath(loci_csv).splitCsv(header:true)
      .map { row -> [ [row.chr as int, row.start as int, row.end as int] ] }
  COMPUTE_LD(ch_loci.combine(genotype))            // genotype = staged bed/bim/fam
  // fan loci LD across every trait:
  ch_in = ch_harmonised.combine(COMPUTE_LD.out.ld) // [meta, hub, locus, ld, bim]
      .map { meta, hub, locus, ld, bim -> [ meta, hub, locus, ld, bim ] }
  SUSIE(ch_in)
  emit: credible_sets = SUSIE.out.credible_sets; versions = ...
  ```
- [ ] **Wire** in `SUMSTATS`: `if ('susie' in mods) { if(!params.loci||!params.genotype) error '...'; SUSIE_SWF(HARMONISE.out.harmonised, PREPARE_REFERENCE.out.genotype, file(params.loci)) }`; mix versions; add a `susie` emit.
- [ ] **Test:** subworkflow on the panel + loci fixtures → one credible-sets file for the locus. **Commit.**

### Task 2.6: Docs + Definition of Done

- [ ] README: SuSiE section (inputs `--genotype` + `--loci`; Mode-A vs Mode-B alignment; the palindrome-drop + LD-vs-GWAS-sample caveats).
- [ ] Plan status banner; **DoD:** `--modules susie` yields per-locus credible sets on the fixtures; the alignment unit proves `keep/flip/drop`; per-locus resilience proven (an empty locus does not fail the run); all new nf-tests green.

---

## Stage 3 — Mode B LD build (`SBayesRC::LDstep`)

**Files:** `modules/local/ldstep1.nf` 🆕, `modules/local/ldstep_block.nf` 🆕, `modules/local/ldstep_merge.nf` 🆕, `PREPARE_REFERENCE` build branch ✏️, `conf/snp_set_build.config` 🆕, tests.

- [ ] `LDSTEP1`: `SBayesRC::LDstep1(mafile=<panel SNP list>, genoPrefix, outDir, blockRef=params.blocks_ref)` → block snplists + `ldm.info`. **Scope SNPs from the genotype panel, not the post-harmonise `.ma`** (keeps the build upstream of HARMONISE; avoids the circular dependency).
- [ ] `LDSTEP_BLOCK` (scatter over `Channel.of(1..N_BLOCK)`): `LDstep2` (full block LD) then `LDstep3` (eigen) per block → `block*.eigen.bin`.
- [ ] `LDSTEP_MERGE`: `LDstep4` → final `snp.info` + `ldm.info`; emit the built `ld_dir`.
- [ ] `PREPARE_REFERENCE` `build` branch: run the three; `ld_dir`/`snp_info` from the merge; `genotype` = the build panel.
- [ ] Handle build-version consistency: `blocks_ref` coordinates must match the genotype build (ties to HARMONISE's liftover handling).
- [ ] Test with a tiny synthetic PLINK panel + a 1–2 block `blocks_ref`.

---

## Open flags (decide at build time, not now)

- **Annotations in Mode B:** a custom panel may have no matching BaselineLD annot → `SBAYESRC_MAIN` runs annotation-free (`annot=[]`). Wired as optional in Task 1.5.
- **`LDstep1` SNP scoping:** feed the panel's SNP list (not the harmonised `.ma`) so the build stays upstream of HARMONISE.
- **SuSiE LD-vs-GWAS-sample match:** a separate statistical concern from allele coding; out of scope here but worth a README note.
- **SBayesRC sampler params** (`sbayesrc_chain_length` etc.) are defined but not passed (parity with source); optionally wire later.
