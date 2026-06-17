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

## Stage 2 — SuSiE subworkflow (port + the part-1 fixes)

**Files:** `bin/susie.R` 🆕, `modules/local/compute_ld.nf` 🆕, `modules/local/susie.nf` 🆕, `subworkflows/local/susie/main.nf` 🆕, params `loci`/`susie_l`/`susie_max_iter` (Task 1.1), fixtures (tiny PLINK panel + `loci.csv`), tests.

- [ ] `COMPUTE_LD`: `plink --r square` per locus from `params.genotype` → `R` + SNP list. PLINK biocontainer.
- [ ] `bin/susie.R` (port of `/tmp/susie_pipeline/scripts/susie.R`), changes: read **hub** columns (`rsid/pos/beta/p/n`, reuse hub `z`); **re-align z to the LD panel `.bim`** (flip where panel A1 = hub `oa`; complement strand; resolve/drop palindromes); honour `--max_iter`; **per-locus resilience** (exit 0 + note on empty/failed locus, never kill the run).
- [ ] `SUSIE` module + subworkflow (`take: ch_harmonised, genotype, loci`): split loci → `COMPUTE_LD` → `SUSIE` → collect credible sets + PIP plots; `emit: credible_sets, versions`.
- [ ] R container with `susieR`+`data.table`+`optparse` (add to the engine image or use a biocontainer).
- [ ] Gate in `SUMSTATS`: `if ('susie' in mods) SUSIE(HARMONISE.out.harmonised, PREPARE_REFERENCE.out.genotype, ch_loci)`.
- [ ] Tests: `COMPUTE_LD` (tiny panel), `bin/susie.R` re-alignment unit (incl. a deliberate flip + a palindrome drop), subworkflow.

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
