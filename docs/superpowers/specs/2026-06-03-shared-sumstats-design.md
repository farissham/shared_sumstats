# shared_sumstats — Design Spec

- **Date:** 2026-06-03
- **Status:** Draft (awaiting review)
- **Target repo:** `bernooi/shared_sumstats` (currently empty)
- **Source repo:** `bernooi/sumstats` (a.k.a. `georgelab/sumstats`) — existing nf-core-style pipeline being ported from

## 1. Context & goal

`bernooi/sumstats` is a single-maintainer Nextflow / nf-core-style pipeline for descriptive analysis of GWAS summary statistics. It munges each trait's raw sumstats to COJO `.ma`, then runs SBayesRC and LDSC and renders a per-trait HTML report.

`shared_sumstats` is a **fresh, shared, multi-contributor** version of the same pipeline. Other lab members will contribute additional analysis tools (e.g. MAGMA, COLOC, fine-mapping, MR, PRS). The goal of this port is therefore **not** simply to copy code, but to **harden the contribution interface** so a new tool can be added without touching shared internals.

The chosen shape (confirmed with the maintainer): **a single shared Nextflow pipeline**, where contributors add **analysis modules** that plug into a common, well-defined interface.

## 2. Design decisions (resolved)

| # | Decision | Choice |
|---|----------|--------|
| D1 | Repo shape | Single shared Nextflow pipeline (not a loose monorepo) |
| D2 | Contribution interface | Canonical **harmonised hub** + per-module **report fragments** (most extensible option) |
| D3 | Module independence | Modules run in **parallel off the hub**; no inter-module dependencies in v1 |
| D4 | Canonical join key | **rsID** (matches the source pipeline's existing build-agnostic approach) |
| D5 | Genome build | **Record-only.** `chr`/`pos` are descriptive columns inherited from the LD ref's `snp.info` (so they carry the panel's build, GRCh37/hg19); build recorded in `meta.build`. **No liftover in v1.** |
| D6 | Report fragment tech | **Self-contained HTML fragment + `metrics.json`** (language-agnostic; not Quarto-coupled) |
| D7 | Module registration | Explicit `--modules` comma-list with `if (module in enabled)` blocks in the workflow (no dynamic registry yet — YAGNI) |

## 3. Architecture

Three layers with one stable interface (the **harmonised hub**) between input handling and analysis.

```
 samplesheet ──►  [1] HARMONISE  ──►  ===== HARMONISED HUB =====  ──►  [2] ANALYSIS MODULES  ──►  [3] REPORT
 (any format)     format-branch +      <id>.harmonised.tsv.gz          self-contained, parallel    assemble
                  preprocessors        (THE contract)                   subworkflows off the hub    fragments

                                              ┌─────────────────────────┼───────────────────────────┐
                                          SBAYESRC                     LDSC                      <contributed>
                                       (hub→.ma internally)      (hub→.sumstats.gz)            (hub→whatever)
                                       results + fragment         results + fragment            results + fragment
                                                                        │
                                                                        ▼  collect fragments per trait
                                                                   <id>.report.html
```

Two distinct branch points (kept conceptually separate):

- **Inside HARMONISE (pre-hub):** branch by *input file format* (`gwas-ssf`/`gwama`/`hail`/`cvdkp`), normalise, converge to one canonical artifact. This is the expensive, fiddly part (allele alignment, rsID mapping, column normalisation, N injection) and is written **once**.
- **After the hub:** fan out by *analysis tool*. Each tool does only a light "hub → my format" conversion, runs, and emits a report fragment.

This means a contributor **never touches format handling** — they always receive the same clean `harmonised.tsv.gz`.

## 4. The harmonised-sumstats hub spec

File: `<id>.harmonised.tsv.gz`, one row per SNP. Deliberate **column superset** so coordinate-based tools (COLOC, fine-mapping, MR) get fields the current COJO `.ma` discards — notably `chr`/`pos`. (Columns vs SNP *coverage* are different things; see the coverage note below.)

| column | required | description |
|--------|----------|-------------|
| `rsid` | yes | dbSNP rsID — **the canonical join key** |
| `chr`  | yes | chromosome (1–22, X) — descriptive, from LD ref `snp.info` |
| `pos`  | yes | base-pair position in `meta.build` — descriptive, from LD ref `snp.info` |
| `ea`   | yes | **effect allele** (== COJO A1); `beta` is its effect |
| `oa`   | yes | other / non-effect allele (== COJO A2) |
| `eaf`  | yes | effect-allele frequency |
| `beta` | yes | per-allele effect size for `ea` |
| `se`   | yes | standard error of `beta` |
| `p`    | yes | p-value |
| `n`    | yes | per-SNP sample size (falls back to `meta.n`) |
| `z`    | no  | signed Z; derived as `beta/se` if absent (convenience for LDSC) |
| `info` | no  | imputation quality if available |

`meta` (per trait) carries: `id, trait, ancestry, build, n, pop_prev, format`.

**Build policy (D4/D5):** matching is by rsID, so the pipeline is build-agnostic at the join. `chr`/`pos` are populated by joining to the LD ref's `snp.info` on rsID (the same trick `preprocess_gwama.R` already uses), so they automatically carry the reference panel's build. `meta.build` records it for downstream consumers. No liftover.

**SNP coverage (v1):** harmonisation aligns to the LD ref's `snp.info` by rsID for **all** input formats (this is how `chr`/`pos` and canonical `ea`/`oa`/`eaf` are attached). The hub is therefore restricted to reference SNPs — exactly what SBayesRC and LDSC already require. Tools needing denser SNP sets (full-density fine-mapping / COLOC) are a future extension that may add a second, fuller harmonised artifact; v1 does not provide it.

Defining `ea`/`oa` explicitly (which allele the effect refers to) removes the single most common sumstats-harmonisation bug.

Example per-module conversions:

- **SBayesRC `.ma`:** `SNP←rsid, A1←ea, A2←oa, freq←eaf, b←beta, se, p, N←n`
- **LDSC:** `SNP←rsid, A1←ea, A2←oa, signed Z, P←p, N←n` + merge-alleles snplist

## 5. The module contract

An analysis module is a **subworkflow** at `subworkflows/local/<name>/main.nf` (workflow name uppercase, e.g. `SBAYESRC`). It must:

**take:**
- `ch_harmonised` = `[ meta, path(harmonised_tsv_gz) ]`
- its own reference channels (declared as params in `nextflow_schema.json`)

**do:**
- convert the hub internally to whatever its tool needs (`.ma`, `.sumstats.gz`, plink, …) — never push its private format onto the hub
- run its tool(s)

**emit (required):**
- `results` — `[ meta, files… ]` → published to `results/<name>/`
- `report` — `[ meta, path("<id>.<name>.fragment.html"), path("<id>.<name>.metrics.json") ]` → consumed by the assembler
- `versions` — `versions.yml` (tool versions, for provenance)

**ship (required):**
- `meta.yml` (description, params, reference requirements)
- an `nf-test`
- a fragment template / renderer

**register:**
- one entry so it is opt-in via `--modules <name>`; the workflow calls it inside an `if (module in enabled)` block

This kills the two coupling problems in the source pipeline: no module touches the report's signature (they emit fragments), and no module is constrained by COJO `.ma` (they start from the richer hub).

**Independence (D3):** modules run in parallel off the hub; none consumes another's output. If a future tool must chain off another (e.g. fine-mapping consuming an LDSC result), add an explicit inter-module dependency mechanism rather than forcing it through the hub. Out of scope for v1.

## 6. Report model

Each module emits two small things instead of feeding a monolithic report:

- `<id>.<name>.fragment.html` — a self-contained section (headings + plots inlined as base64). The contributor owns the visualisation in any language that can emit HTML.
- `<id>.<name>.metrics.json` — flat headline numbers (e.g. `{"h2": 0.21, "h2_se": 0.03}`) for a uniform cross-module summary table and cross-trait aggregation.

The **assembler** (`subworkflows/local/report/`) is deliberately "dumb" and stable:

1. `groupTuple(by: 0)` all fragments for a trait
2. render a trait header + a summary table built from the `metrics.json`s
3. concatenate fragments in a defined module order (following `--modules`)
4. wrap in one themed HTML shell → `<id>.report.html`

Adding a module never edits the assembler. **This removes the `NO_FILE` sentinel hack** in the source `REPORT`: a trait simply has fewer fragments when fewer modules ran.

**Cohort-level outputs:** LDSC `rg` (genetic correlation) is cross-trait, so it does **not** fit the per-trait fragment model. It becomes a cohort-level artifact: `rg_all_pairs.tsv` plus an optional cohort report page, separate from the per-trait reports. Future cross-trait tools follow the same pattern.

## 7. Workflow orchestration

`workflows/sumstats.nf` shrinks to a thin orchestrator:

```nextflow
HARMONISE(ch_samplesheet)                       // → ch_harmonised [meta, harmonised.tsv.gz]

def enabled = params.modules.tokenize(',')*.trim()
ch_fragments = Channel.empty()
ch_versions  = Channel.empty()

if ('sbayesrc' in enabled) {
    SBAYESRC(ch_harmonised, ch_ld_ref, ch_annot)
    ch_fragments = ch_fragments.mix(SBAYESRC.out.report)
    ch_versions  = ch_versions.mix(SBAYESRC.out.versions)
}
if ('ldsc' in enabled) {
    LDSC(ch_harmonised, ch_ldsc_ld, ch_ldsc_snps)   // per-trait h2 fragment + cohort rg
    ch_fragments = ch_fragments.mix(LDSC.out.report)
    ch_versions  = ch_versions.mix(LDSC.out.versions)
}
// ... one block per registered module ...

REPORT(ch_fragments.groupTuple(by: 0), assembler_template)   // → per-trait html
```

Default: `params.modules = 'sbayesrc,ldsc'`.

## 8. Repo layout

```
shared_sumstats/
├── main.nf, nextflow.config, nextflow_schema.json, modules.json
├── conf/            base, modules, test (tiny portable data), snp_set_{hm3,7m}
├── workflows/sumstats.nf          # HARMONISE → enabled modules → REPORT
├── subworkflows/
│   ├── local/
│   │   ├── harmonise/             # [1] input → harmonised hub
│   │   ├── sbayesrc/              # [2] analysis module
│   │   ├── ldsc/                  # [2] analysis module
│   │   └── report/               # [3] fragment collect + assemble
│   └── nf-core/utils_*           # ported as-is
├── modules/local/                 # leaf processes the subworkflows call
├── bin/                           # R/Python scripts
├── assets/
│   ├── schema_input.json
│   ├── report/                    # assembler template + fragment templates
│   └── module_template/           # cookiecutter for "add a new module"
├── docker/                        # or per-module container definitions
├── docs/
│   ├── usage.md
│   ├── output.md
│   └── contributing_a_module.md   # THE doc collaborators read
└── .github/workflows/             # nf-test (per-module shard) + lint
```

## 9. What to port

**Port wholesale (the backbone — already good):**
- `main.nf`, `workflows/sumstats.nf` (refactored per §7), `subworkflows/` (nf-core `utils_*` + local `PIPELINE_INITIALISATION`/`COMPLETION`)
- `nextflow.config`, `conf/base.config`, `conf/modules.config`, `nextflow_schema.json`, `assets/schema_input.json`, `modules.json`
- CI + hygiene: `.github/workflows/` (nf-test, linting), `.nf-core.yml`, `.pre-commit-config.yaml`, `nf-test.config`, `docs/` scaffold

**Port the science (the tools):**
- 4 preprocessors, SBayesRC ×3, LDSC ×3, the report renderer
- supporting `bin/` scripts + `docker/`

**Port but generalise:** see §11.

**Leave behind (do not seed):**
- checked-in run artifacts: `.nextflow*.log`, `nf-*-reports.tsv`, `results/`, `work/`, `null/` — `.gitignore` them

**New scaffolding (for contributors):**
- `assets/module_template/` cookiecutter (subworkflow + leaf module + `meta.yml` + `nf-test` + fragment template)
- `docs/contributing_a_module.md`
- PR template updated with a module checklist

## 10. Refactor map (source → shared_sumstats)

| source | becomes |
|--------|---------|
| `branch{}` inside `workflows/sumstats.nf` | moves into `subworkflows/local/harmonise/` |
| `munge_sumstats` + `preprocess_{gwama,hail,cvdkp}` → emit `.ma` | same processes, but emit **`harmonised.tsv.gz`** (bin scripts updated); wrapped by harmonise, which now passes `snp.info` to **all** format paths so each can attach `chr`/`pos` and align alleles/freq to the ref |
| `SBAYESRC_TIDY/IMPUTE/MAIN` consume `.ma` | `subworkflows/local/sbayesrc/` consumes hub; adds hub→`.ma` convert; emits `results` + `report` fragment + `versions` |
| `LDSC_MUNGE` (`ma_to_ldsc.py`) | `subworkflows/local/ldsc/`; hub→ldsc convert; `h2`→per-trait fragment; `rg`→cohort output |
| `REPORT` (fixed tuple + `NO_FILE`) | `subworkflows/local/report/`: `groupTuple` fragments → `ASSEMBLE_REPORT` |
| `report.Rmd` (monolith) | split into per-module fragment renderers + a thin assembler template |
| `maxForks 1` on `SBAYESRC_MAIN` | preserved (each MCMC peaks ~5 GB; serialise across traits) |

## 11. Generalisation checklist (de-"bernooi"-fy)

- **Hardcoded paths:** `conf/snp_set_7m.config` has an absolute `/Volumes/bernard_ssd/...` annot path and a `${projectDir}` cache; `conf/test.config` points at the maintainer's local data. Route everything through `--data_dir`.
- **Container images:** `ghcr.io/bernooi/gctb-sbayesrc:dev` and `ghcr.io/bernooi/ldsc:dev` → shared/org namespace, and **pin versions** (no `:dev`) for reproducibility.
- **Manifest / URLs:** source `main.nf` header says `github.com/georgelab/sumstats`; set manifest `name`/`homePage` to the real repo (`bernooi/shared_sumstats`, or a `georgelab` org if one is created — confirm at implementation).
- **Dangling `FETCH_LD_REF`:** referenced in the source README + `snp_set_7m.config` but never implemented. On this clean slate, either implement the download-and-cache process or remove the references.
- **Imperial CX3 profile / `imperial_launch.pbs`:** keep (lab-useful) but parameterise; consider moving institutional bits toward nf-core/configs.

## 12. Testing, CI, governance

- **nf-test per module** (each subworkflow ships its own test); CI shards across modules (the source already has a get-shards action).
- **Portable `test` profile** with tiny bundled data so CI and new contributors can run without the maintainer's local files.
- **Linting:** adopt nf-core lint (source README notes it is not yet compliant) + the existing pre-commit / prettier config.
- **Contribution governance:** `docs/contributing_a_module.md` + `assets/module_template/` define the paved path; PR template enforces the module checklist (contract emits, `meta.yml`, nf-test, fragment).

## 13. Out of scope / future

- Liftover / multi-build harmonisation (v1 is record-only)
- A denser, non-reference-restricted harmonised artifact for full-density tools (fine-mapping / COLOC)
- Inter-module dependencies (v1 modules are independent)
- Dynamic module registry (v1 uses explicit `--modules` blocks)
- Cross-ancestry SBayesRC-multi
- Additional analysis modules themselves (MAGMA, COLOC, MR, PRS) — they are *enabled* by this design, not built by it

## 14. Open questions

- **Org vs personal namespace:** is there (or will there be) a `georgelab` GitHub org for manifest/container naming, or stay under `bernooi`? (Affects §11 only; does not change the architecture.)
- **Container strategy:** one shared image vs per-module images. The source uses two (`gctb-sbayesrc`, `ldsc`). Recommendation: port those two as-is for v1; move toward per-module images as contributors bring tools with conflicting dependencies.
