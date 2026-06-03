# shared_sumstats Phase 1 — Foundation + Harmonise Hub Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the `shared_sumstats` Nextflow pipeline skeleton (ported + generalised from `bernooi/sumstats`) and implement the HARMONISE layer for the `gwas-ssf` format, producing the canonical harmonised-sumstats hub (`<id>.harmonised.tsv.gz`) end-to-end on a tiny bundled fixture.

**Architecture:** Three-layer pipeline (HARMONISE → analysis modules → REPORT). Phase 1 builds the scaffold and the first layer only. Input is branched by `meta.format`; the `gwas-ssf` path joins the input to the LD reference `snp.info` by rsID to attach `chr`/`pos` and orient alleles to a common reference convention, emitting the hub. The other three formats (`gwama`/`hail`/`cvdkp`) are guarded with a loud error and ported in Phase 2.

**Tech Stack:** Nextflow DSL2 (nf-core conventions, `nf-schema` plugin), R (`data.table`, `optparse`), nf-test, Docker/Singularity containers.

**Scope boundary:** Phase 1 = scaffold + `gwas-ssf` harmonise path. Out of scope (later phases): the `gwama`/`hail`/`cvdkp` harmonise paths (Phase 2a — mechanical repeats of Task 8/9 using `bin/preprocess_{gwama,hail,cvdkp}.R` as the basis), the SBayesRC module, the LDSC module, the report assembler, contribution scaffolding/CI.

**Path constants used throughout:**
- `SRC` = `/Volumes/bernard_ssd/extracurricular/bioinformatics/georgelab/nf-core-sumstats` (source repo to port from)
- `DEST` = `/Volumes/bernard_ssd/extracurricular/bioinformatics/georgelab/shared_sumstats` (this repo; already cloned, contains only the committed spec)

All `git` commands run with `-C "$DEST"`. Run `export SRC=... DEST=...` once at the start of each session.

---

### Task 1: Repo skeleton + .gitignore

**Files:**
- Create: `$DEST/.gitignore`
- Create: directory tree (empty dirs are not tracked; created as files land)

- [ ] **Step 1: Set path vars for the session**

```bash
export SRC="/Volumes/bernard_ssd/extracurricular/bioinformatics/georgelab/nf-core-sumstats"
export DEST="/Volumes/bernard_ssd/extracurricular/bioinformatics/georgelab/shared_sumstats"
```

- [ ] **Step 2: Write `.gitignore`** (keeps the leave-behind run artifacts out — spec §9)

Create `$DEST/.gitignore`:

```gitignore
# Nextflow run artifacts
.nextflow*
work/
results/
null/
*.log
nf-*-reports.tsv

# LD reference download cache (imputed_7m profile)
ld_ref_cache/

# OS / editor
.DS_Store
*.swp
```

- [ ] **Step 3: Verify and commit**

Run: `git -C "$DEST" add .gitignore && git -C "$DEST" status -s`
Expected: `A  .gitignore`

```bash
git -C "$DEST" commit -m "chore: add .gitignore for nextflow run artifacts"
```

---

### Task 2: Port nf-core utility subworkflows + modules.json (verbatim)

These are stable nf-core boilerplate (`utils_nextflow_pipeline`, `utils_nfcore_pipeline`, `utils_nfschema_plugin`) that the pipeline-init subworkflow depends on. Copy them unchanged.

**Files:**
- Create: `$DEST/subworkflows/nf-core/` (copied tree)
- Create: `$DEST/modules.json`

- [ ] **Step 1: Copy the nf-core subworkflows tree**

```bash
mkdir -p "$DEST/subworkflows"
cp -R "$SRC/subworkflows/nf-core" "$DEST/subworkflows/nf-core"
```

- [ ] **Step 2: Copy modules.json**

```bash
cp "$SRC/modules.json" "$DEST/modules.json"
```

- [ ] **Step 3: Verify the tree exists**

Run: `find "$DEST/subworkflows/nf-core" -name main.nf | sort`
Expected: three `main.nf` files (utils_nextflow_pipeline, utils_nfcore_pipeline, utils_nfschema_plugin).

- [ ] **Step 4: Commit**

```bash
git -C "$DEST" add subworkflows/nf-core modules.json
git -C "$DEST" commit -m "chore: port nf-core utility subworkflows + modules.json"
```

---

### Task 3: Port + generalise `nextflow.config`, `conf/base.config`, `conf/modules.config`

**Files:**
- Create: `$DEST/nextflow.config` (copied, then edited)
- Create: `$DEST/conf/base.config` (verbatim)
- Create: `$DEST/conf/modules.config` (verbatim)

- [ ] **Step 1: Copy base.config and modules.config verbatim**

```bash
mkdir -p "$DEST/conf"
cp "$SRC/conf/base.config" "$DEST/conf/base.config"
cp "$SRC/conf/modules.config" "$DEST/conf/modules.config"
```

- [ ] **Step 2: Copy nextflow.config**

```bash
cp "$SRC/nextflow.config" "$DEST/nextflow.config"
```

- [ ] **Step 3: Add the two new params** (after the existing `ldsc_snplist` line in the `params {}` block)

In `$DEST/nextflow.config`, find:

```groovy
    // LDSC inputs (set by snp_set profiles). When both are null the LDSC branch is skipped.
    ldsc_ld_dir                = null   // dir of eur_w_ld_chr LD scores
    ldsc_snplist               = null   // w_hm3.snplist for munge_sumstats --merge-alleles
```

Replace with (adds `genome_build` and `modules`):

```groovy
    // LDSC inputs (set by snp_set profiles). When both are null the LDSC branch is skipped.
    ldsc_ld_dir                = null   // dir of eur_w_ld_chr LD scores
    ldsc_snplist               = null   // w_hm3.snplist for munge_sumstats --merge-alleles

    // Genome build recorded in meta.build when a samplesheet row omits it (spec D5: record-only)
    genome_build               = 'GRCh37'

    // Comma-list of analysis modules to run off the harmonised hub. Empty in Phase 1
    // (no analysis modules ported yet). Example later: 'sbayesrc,ldsc'.
    modules                    = ''
```

- [ ] **Step 4: Update the manifest** (point at the new repo)

In `$DEST/nextflow.config`, find `name            = 'georgelab/sumstats'` and replace with:

```groovy
    name            = 'bernooi/shared_sumstats'
```

Find `homePage        = 'https://github.com/bernooi/sumstats'` and replace with:

```groovy
    homePage        = 'https://github.com/bernooi/shared_sumstats'
```

Find the `description     = """nf-core-style pipeline for descriptive analysis of GWAS summary statistics"""` line and replace with:

```groovy
    description     = """Shared, multi-contributor pipeline for descriptive analysis of GWAS summary statistics"""
```

> Note: if a `georgelab` GitHub org is later created, change `name`/`homePage` accordingly (spec §11, open question).

- [ ] **Step 5: Verify config parses**

Run: `cd "$DEST" && nextflow config -profile test 2>&1 | head -20`
Expected: prints resolved config without a parse error. (It may warn that `conf/test.config` does not exist yet — that is fixed in Task 7; ignore for now, or run without `-profile test`: `nextflow config 2>&1 | head -20`.)

- [ ] **Step 6: Commit**

```bash
git -C "$DEST" add nextflow.config conf/base.config conf/modules.config
git -C "$DEST" commit -m "feat: port + generalise nextflow.config, base.config, modules.config

- manifest -> bernooi/shared_sumstats
- add params.genome_build (record-only) and params.modules"
```

---

### Task 4: Port + generalise `nextflow_schema.json`

**Files:**
- Create: `$DEST/nextflow_schema.json` (copied, then edited)

- [ ] **Step 1: Copy the schema**

```bash
cp "$SRC/nextflow_schema.json" "$DEST/nextflow_schema.json"
```

- [ ] **Step 2: Update `$id` and `title`**

In `$DEST/nextflow_schema.json`, replace any occurrence of `georgelab/sumstats` with `bernooi/shared_sumstats` (there are typically two: `$id` URL and `title`).

Run to check: `grep -n "georgelab/sumstats" "$DEST/nextflow_schema.json"` → expected: no output after editing.

- [ ] **Step 3: Add `genome_build` and `modules` to the schema**

The source schema groups params under `$defs`. Locate the group that contains `data_dir`/`snp_set` (the reference/input options group) and add these two properties inside that group's `properties` object:

```json
            "genome_build": {
                "type": "string",
                "enum": ["GRCh37", "GRCh38"],
                "default": "GRCh37",
                "description": "Genome build recorded in meta.build when a samplesheet row omits it (record-only; no liftover)."
            },
            "modules": {
                "type": "string",
                "default": "",
                "description": "Comma-separated analysis modules to run off the harmonised hub (e.g. 'sbayesrc,ldsc'). Empty runs harmonise only."
            }
```

- [ ] **Step 4: Validate the JSON is well-formed**

Run: `python3 -c "import json,sys; json.load(open('$DEST/nextflow_schema.json')); print('valid json')"`
Expected: `valid json`

- [ ] **Step 5: Commit**

```bash
git -C "$DEST" add nextflow_schema.json
git -C "$DEST" commit -m "feat: port nextflow_schema.json; add genome_build + modules params"
```

---

### Task 5: Port `assets/schema_input.json` + add optional `build` column

**Files:**
- Create: `$DEST/assets/schema_input.json` (copied, then edited)

- [ ] **Step 1: Copy the samplesheet schema**

```bash
mkdir -p "$DEST/assets"
cp "$SRC/assets/schema_input.json" "$DEST/assets/schema_input.json"
```

- [ ] **Step 2: Update `$id`**

Replace `georgelab/sumstats` with `bernooi/shared_sumstats` in `$DEST/assets/schema_input.json`.

- [ ] **Step 3: Add the optional `build` property**

In `$DEST/assets/schema_input.json`, find the `pop_prev` property block:

```json
            "pop_prev": {
                "type": "number",
                "minimum": 0,
                "maximum": 1,
                "errorMessage": "Population prevalence (pop_prev) must be a number between 0 and 1; leave blank for quantitative traits",
                "meta": ["pop_prev"]
            }
```

Replace it with (adds `build` after it; note the added comma):

```json
            "pop_prev": {
                "type": "number",
                "minimum": 0,
                "maximum": 1,
                "errorMessage": "Population prevalence (pop_prev) must be a number between 0 and 1; leave blank for quantitative traits",
                "meta": ["pop_prev"]
            },
            "build": {
                "type": "string",
                "enum": ["GRCh37", "GRCh38"],
                "errorMessage": "build must be GRCh37 or GRCh38; leave blank to inherit --genome_build",
                "meta": ["build"]
            }
```

(`build` is NOT added to the `required` array — it is optional and defaults to `params.genome_build`, applied in Task 6.)

- [ ] **Step 4: Validate JSON**

Run: `python3 -c "import json; json.load(open('$DEST/assets/schema_input.json')); print('valid json')"`
Expected: `valid json`

- [ ] **Step 5: Commit**

```bash
git -C "$DEST" add assets/schema_input.json
git -C "$DEST" commit -m "feat: port schema_input.json; add optional build column"
```

---

### Task 6: Port pipeline init/completion subworkflow; default `meta.build`

**Files:**
- Create: `$DEST/subworkflows/local/utils_nfcore_shared_sumstats_pipeline/main.nf` (copied from source, then edited)

- [ ] **Step 1: Copy the local pipeline subworkflow (renamed dir)**

```bash
mkdir -p "$DEST/subworkflows/local"
cp -R "$SRC/subworkflows/local/utils_nfcore_sumstats_pipeline" \
      "$DEST/subworkflows/local/utils_nfcore_shared_sumstats_pipeline"
```

- [ ] **Step 2: Default `meta.build` from `params.genome_build`**

In `$DEST/subworkflows/local/utils_nfcore_shared_sumstats_pipeline/main.nf`, find:

```groovy
    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, sumstats ->
            [ meta, file(sumstats, checkIfExists: true) ]
        }
        .set { ch_samplesheet }
```

Replace with (injects build default):

```groovy
    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, sumstats ->
            meta.build = meta.build ?: params.genome_build
            [ meta, file(sumstats, checkIfExists: true) ]
        }
        .set { ch_samplesheet }
```

- [ ] **Step 3: Verify no stale path references**

Run: `grep -n "sumstats" "$DEST/subworkflows/local/utils_nfcore_shared_sumstats_pipeline/main.nf" | grep -i georgelab`
Expected: no output (the file references `workflow.manifest.name` dynamically, so no hardcoded repo name needs changing).

- [ ] **Step 4: Commit**

```bash
git -C "$DEST" add subworkflows/local/utils_nfcore_shared_sumstats_pipeline
git -C "$DEST" commit -m "feat: port pipeline init/completion subworkflow; default meta.build"
```

---

### Task 7: SNP-set (hm3) profile, portable test profile, and tiny fixtures

**Files:**
- Create: `$DEST/conf/snp_set_hm3.config`
- Create: `$DEST/conf/test.config`
- Create: `$DEST/assets/test/snp.info`
- Create: `$DEST/assets/test/gwas_ssf_demo.tsv`
- Create: `$DEST/assets/test/samplesheet.csv`

- [ ] **Step 1: Port the hm3 profile (already data_dir-based, copy verbatim)**

```bash
cp "$SRC/conf/snp_set_hm3.config" "$DEST/conf/snp_set_hm3.config"
```

- [ ] **Step 2: Write a portable `conf/test.config`** (no maintainer-local paths; uses bundled fixtures)

Create `$DEST/conf/test.config`:

```groovy
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Portable test profile: tiny bundled gwas-ssf fixture + tiny snp.info.
    Runs the HARMONISE layer only (no analysis modules).
        nextflow run . -profile test,docker --outdir results
----------------------------------------------------------------------------------------
*/

params {
    config_profile_name        = 'Phase 1 harmonise test'
    config_profile_description = 'Tiny bundled gwas-ssf fixture -> harmonised hub'

    input        = "${projectDir}/assets/test/samplesheet.csv"
    outdir       = "${projectDir}/results"
    genome_build = 'GRCh37'

    // Point the SBayesRC LD-ref dir at the bundled tiny snp.info's directory.
    ld_ref_dir   = "${projectDir}/assets/test"
    modules      = ''
}
```

- [ ] **Step 3: Write the tiny `snp.info` fixture** (same schema as the real file: `Chrom ID Index GenPos PhysPos A1 A2 A1Freq N Block`, tab-separated)

Create `$DEST/assets/test/snp.info` (8 SNPs; `rs7` is intentionally absent so the demo's `rs7` is dropped on join; `rs3` has alleles swapped relative to the demo to exercise the flip path):

```text
Chrom	ID	Index	GenPos	PhysPos	A1	A2	A1Freq	N	Block
1	rs1	0	0	1001	A	G	0.30	20000	1
1	rs2	1	0	2002	C	T	0.45	20000	1
1	rs3	2	0	3003	G	A	0.20	20000	1
1	rs4	3	0	4004	T	C	0.50	20000	1
2	rs5	4	0	5005	A	T	0.10	20000	2
2	rs6	5	0	6006	C	G	0.40	20000	2
2	rs8	6	0	8008	G	C	0.35	20000	2
2	rs9	7	0	9009	A	G	0.25	20000	2
```

- [ ] **Step 4: Write the tiny `gwas-ssf` fixture** (GWAS-SSF columns; tab-separated). `rs1`,`rs2`,`rs4`,`rs5`,`rs6`,`rs8` match cleanly; `rs3` has effect/other swapped vs `snp.info` (tests flip); `rs7` is not in `snp.info` (tests drop).

Create `$DEST/assets/test/gwas_ssf_demo.tsv`:

```text
rsid	chromosome	base_pair_location	effect_allele	other_allele	effect_allele_frequency	beta	standard_error	p_value	n_total
rs1	1	1001	A	G	0.31	0.025	0.004	1.0e-09	10000
rs2	1	2002	C	T	0.44	-0.011	0.005	2.8e-02	10000
rs3	1	3003	A	G	0.81	0.040	0.006	1.0e-11	10000
rs4	1	4004	T	C	0.49	0.003	0.005	5.5e-01	10000
rs5	2	5005	A	T	0.09	-0.030	0.007	1.9e-05	10000
rs6	2	6006	C	G	0.41	0.015	0.004	1.7e-04	10000
rs7	7	7007	A	G	0.22	0.050	0.010	5.0e-07	10000
rs8	2	8008	G	C	0.36	-0.020	0.006	8.4e-04	10000
```

> Expected harmonise outcome on this fixture: 8 input SNPs; `rs7` is absent from `snp.info` and drops on the rsID join → **7 rows out**. `rs3` flips: demo effect allele `A` == ref `A2`, so `ea` stays ref `A1=G`, `beta` sign flips `0.040 -> -0.040`, `eaf` flips `0.81 -> 0.19`.

- [ ] **Step 5: Write the test samplesheet** (absolute path templating handled by the nf-tests in Tasks 8–10; this CSV is for the manual `-profile test` smoke run)

Create `$DEST/assets/test/samplesheet.csv`:

```csv
id,sumstats,format,n,ancestry,trait,pop_prev,build
demo,assets/test/gwas_ssf_demo.tsv,gwas-ssf,10000,EUR,Demo trait,,GRCh37
```

> The relative `sumstats` path resolves against the launch directory for the manual smoke run (run from `$DEST`). Automated tests (Tasks 8–10) use `file("${projectDir}/assets/test/...")` directly and do not depend on this CSV.

- [ ] **Step 6: Commit**

```bash
git -C "$DEST" add conf/snp_set_hm3.config conf/test.config assets/test
git -C "$DEST" commit -m "test: add hm3 profile, portable test profile, tiny harmonise fixtures"
```

---

### Task 8: HARMONISE_GWAS_SSF process + R script (TDD)

**Files:**
- Create: `$DEST/bin/harmonise_gwas_ssf.R`
- Create: `$DEST/modules/local/harmonise_gwas_ssf.nf`
- Test: `$DEST/tests/modules/local/harmonise_gwas_ssf.nf.test`

- [ ] **Step 1: Write the failing nf-test**

Create `$DEST/tests/modules/local/harmonise_gwas_ssf.nf.test`:

```groovy
nextflow_process {

    name "Test HARMONISE_GWAS_SSF"
    script "modules/local/harmonise_gwas_ssf.nf"
    process "HARMONISE_GWAS_SSF"

    test("gwas-ssf demo harmonises to the hub columns") {
        when {
            process {
                """
                input[0] = [
                    [ id:'demo', format:'gwas-ssf', n:10000, ancestry:'EUR', trait:'Demo', build:'GRCh37' ],
                    file("${projectDir}/assets/test/gwas_ssf_demo.tsv")
                ]
                input[1] = file("${projectDir}/assets/test/snp.info")
                """
            }
        }
        then {
            assert process.success
            def lines = path(process.out.harmonised[0][1]).linesGzip
            // canonical hub header
            assert lines[0] == "rsid\tchr\tpos\tea\toa\teaf\tbeta\tse\tp\tn\tz"
            // rs7 is absent from snp.info -> dropped on join (8 input -> 7 data rows)
            def data = lines.findAll { it.startsWith("rs") }
            assert data.size() == 7
            // rs3 flip: ea stays ref A1=G, beta sign flipped to negative
            def rs3 = data.find { it.startsWith("rs3\t") }.split("\t")
            assert rs3[3] == "G"            // ea == ref A1
            assert rs3[4] == "A"            // oa == ref A2
            assert (rs3[6] as double) < 0   // beta flipped negative
        }
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd "$DEST" && nf-test test tests/modules/local/harmonise_gwas_ssf.nf.test`
Expected: FAIL — `modules/local/harmonise_gwas_ssf.nf` does not exist.

- [ ] **Step 3: Write the R script**

Create `$DEST/bin/harmonise_gwas_ssf.R`:

```r
#!/usr/bin/env Rscript
# harmonise_gwas_ssf.R
# Convert a GWAS-SSF / GWAS Catalog harmonised .tsv(.gz) into the canonical
# harmonised-sumstats hub:
#   rsid chr pos ea oa eaf beta se p n z   (tab-separated, gzipped)
#
# Alleles/coords are aligned to the LD reference snp.info (joined by rsID):
#   - ea/oa  = snp.info A1/A2 (one common orientation across all traits)
#   - beta/eaf flipped to match ea = A1
#   - chr/pos = snp.info Chrom/PhysPos (carry the panel's build; recorded in meta.build)

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--input",    type = "character", help = "GWAS-SSF .tsv(.gz) sumstats"),
    make_option("--snp_info", type = "character", help = "LD reference snp.info"),
    make_option("--n",        type = "integer",   help = "Sample size fallback if file has no N"),
    make_option("--output",   type = "character", help = "Output harmonised .tsv.gz path")
)))

stopifnot(file.exists(opt$input), file.exists(opt$snp_info),
          !is.null(opt$n), !is.null(opt$output))

dt <- if (grepl("\\.gz$", opt$input, ignore.case = TRUE)) {
    fread(cmd = sprintf("zcat %s", shQuote(opt$input)))
} else fread(opt$input)

# GWAS-SSF permits either `rsid` or `variant_id` for the rsID column.
if (!"rsid" %in% names(dt)) {
    alt <- intersect(c("variant_id", "rs_id", "SNP", "snp"), names(dt))
    if (length(alt)) setnames(dt, alt[1], "rsid")
}

needed <- c("rsid", "effect_allele", "other_allele",
            "effect_allele_frequency", "beta", "standard_error", "p_value")
missing <- setdiff(needed, names(dt))
if (length(missing)) stop("Sumstats missing required columns: ",
                          paste(missing, collapse = ", "))

# ---- per-SNP N: effective N for case/control, else N_total, else fixed --n ----
n_case_col  <- grep("^n_case$",  names(dt), ignore.case = TRUE, value = TRUE)[1]
n_total_col <- grep("^n_total$", names(dt), ignore.case = TRUE, value = TRUE)[1]
n_col       <- grep("^n$",       names(dt), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(n_case_col) && !is.na(n_total_col)) {
    nc <- as.numeric(dt[[n_case_col]]); nt <- as.numeric(dt[[n_total_col]])
    n_eff <- 4 * nc * (nt - nc) / nt
    n_eff[!is.finite(n_eff) | n_eff <= 0] <- NA_real_
    dt[, n_pipeline := n_eff]
} else if (!is.na(n_total_col)) {
    dt[, n_pipeline := as.numeric(dt[[n_total_col]])]
} else if (!is.na(n_col)) {
    dt[, n_pipeline := as.numeric(dt[[n_col]])]
} else {
    dt[, n_pipeline := opt$n]
}

# ---- join to LD ref snp.info on rsID -> attach chr/pos + ref alleles ----
snpinfo <- fread(opt$snp_info, select = c("Chrom", "ID", "PhysPos", "A1", "A2"))
setnames(snpinfo, c("Chrom", "ID", "PhysPos"), c("chr", "rsid", "pos"))
setkey(snpinfo, rsid); setkey(dt, rsid)
m <- snpinfo[dt, nomatch = NULL]

# ---- orient to ea = ref A1 (flip beta + eaf when input effect allele is A2) ----
m[, flip := fcase(
    effect_allele == A1 & other_allele == A2, FALSE,
    effect_allele == A2 & other_allele == A1, TRUE
)]
m <- m[!is.na(flip)]                       # drop ambiguous / non-matching alleles
m[, beta_out := fifelse(flip, -beta, beta)]
m[, eaf_out  := fifelse(flip, 1 - effect_allele_frequency, effect_allele_frequency)]

out <- m[, .(rsid = rsid,
             chr  = chr,
             pos  = pos,
             ea   = A1,
             oa   = A2,
             eaf  = eaf_out,
             beta = beta_out,
             se   = standard_error,
             p    = p_value,
             n    = n_pipeline)]
out[, z := beta / se]

out <- out[!is.na(rsid) & nzchar(rsid) & rsid != "."]
out <- out[!is.na(beta) & !is.na(se) & se > 0 & !is.na(n) & n > 0]
setorder(out, chr, pos)

fwrite(out, opt$output, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[harmonise:gwas-ssf] wrote %d SNPs to %s\n", nrow(out), opt$output))
```

- [ ] **Step 4: Make the script executable**

```bash
chmod +x "$DEST/bin/harmonise_gwas_ssf.R"
```

- [ ] **Step 5: Write the process**

Create `$DEST/modules/local/harmonise_gwas_ssf.nf`:

```groovy
process HARMONISE_GWAS_SSF {
    tag "$meta.id"
    label 'process_low'
    container 'ghcr.io/bernooi/gctb-sbayesrc:dev'   // TODO Phase 4: shared namespace + pinned tag

    publishDir "${params.outdir}/harmonise", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(sumstats)
    path snp_info

    output:
    tuple val(meta), path("${meta.id}.harmonised.tsv.gz"), emit: harmonised
    path "versions.yml",                                   emit: versions

    script:
    """
    harmonise_gwas_ssf.R \\
        --input ${sumstats} \\
        --snp_info ${snp_info} \\
        --n ${meta.n} \\
        --output ${meta.id}.harmonised.tsv.gz

    printf '"%s":\\n    r-base: %s\\n    data.table: %s\\n' \\
        "${task.process}" \\
        "\$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')" \\
        "\$(Rscript -e 'cat(as.character(packageVersion("data.table")))')" \\
        > versions.yml
    """
}
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `cd "$DEST" && nf-test test tests/modules/local/harmonise_gwas_ssf.nf.test`
Expected: PASS (1 test). Requires Docker available and the container image pullable; if running locally without Docker, run with the `-profile` that exposes a local R with `data.table` installed, or run the script directly as the fast check in Step 7.

- [ ] **Step 7: (Fast local check, optional) run the script directly**

Run:
```bash
cd "$DEST" && Rscript bin/harmonise_gwas_ssf.R \
  --input assets/test/gwas_ssf_demo.tsv \
  --snp_info assets/test/snp.info \
  --n 10000 \
  --output /tmp/demo.harmonised.tsv.gz && zcat /tmp/demo.harmonised.tsv.gz
```
Expected: header `rsid chr pos ea oa eaf beta se p n z` + 7 data rows; `rs3` line shows `G A` alleles and a negative `beta`.

- [ ] **Step 8: Commit**

```bash
git -C "$DEST" add bin/harmonise_gwas_ssf.R modules/local/harmonise_gwas_ssf.nf tests/modules/local/harmonise_gwas_ssf.nf.test
git -C "$DEST" commit -m "feat: HARMONISE_GWAS_SSF process + R script (input -> harmonised hub)"
```

---

### Task 9: HARMONISE subworkflow (branch by format; guard unported formats) (TDD)

**Files:**
- Create: `$DEST/subworkflows/local/harmonise/main.nf`
- Test: `$DEST/tests/subworkflows/local/harmonise.nf.test`

- [ ] **Step 1: Write the failing nf-test**

Create `$DEST/tests/subworkflows/local/harmonise.nf.test`:

```groovy
nextflow_workflow {

    name "Test HARMONISE subworkflow"
    script "subworkflows/local/harmonise/main.nf"
    workflow "HARMONISE"

    test("routes gwas-ssf to the hub") {
        when {
            workflow {
                """
                input[0] = Channel.of([
                    [ id:'demo', format:'gwas-ssf', n:10000, ancestry:'EUR', trait:'Demo', build:'GRCh37' ],
                    file("${projectDir}/assets/test/gwas_ssf_demo.tsv")
                ])
                input[1] = file("${projectDir}/assets/test/snp.info")
                """
            }
        }
        then {
            assert workflow.success
            assert workflow.out.harmonised.size() == 1
            def lines = path(workflow.out.harmonised[0][1]).linesGzip
            assert lines[0] == "rsid\tchr\tpos\tea\toa\teaf\tbeta\tse\tp\tn\tz"
        }
    }
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd "$DEST" && nf-test test tests/subworkflows/local/harmonise.nf.test`
Expected: FAIL — `subworkflows/local/harmonise/main.nf` does not exist.

- [ ] **Step 3: Write the subworkflow**

Create `$DEST/subworkflows/local/harmonise/main.nf`:

```groovy
//
// HARMONISE: branch raw sumstats by input format, normalise each to the
// canonical harmonised-sumstats hub (<id>.harmonised.tsv.gz).
// Phase 1 implements the gwas-ssf path; other formats fail loudly until ported.
//

include { HARMONISE_GWAS_SSF } from '../../../modules/local/harmonise_gwas_ssf'

workflow HARMONISE {

    take:
    ch_samplesheet   // channel: [ meta, sumstats ]
    snp_info         // value:   path to LD reference snp.info

    main:
    ch_versions = Channel.empty()

    def ch_branched = ch_samplesheet.branch {
        gwas_ssf: it[0].format == 'gwas-ssf'
        other:    true
    }

    // Phase 1 supports gwas-ssf only. Surface anything else as a hard error
    // instead of silently dropping the trait.
    ch_branched.other.map { meta, sumstats ->
        error "HARMONISE: format '${meta.format}' is not yet ported (Phase 1 = gwas-ssf only) for id='${meta.id}'"
    }

    HARMONISE_GWAS_SSF(ch_branched.gwas_ssf, snp_info)
    ch_versions = ch_versions.mix(HARMONISE_GWAS_SSF.out.versions)

    emit:
    harmonised = HARMONISE_GWAS_SSF.out.harmonised   // [ meta, harmonised.tsv.gz ]
    versions   = ch_versions                          // path versions.yml
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `cd "$DEST" && nf-test test tests/subworkflows/local/harmonise.nf.test`
Expected: PASS (1 test).

- [ ] **Step 5: Commit**

```bash
git -C "$DEST" add subworkflows/local/harmonise tests/subworkflows/local/harmonise.nf.test
git -C "$DEST" commit -m "feat: HARMONISE subworkflow (format branch; gwas-ssf wired, others guarded)"
```

---

### Task 10: Wire `workflows/sumstats.nf` + `main.nf`; end-to-end smoke

**Files:**
- Create: `$DEST/workflows/sumstats.nf`
- Create: `$DEST/main.nf` (copied from source, then edited)
- Test: `$DEST/tests/pipeline/main.nf.test`
- Create: `$DEST/nf-test.config` (copied verbatim)

- [ ] **Step 1: Copy nf-test config**

```bash
cp "$SRC/nf-test.config" "$DEST/nf-test.config"
```

- [ ] **Step 2: Write the Phase 1 `workflows/sumstats.nf`** (thin orchestrator: HARMONISE only)

Create `$DEST/workflows/sumstats.nf`:

```groovy
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { HARMONISE } from '../subworkflows/local/harmonise'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow SUMSTATS {

    take:
    ch_samplesheet // channel: [ meta, sumstats ]

    main:
    if (!params.ld_ref_dir) {
        error "params.ld_ref_dir must be set (select an SNP-set profile, e.g. -profile hm3, or pass --ld_ref_dir)"
    }
    def ch_snp_info = channel.fromPath("${params.ld_ref_dir}/snp.info", checkIfExists: true).first()

    HARMONISE(ch_samplesheet, ch_snp_info)

    emit:
    harmonised = HARMONISE.out.harmonised
    versions   = HARMONISE.out.versions
}
```

- [ ] **Step 3: Copy and adapt `main.nf`**

```bash
cp "$SRC/main.nf" "$DEST/main.nf"
```

Then in `$DEST/main.nf`:

(a) Update the include paths for the renamed init subworkflow. Find:

```groovy
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_sumstats_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_sumstats_pipeline'
```

Replace with:

```groovy
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_shared_sumstats_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_shared_sumstats_pipeline'
```

(b) Update the named wrapper workflow to drop the unused `outdir` arg. Find:

```groovy
workflow GEORGELAB_SUMSTATS {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SUMSTATS (
        samplesheet,
        params.outdir,
    )
}
```

Replace with:

```groovy
workflow SHARED_SUMSTATS {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    SUMSTATS (
        samplesheet
    )
}
```

(c) Update the call site in the entry `workflow {}`. Find:

```groovy
    GEORGELAB_SUMSTATS (
        PIPELINE_INITIALISATION.out.samplesheet
    )
```

Replace with:

```groovy
    SHARED_SUMSTATS (
        PIPELINE_INITIALISATION.out.samplesheet
    )
```

- [ ] **Step 4: Write the end-to-end pipeline nf-test**

Create `$DEST/tests/pipeline/main.nf.test`:

```groovy
nextflow_pipeline {

    name "Phase 1 end-to-end: samplesheet -> harmonised hub"
    script "main.nf"

    test("test profile produces one harmonised file") {
        when {
            params {
                input        = "${projectDir}/assets/test/samplesheet.csv"
                ld_ref_dir   = "${projectDir}/assets/test"
                outdir       = "${outputDir}"
                genome_build = "GRCh37"
                modules      = ""
            }
        }
        then {
            assert workflow.success
            // harmonise publishDir contains the demo hub file
            assert path("${outputDir}/harmonise/demo.harmonised.tsv.gz").exists()
        }
    }
}
```

> Note: this test uses `assets/test/samplesheet.csv`, whose `sumstats` path is relative. nf-test sets the launch dir to the project root, so the relative path resolves. If your nf-schema version resolves `file-path` relative to the samplesheet instead, change the CSV `sumstats` cell to an absolute path or `${projectDir}`-templated value generated in an nf-test `setup {}` block.

- [ ] **Step 5: Run the end-to-end test**

Run: `cd "$DEST" && nf-test test tests/pipeline/main.nf.test`
Expected: PASS (requires Docker + the container image, or local R with data.table).

- [ ] **Step 6: (Optional) manual smoke run**

Run: `cd "$DEST" && nextflow run . -profile test,docker --outdir results && zcat results/harmonise/demo.harmonised.tsv.gz | head`
Expected: pipeline completes; the harmonised file has the canonical header and 7 data rows.

- [ ] **Step 7: Commit**

```bash
git -C "$DEST" add main.nf workflows/sumstats.nf nf-test.config tests/pipeline/main.nf.test
git -C "$DEST" commit -m "feat: wire main.nf + workflows/sumstats.nf to run HARMONISE end-to-end"
```

---

### Task 11: README + usage stub for Phase 1

**Files:**
- Create: `$DEST/README.md`

- [ ] **Step 1: Write a minimal README** describing the current (Phase 1) state

Create `$DEST/README.md`:

```markdown
# bernooi/shared_sumstats

Shared, multi-contributor Nextflow / nf-core-style pipeline for descriptive
analysis of GWAS summary statistics. Ported from `bernooi/sumstats`.

**Design:** `docs/superpowers/specs/2026-06-03-shared-sumstats-design.md`
**Phase 1 plan:** `docs/superpowers/plans/2026-06-03-shared-sumstats-phase1-foundation.md`

## Status

- [x] Phase 1: scaffold + HARMONISE layer (gwas-ssf) -> harmonised hub
- [ ] Phase 2: gwama/hail/cvdkp harmonise paths
- [ ] Phase 3: SBayesRC module
- [ ] Phase 4: LDSC module
- [ ] Phase 5: report assembler (fragments)
- [ ] Phase 6: contribution scaffolding + CI

## Quick test (Phase 1)

```bash
nextflow run . -profile test,docker --outdir results
zcat results/harmonise/demo.harmonised.tsv.gz | head
```

## The harmonised hub

Every input format is normalised to one canonical artifact,
`<id>.harmonised.tsv.gz`, with columns:

`rsid chr pos ea oa eaf beta se p n z`

rsID is the canonical join key; `chr`/`pos` are attached from the LD reference
`snp.info` (record-only build, see the design spec).
```

- [ ] **Step 2: Commit**

```bash
git -C "$DEST" add README.md
git -C "$DEST" commit -m "docs: add README with Phase 1 status and quick test"
```

---

## Phase 1 Definition of Done

- [ ] `nf-test test` passes all three tests (process, subworkflow, pipeline).
- [ ] `nextflow run . -profile test,docker --outdir results` completes and writes `results/harmonise/demo.harmonised.tsv.gz` with the canonical header and 7 data rows.
- [ ] No maintainer-local absolute paths remain in committed configs (`grep -rn "/Volumes/" "$DEST" --include=*.config` returns nothing).
- [ ] Manifest points at `bernooi/shared_sumstats`.
- [ ] Run artifacts (`work/`, `results/`, `.nextflow*`) are git-ignored and untracked.

## Notes for later phases (not Phase 1 work)

- **Phase 2a (gwama/hail/cvdkp):** repeat Task 8/9 for each, basing the R scripts on `$SRC/bin/preprocess_{gwama,hail,cvdkp}.R` but emitting the hub columns (add `chr`/`pos` from `snp.info`, `z = beta/se`) instead of COJO `.ma`. Wire each into the HARMONISE branch and remove its guard. Add a fixture per format.
- **Phase 3 (SBayesRC):** new `subworkflows/local/sbayesrc/` consuming `ch_harmonised`; internal `harmonised -> .ma` converter; emits `results` + `report` fragment + `versions`. Port `SBAYESRC_TIDY/IMPUTE/MAIN` (keep `maxForks 1`).
- **Phase 5 (report):** `subworkflows/local/report/` collects `report` fragments via `groupTuple` and assembles per-trait HTML; removes the `NO_FILE` sentinel pattern from the source.
