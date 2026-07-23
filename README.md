# farissham/shared_sumstats

Shared, multi-contributor Nextflow / nf-core-style pipeline for descriptive
analysis of GWAS summary statistics. Ported from `bernooi/sumstats`.

**Design:** `docs/superpowers/specs/2026-06-03-shared-sumstats-design.md`
**Phase 1 plan:** `docs/superpowers/plans/2026-06-03-shared-sumstats-phase1-foundation.md`

## Status

- [x] Phase 1: scaffold + HARMONISE layer (gwas-ssf) -> harmonised hub
- [x] Phase 2: harmonise paths — gwas-ssf, gwama, hail, cvdkp (unified `bin/harmonise.R`)
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

## Harmonisers, the build matcher, and orientation

All formats are parsed by per-format functions in `bin/harmonise.R`, which normalise
each layout to a common intermediate and then run the **same** shared steps:
resolve rsID (build matcher) -> panel join -> allele orientation -> hub.

| format | identifier | effect allele | notes |
|--------|-----------|---------------|-------|
| `gwas-ssf` | `rsid`/`variant_id` | `effect_allele` | per-SNP N (eff/total/n) |
| `gwama` | `rs_number` (rsID *or* chr:pos) | `reference_allele` | fixed N; `eaf` may be `-9` |
| `hail` | `locus` (chr:pos) + JSON `alleles` | ALT | always chr:pos -> build matcher |
| `cvdkp` | `rsid` (rsID rows kept; chr:pos resolved/dropped) | `Allele1` | `log10(p)`, per-SNP eff N |

Markers that are already rsIDs are used directly (build-agnostic). For `chr:pos`
markers the rsID join is build-sensitive, so the **build matcher** resolves it:

1. an explicit `rsid` column, if present;
2. else, `rs_number` values that look like rsIDs (`rs#######`) are used directly;
3. else (`chr:pos`), the **build matcher** maps positions to rsIDs by genome build:
   - if the input's `build` (samplesheet column, default `params.genome_build`)
     equals the panel build (`params.ld_ref_build`, default `GRCh37`), the map is
     derived directly from `snp.info` — no external file needed;
   - else, if the row supplies a UCSC `chain` (input build -> panel build), the
     `chr:pos` are **lifted over** to the panel build (inline, pure-R, no external
     tools) and then resolved from `snp.info`;
   - otherwise, supply a build-matched `chr,pos,rsid` map via the row's `rsid_map`.

If the matched fraction is implausibly low (default `< 1%`), the run fails with a
"declared build is likely wrong" error rather than silently dropping ~all variants.

Liftover runs only for cross-build `chr:pos` markers, immediately before the rsID
lookup; rsID and same-build markers skip it. Positions in chain gaps or on
chromosomes the chain doesn't cover are dropped. (`meta.chain` / `meta.rsid_map`
paths must be readable by the task container.)

There are two liftover options for cross-build `chr:pos` inputs:

- **inline** (`meta.chain`): a dependency-free pure-R chain applier in `harmonise.R`
  lifts the input's coordinates. Correct for the common single-block SNV case; does
  not replicate liftOver's multi-chain/`minMatch` edge-case handling.
- **`LIFTOVER_PANEL` module** (real UCSC `liftOver`, own container): lifts the
  *panel* `snp.info` once to a target build, producing a build-matched `chr,pos,rsid`
  map that the build matcher consumes via `rsid_map`. Parse-free, reusable across all
  inputs of that build, and only ever lifts clean panel SNVs (liftOver's ~99.99% case).
  Limited to panel SNPs (already the pipeline's scope) and to position-only matching.

The `LIFTOVER_PANEL` route is **auto-wired**: set `params.liftover_chains` to a
`build -> chain` map, e.g.

```groovy
liftover_chains = [ GRCh38: '/refs/chains/hg38ToHg19.over.chain.gz' ]
```

and any samplesheet row whose `build` differs from `ld_ref_build` (and that has no
`rsid_map`/`chain` of its own) triggers a one-off panel lift for that build, whose
map is fed to the harmoniser automatically. Builds with no configured chain fall
through to the normal resolution (rsID resolves; unmatched chr:pos errors loudly).

### Allele orientation (strand + palindromes)

Alleles are oriented to the panel (`ea = A1`, flipping `beta`/`eaf` as needed):

- **exact match** (same strand): `ra`/`oa` equal `A1`/`A2` in some order;
- **opposite strand**: the complement of `ra`/`oa` equals `A1`/`A2` — reverse-strand
  inputs are recovered instead of being silently dropped;
- **palindromic SNPs** (`A/T`, `C/G`): strand is unknowable from the alleles, so
  orientation is inferred from allele frequency (study `eaf` vs reference `A1Freq`).
  Palindromes are **dropped** when the MAF is too close to 0.5 to call (default
  `--palindrome_maf 0.42`) or when `eaf` is missing — pass `--assume_forward_strand`
  to keep them on the forward-strand assumption instead.

## Monitoring with Seqera Platform (Tower)

Live run monitoring is **opt-in** and activates automatically when a
`TOWER_ACCESS_TOKEN` is present in the environment — no config edits, no token
in git. Create a token at <https://cloud.seqera.io> (avatar → *Access tokens*),
then:

```bash
export TOWER_ACCESS_TOKEN=<your-token>
# export TOWER_WORKSPACE_ID=<id>   # optional: stream into a shared workspace
nextflow run . -profile test,docker --outdir results
```

Without `TOWER_ACCESS_TOKEN` set, runs behave exactly as before. Prefer not to
keep the token in your shell? Use Nextflow secrets instead:
`nextflow secrets set TOWER_ACCESS_TOKEN <token>`. Never commit the token.
