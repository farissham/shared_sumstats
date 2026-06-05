# bernooi/shared_sumstats

Shared, multi-contributor Nextflow / nf-core-style pipeline for descriptive
analysis of GWAS summary statistics. Ported from `bernooi/sumstats`.

**Design:** `docs/superpowers/specs/2026-06-03-shared-sumstats-design.md`
**Phase 1 plan:** `docs/superpowers/plans/2026-06-03-shared-sumstats-phase1-foundation.md`

## Status

- [x] Phase 1: scaffold + HARMONISE layer (gwas-ssf) -> harmonised hub
- [~] Phase 2: harmonise paths — gwama done; hail/cvdkp pending
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

## GWAMA inputs and the build matcher

GWAMA's `rs_number` column may hold actual rsIDs *or* a `chr:pos` string, depending
on the input studies. The rsID join to the panel is build-agnostic, but the
`chr:pos -> rsID` step is build-sensitive, so the gwama path resolves it as follows:

1. an explicit `rsid` column, if present;
2. else, `rs_number` values that look like rsIDs (`rs#######`) are used directly;
3. else (`chr:pos`), the **build matcher** maps positions to rsIDs by genome build:
   - if the input's `build` (samplesheet column, default `params.genome_build`)
     equals the panel build (`params.ld_ref_build`, default `GRCh37`), the map is
     derived directly from `snp.info` — no external file needed;
   - otherwise, supply a build-matched `chr,pos,rsid` map via the row's `rsid_map`
     (or liftover the input to the panel build first).

If the matched fraction is implausibly low (default `< 1%`), the run fails with a
"declared build is likely wrong" error rather than silently dropping ~all variants.

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
