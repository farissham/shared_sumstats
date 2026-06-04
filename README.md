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
