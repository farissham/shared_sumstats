//
// HARMONISE: normalise raw sumstats of any supported format to the canonical
// harmonised-sumstats hub (<id>.harmonised.tsv.gz). Cross-build chr:pos inputs
// with a configured chain are resolved by lifting the panel once per build
// (LIFTOVER_PANEL) and feeding the generated map to the harmoniser. Per-format
// parsing, the build matcher, and allele orientation live in bin/harmonise.R.
//

include { LIFTOVER_PANEL     } from '../../../modules/local/liftover_panel'
include { HARMONISE_SUMSTATS } from '../../../modules/local/harmonise'

workflow HARMONISE {

    take:
    ch_samplesheet   // channel: [ meta, sumstats ]
    snp_info         // value:   path to LD reference snp.info

    main:
    ch_versions = Channel.empty()

    def supported   = ['gwas-ssf', 'gwama', 'hail', 'cvdkp']
    def panel_build = params.ld_ref_build
    def chains      = params.liftover_chains ?: [:]

    def ch_in = ch_samplesheet.map { meta, sumstats ->
        if (!(meta.format in supported))
            error "HARMONISE: unsupported format '${meta.format}' for id='${meta.id}' (supported: ${supported.join(', ')})"
        [ meta, sumstats ]
    }

    // Route cross-build rows that have a configured chain (and no user-supplied
    // map/chain of their own) through a one-off panel liftover; everything else
    // (rsID inputs, same-build, or user-managed resolution) goes straight through.
    def ch_routed = ch_in.branch { meta, sumstats ->
        lift:   meta.build && meta.build != panel_build && !meta.rsid_map && !meta.chain && chains.containsKey(meta.build)
        direct: true
    }

    // Lift the panel once per distinct cross-build build -> [ build, rsid_map ].
    def ch_build_chain = ch_routed.lift
        .map { meta, sumstats -> meta.build }
        .unique()
        .map { build -> [ build, file(chains[build], checkIfExists: true) ] }

    LIFTOVER_PANEL(ch_build_chain, snp_info)
    ch_versions = ch_versions.mix(LIFTOVER_PANEL.out.versions)

    // Attach each generated map to its rows (join by build); direct rows carry none.
    def ch_lift = ch_routed.lift
        .map { meta, sumstats -> [ meta.build, meta, sumstats ] }
        .combine(LIFTOVER_PANEL.out.rsid_map, by: 0)
        .map { build, meta, sumstats, rmap -> [ meta, sumstats, rmap ] }

    def ch_direct = ch_routed.direct.map { meta, sumstats -> [ meta, sumstats, [] ] }

    HARMONISE_SUMSTATS(ch_lift.mix(ch_direct), snp_info)
    ch_versions = ch_versions.mix(HARMONISE_SUMSTATS.out.versions)

    emit:
    harmonised = HARMONISE_SUMSTATS.out.harmonised   // [ meta, harmonised.tsv.gz ]
    versions   = ch_versions                          // path versions.yml
}
