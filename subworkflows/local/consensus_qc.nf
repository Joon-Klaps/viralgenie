include { CHECKV_ENDTOEND } from '../../modules/nf-core/checkv/endtoend/main'
include { QUAST           } from '../../modules/nf-core/quast/main'

workflow CONSENSUS_QC  {

    take:
    ch_genome // channel: [ val(meta), [ genome ] ]
    checkv_db // channel: [ checkv_db ]
    skip_checkv // boolean
    skip_quast // boolean
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( !skip_checkv ) {
        CHECKV_ENDTOEND ( ch_genome, checkv_db )
        ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions)
    }

    if ( !skip_quast ) {
        QUAST (
            ch_genome,
            [],
            []
        )
        ch_versions = ch_versions.mix(QUAST.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv)
    }

    emit:
    mqc      = ch_multiqc_files // channel: [ tsv ]
    versions = ch_versions      // channel: [ versions.yml ]
}

