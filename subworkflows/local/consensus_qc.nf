include { CHECKV_ENDTOEND              } from '../../modules/nf-core/checkv/endtoend/main'
include { CAT_CAT as CAT_CAT_CONSENSUS } from '../../modules/nf-core/cat/cat/main'
include { QUAST                        } from '../../modules/nf-core/quast/main'

workflow CONSENSUS_QC  {

    take:
    ch_genome   // channel: [ val(meta), [ genome ] ]
    checkv_db   // channel: [ checkv_db ]
    skip_checkv // boolean
    skip_quast  // boolean
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //TODO: Check if we can distinguish which genome came from which cluster group, sample, ...
    ch_genome
        .collect{it[1]}
        .map{it -> [[id: 'genomes_combined'], it]}
        .set { ch_genome_combined }
    CAT_CAT_CONSENSUS(ch_genome_combined)
    ch_versions = ch_versions.mix(CAT_CAT_CONSENSUS.out.versions)

    if ( !skip_checkv ) {
        CHECKV_ENDTOEND ( CAT_CAT_CONSENSUS.out.file_out, checkv_db )
        ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions)
    }

    if ( !skip_quast ) {
        QUAST (
            CAT_CAT_CONSENSUS.out.file_out,
            [[:],[]],
            [[:],[]]
        )
        ch_versions = ch_versions.mix(QUAST.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv)
    }

    emit:
    mqc      = ch_multiqc_files // channel: [ tsv ]
    versions = ch_versions      // channel: [ versions.yml ]
}

