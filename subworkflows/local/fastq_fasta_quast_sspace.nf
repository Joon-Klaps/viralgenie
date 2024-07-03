//
// Create contigs using
//
include { QUAST         } from '../../modules/nf-core/quast/main'
include { SSPACE_BASIC  } from '../../modules/local/sspace_basic/main'

workflow FASTQ_FASTA_QUAST_SSPACE {

    take:
    reads       // channel: [ val(meta), [ reads ] ]
    scaffolds   // channel: [ val(meta), [ scaffolds ] ]
    name        // value 'spades','trinity','megahit'

    main:
    ch_versions   = Channel.empty()
    ch_scaffolds  = Channel.empty()
    ch_multiqc    = Channel.empty()

    scaffolds
        .filter(meta, contigs -> contigs != null)
        .filter(meta, contigs -> contigs.countFasta() > 0)
        .set{ch_scaffolds}

    // QUAST
    QUAST(ch_scaffolds, [[:],[]], [[:],[]])
    ch_versions = ch_versions.mix(QUAST.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(QUAST.out.tsv)

    // SSPACE_BASIC
    if (!params.skip_sspace_basic){
        ch_scaffolds_filter
            .join(reads)
            .multiMap { meta, scaffolds, reads ->
                reads : [meta, reads]
                scaffolds : [meta, scaffolds]
                settings: [params.read_distance, params.read_distance_sd, params.read_orientation]
                name: [name]
            }
            .set{ch_sspace_input}

        SSPACE_BASIC(
            ch_sspace_input.reads,
            ch_sspace_input.scaffolds,
            ch_sspace_input.settings,
            ch_sspace_input.name
        )
        ch_versions = ch_versions.mix(SSPACE_BASIC.out.versions.first())

        ch_scaffolds = SSPACE_BASIC.out.fasta
    }


    emit:
    scaffolds            = ch_scaffolds  // channel: [ val(meta), [ scaffolds] ]
    mqc                  = ch_multiqc    // channel: [ val(meta), [ mqc ] ]
    versions             = ch_versions   // channel: [ versions.yml ]
}

