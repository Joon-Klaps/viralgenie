//
// Create contigs using
//
include { SPADES   } from '../../modules/nf-core/spades/main'
include { TRINITY  } from '../../modules/nf-core/trinity/main'
include { MEGAHIT  } from '../../modules/nf-core/megahit/main'
include { CAT_CAT  } from '../../modules/nf-core/cat/cat/main'

workflow FASTQ_SPADES_TRINITY_MEGAHIT  {

    take:
    reads           // channel: [ val(meta), [ reads ] ]
    assemblers      // value ['spades','trinity','megahit']
    ch_spades_yml   // channel: ['path/to/yml']
    ch_spades_hmm   // channel: ['path/to/hmm']

    main:
    ch_versions          = Channel.empty()
    ch_scaffolds         = Channel.empty()

    // SPADES
    if ('spades' in assemblers) {

        SPADES(
            reads.map {meta, reads -> [meta, reads, [], []]},
            ch_spades_yml,
            ch_spades_hmm
            )

        ch_versions         = ch_versions.mix(SPADES.out.versions.first())
        ch_scaffolds        = ch_scaffolds.mix( SPADES.out.scaffolds)
    }

    // TRINITY
    if ('trinity' in assemblers) {
        TRINITY(reads)

        TRINITY
            .out
            .transcript_fasta
            .filter{ meta, contigs -> contigs != null } // filter out empty contigs check issue #21
            .set{ch_scaffolds_trinity}

        ch_versions          = ch_versions.mix(TRINITY.out.versions.first())
        ch_scaffolds         = ch_scaffolds.mix(ch_scaffolds_trinity)
    }

    // MEGAHIT
    if ('megahit' in assemblers) {
        MEGAHIT(reads)
        ch_versions          = ch_versions.mix(MEGAHIT.out.versions.first())
        ch_scaffolds         = ch_scaffolds.mix(MEGAHIT.out.contigs)
    }

    // ch_scaffolds, go from [meta,scaffold1,scaffold2] to [meta,[scaffolds]]
    ch_scaffolds
        .groupTuple()
        .set{ch_scaffolds_combined}

    CAT_CAT(ch_scaffolds_combined)
    ch_versions = CAT_CAT.out.versions.first()



    emit:
    scaffolds            = CAT_CAT.out.file_out      // channel: [ val(meta), [ scaffolds] ]

    versions             = ch_versions               // channel: [ versions.yml ]
    // there are not any MQC files available for spades, trinity and megahit
}

