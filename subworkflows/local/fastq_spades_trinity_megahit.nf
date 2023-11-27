//
// Create contigs using
//
include { SPADES                 } from '../../modules/nf-core/spades/main'
include { QUAST as QUAST_SPADES  } from '../../modules/nf-core/quast/main'
include { TRINITY                } from '../../modules/nf-core/trinity/main'
include { QUAST as QUAST_TRINITY } from '../../modules/nf-core/quast/main'
include { MEGAHIT                } from '../../modules/nf-core/megahit/main'
include { QUAST as QUAST_MEGAHIT } from '../../modules/nf-core/quast/main'
include { CAT_CAT                } from '../../modules/nf-core/cat/cat/main'

workflow FASTQ_SPADES_TRINITY_MEGAHIT  {

    take:
    reads           // channel: [ val(meta), [ reads ] ]
    assemblers      // value ['spades','trinity','megahit']
    ch_spades_yml   // channel: ['path/to/yml']
    ch_spades_hmm   // channel: ['path/to/hmm']

    main:
    ch_versions   = Channel.empty()
    ch_scaffolds  = Channel.empty()
    ch_multiqc    = Channel.empty()

    // SPADES
    if ('spades' in assemblers) {

        SPADES(
            reads.map {meta, reads -> [meta, reads, [], []]},
            ch_spades_yml,
            ch_spades_hmm
            )

        ch_versions         = ch_versions.mix(SPADES.out.versions.first())
        ch_scaffolds        = ch_scaffolds.mix( SPADES.out.scaffolds)
        spades_filtered     = SPADES.out.scaffolds.filter{ meta, contigs -> contigs.countFasta() > 0 }

        QUAST_SPADES (
            SPADES.out.scaffolds.filter(),
            [[:],[]],
            [[:],[]]
        )
        ch_versions         = ch_versions.mix(QUAST_SPADES.out.versions.first())
        ch_multiqc          = ch_multiqc.mix(QUAST_SPADES.out.tsv)
    }

    // TRINITY
    if ('trinity' in assemblers) {
        TRINITY(reads)

        TRINITY
            .out
            .transcript_fasta
            .filter{ meta, contigs -> contigs != null } // filter out empty contigs check issue #21
            .set{trinity_filtered}

        ch_scaffolds         = ch_scaffolds.mix(trinity_filtered)
        ch_versions          = ch_versions.mix(TRINITY.out.versions.first())

        QUAST_TRINITY (
            trinity_filtered,
            [[:],[]],
            [[:],[]]
        )

        ch_versions          = ch_versions.mix(QUAST_TRINITY.out.versions.first())
        ch_multiqc           = ch_multiqc.mix(QUAST_TRINITY.out.tsv)
    }

    // MEGAHIT
    if ('megahit' in assemblers) {
        MEGAHIT(reads)
        ch_versions          = ch_versions.mix(MEGAHIT.out.versions.first())
        ch_scaffolds         = ch_scaffolds.mix(MEGAHIT.out.contigs)
        megahit_filtered     = MEGAHIT.out.contigs.filter{ meta, contigs -> contigs.countFasta() > 0 }


        QUAST_MEGAHIT (
            megahit_filtered,
            [[:],[]],
            [[:],[]]
        )
        ch_versions          = ch_versions.mix(QUAST_MEGAHIT.out.versions.first())
        ch_multiqc           = ch_multiqc.mix(QUAST_MEGAHIT.out.tsv)
    }

    // ch_scaffolds, go from [[meta,scaffold1],[meta,scaffold2], ...] to [meta,[scaffolds]]
    ch_scaffolds
        .map { meta, scaffold  -> tuple( groupKey(meta, assemblers.size()), scaffold ) }
        .groupTuple(remainder: true)
        .set{ch_scaffolds_combined}

    CAT_CAT(ch_scaffolds_combined)
    ch_versions = CAT_CAT.out.versions.first()



    emit:
    scaffolds            = CAT_CAT.out.file_out      // channel: [ val(meta), [ scaffolds] ]
    mqc                  = ch_multiqc                // channel: [ val(meta), [ mqc ] ]
    versions             = ch_versions               // channel: [ versions.yml ]
    // there are not any MQC files available for spades, trinity and megahit
}

