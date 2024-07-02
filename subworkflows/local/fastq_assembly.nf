//
// Create contigs using
//
include { SPADES                                     } from '../../modules/nf-core/spades/main'
include { TRINITY                                    } from '../../modules/nf-core/trinity/main'
include { MEGAHIT                                    } from '../../modules/nf-core/megahit/main'
include { FASTQ_FASTA_QUAST_SSPACE as EXTEND_SPADES  } from './fastq_fasta_quast_sspace.nf'
include { FASTQ_FASTA_QUAST_SSPACE as EXTEND_TRINITY } from './fastq_fasta_quast_sspace.nf'
include { FASTQ_FASTA_QUAST_SSPACE as EXTEND_MEGAHIT } from './fastq_fasta_quast_sspace.nf'
include { CAT_CAT as CAT_ASSEMBLERS                  } from '../../modules/nf-core/cat/cat/main'


workflow FASTQ_ASSEMBLY {

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
        ch_versions          = ch_versions.mix(SPADES.out.versions.first())

        EXTEND_SPADES( reads, SPADES.out.scaffolds, "spades")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_SPADES.out.scaffolds)
        ch_versions          = ch_versions.mix(EXTEND_SPADES.out.versions)
        ch_multiqc           = ch_multiqc.mix(EXTEND_SPADES.out.mqc)
    }


    // TRINITY
    if ('trinity' in assemblers) {
        TRINITY(reads)
        ch_versions          = ch_versions.mix(TRINITY.out.versions.first())

        EXTEND_TRINITY( reads, TRINITY.out.transcript_fasta, "trinity")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_TRINITY.out.scaffolds)
        ch_versions          = ch_versions.mix(EXTEND_TRINITY.out.versions)
        ch_multiqc           = ch_multiqc.mix(EXTEND_TRINITY.out.mqc)
    }

    // MEGAHIT
    if ('megahit' in assemblers) {
        MEGAHIT(reads)
        ch_versions          = ch_versions.mix(MEGAHIT.out.versions.first())

        EXTEND_MEGAHIT( reads, MEGAHIT.out.scaffolds, "megahit")
        ch_scaffolds         = ch_scaffolds.mix(EXTEND_MEGAHIT.out.scaffolds)
        ch_versions          = ch_versions.mix(EXTEND_MEGAHIT.out.versions)
        ch_multiqc           = ch_multiqc.mix(EXTEND_MEGAHIT.out.mqc)
    }

    // ch_scaffolds, go from [[meta,scaffold1],[meta,scaffold2], ...] to [meta,[scaffolds]]
    ch_scaffolds
        .map { meta, scaffold  -> tuple( groupKey(meta, assemblers.size()), scaffold ) }
        .groupTuple(remainder: true)
        .set{ch_scaffolds_combined}

    CAT_ASSEMBLERS(ch_scaffolds_combined)
    ch_scaffolds = CAT_ASSEMBLERS.out.file_out
    ch_versions  = CAT_ASSEMBLERS.out.versions.first()

    emit:
    scaffolds            = ch_scaffolds  // channel: [ val(meta), [ scaffolds] ]
    mqc                  = ch_multiqc    // channel: [ val(meta), [ mqc ] ]
    versions             = ch_versions   // channel: [ versions.yml ]
    // there are not any MQC files available for spades, trinity and megahit
}

