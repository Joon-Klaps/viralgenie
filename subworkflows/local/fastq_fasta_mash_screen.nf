//
// Identify the best hits of a set from a set of reference genomes using MASH
//

include { CAT_CAT as CAT_CAT_READS } from '../../modules/nf-core/cat/cat/main'
include { MASH_SKETCH              } from '../../modules/nf-core/mash/sketch/main'
include { MASH_SCREEN              } from '../../modules/nf-core/mash/screen/main'

workflow FASTQ_FASTA_MASH_SCREEN {

    take:
    fasta_reads // channel of [[meta], [multi-fasta], [read1, read2]]

    main:
    ch_versions = Channel.empty()

    //
    // Join reads
    //
    ch_input_cat = fasta_reads.map{meta, fasta, reads -> [meta, reads]}
    CAT_CAT_READS ( ch_input_cat )
    ch_versions = ch_versions.mix(CAT_CAT_READS.out.versions)


    //
    // Sketch the input sequences
    //
    ch_input_sketch = fasta_reads.map{meta, fasta, reads -> [meta, fasta]}
    MASH_SKETCH ( ch_input_sketch )
    ch_versions = ch_versions.mix(MASH_SKETCH.out.versions)

    ch_input_screen = CAT_CAT_READS.out.file_out
        .join(MASH_SKETCH.out.mash)
        .multiMap{
            meta, fastq, sketch ->
            query: [meta, fastq]
            sequences: [meta, sketch]
        }

    //
    // Identify best hits from the reference sketch.
    //
    MASH_SCREEN ( ch_input_screen.query, ch_input_screen.sequences )
    ch_versions = ch_versions.mix(MASH_SCREEN.out.versions)


    // TODO: write a script that extracts the top hit from the screen output
    out = MASH_SCREEN.out.screen

    emit:
    out      = out
    versions = ch_versions                     // channel: [ versions.yml ]
}

