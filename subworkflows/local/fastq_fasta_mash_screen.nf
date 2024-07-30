//
// Identify the best hits of a set from a set of reference genomes using MASH
//

include { CAT_CAT as CAT_CAT_READS } from '../../modules/nf-core/cat/cat/main'
include { MASH_SKETCH              } from '../../modules/nf-core/mash/sketch/main'
include { MASH_SCREEN              } from '../../modules/nf-core/mash/screen/main'
include { SELECT_REFERENCE         } from '../../modules/local/select_reference/main'

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

    //
    // Isolate/extract the best hit from mash screen using custom script select_reference
    //
    ch_input_select_reference = MASH_SCREEN.out.screen.join(fasta_reads)
    SELECT_REFERENCE ( ch_input_select_reference )
    ch_versions = ch_versions.mix(SELECT_REFERENCE.out.versions)

    reference_fastq = SELECT_REFERENCE.out.fasta_reads
        .filter{
            meta, json, fasta, reads ->
            fasta.countFasta() > 0
        }
        .map{
            meta, json, fasta, reads ->
            json = WorkflowCommons.getMapFromJson(json)
            return [meta + json, fasta, reads]
        }

    emit:
    reference_fastq = reference_fastq   // channel: [meta, fasta, reads]
    versions        = ch_versions       // channel: [ versions.yml ]
}

