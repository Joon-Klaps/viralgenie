//
// Identify the best hits of a set from a set of reference genomes using MASH
//

include { FIND_CONCATENATE as CAT_CAT_READS } from '../../../modules/nf-core/find/concatenate/main'
include { MASH_SKETCH              } from '../../../modules/nf-core/mash/sketch/main'
include { MASH_SCREEN              } from '../../../modules/nf-core/mash/screen/main'
include { SELECT_REFERENCE         } from '../../../modules/local/select_reference/main'

workflow FASTQ_FASTA_MASH_SCREEN {
    take:
    ch_fasta_reads // channel of [[meta], [multi-fasta], [read1, read2]]

    main:

    //
    // Join reads
    //
    ch_input_cat = ch_fasta_reads.map { meta, _fasta, reads -> [meta, reads] }
    CAT_CAT_READS(ch_input_cat)


    //
    // Sketch the input sequences
    //
    ch_input_sketch = ch_fasta_reads.map { meta, fasta, _reads -> [meta, fasta] }
    MASH_SKETCH(ch_input_sketch)

    ch_input_screen = CAT_CAT_READS.out.file_out
        .join(MASH_SKETCH.out.mash)
        .multiMap { meta, fastq, sketch ->
            query: [meta, fastq]
            sequences: [meta, sketch]
        }

    //
    // Identify best hits from the reference sketch.
    //
    MASH_SCREEN(ch_input_screen.query, ch_input_screen.sequences)

    //
    // Isolate/extract the best hit from mash screen using custom script select_reference
    //
    ch_input_select_reference = MASH_SCREEN.out.screen.join(ch_fasta_reads)
    SELECT_REFERENCE(ch_input_select_reference)

    ch_reference_fastq = SELECT_REFERENCE.out.fasta_reads
        .filter { _meta, _json, fasta, _reads ->
            fasta.countFasta() > 0
        }
        .map { meta, _json, fasta, reads ->
            [meta, fasta, reads]
        }

    ch_json = SELECT_REFERENCE.out.fasta_reads.map { meta, json, _fasta, _reads -> [meta, json] }

    emit:
    reference_fastq = ch_reference_fastq // channel: [meta, fasta, reads ]
    json            = ch_json            // channel: [meta, json ]
}
