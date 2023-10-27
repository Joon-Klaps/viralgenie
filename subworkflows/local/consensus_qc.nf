include { CHECKV_ENDTOEND                 } from '../../modules/nf-core/checkv/endtoend/main'
include { CAT_CAT as CAT_CAT_QC           } from '../../modules/nf-core/cat/cat/main'
include { QUAST                           } from '../../modules/nf-core/quast/main'
include { BLAST_BLASTN as BLAST_BLASTN_QC } from '../../modules/nf-core/blast/blastn/main'

workflow CONSENSUS_QC  {

    take:
    ch_genome      // channel: [ val(meta), [ genome ] ]
    checkv_db      // channel: [ checkv_db ]
    blast_db       // channel: [ blast_db ]
    skip_checkv    // boolean
    skip_quast     // boolean
    skip_blast_qc  // boolean
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( !skip_checkv ) {
        ch_genome
            .collect{it[1]}
            .map{it -> [ [id: "final_seq_combined"], it ] }
            .set{ch_genome_collapsed}
        CAT_CAT_QC(ch_genome_collapsed)
        ch_versions = ch_versions.mix(CAT_CAT_QC.out.versions)

        CHECKV_ENDTOEND ( CAT_CAT_QC.out.file_out, checkv_db )
        ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions)
    }

    if ( !skip_quast ) {
        QUAST (
            ch_genome,
            [[:],[]],
            [[:],[]]
        )
        ch_versions = ch_versions.mix(QUAST.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv)
    }

    if ( !skip_blast_qc ){
        BLAST_BLASTN_QC (
            ch_genome,
            blast_db
        )
    }


    emit:
    blast_txt = BLAST_BLASTN_QC.out.txt // channel: [ val(meta), [ txt ] ]
    mqc       = ch_multiqc_files        // channel: [ tsv ]
    versions  = ch_versions             // channel: [ versions.yml ]
}

