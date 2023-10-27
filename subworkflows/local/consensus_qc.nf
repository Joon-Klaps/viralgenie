include { CHECKV_ENDTOEND                 } from '../../modules/nf-core/checkv/endtoend/main'
include { CAT_CAT as CAT_CAT_QC           } from '../../modules/nf-core/cat/cat/main'
include { QUAST  as QUAST_QC              } from '../../modules/nf-core/quast/main'
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

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( !skip_checkv ) {
        ch_genome
            .map{meta, genome -> [meta.subMap('id','cluster_id','sample'), genome]}
            .groupTuple()
            .set{ch_genome_collapsed}

        CAT_CAT_QC(ch_genome_collapsed)
        ch_versions = ch_versions.mix(CAT_CAT_QC.out.versions)

        CHECKV_ENDTOEND ( CAT_CAT_QC.out.file_out, checkv_db )
        ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions)

        //insert alignment of everything here
    }

    if ( !skip_quast ) {
        QUAST_QC (
            ch_genome,
            [[:],[]],
            [[:],[]]
        )
        ch_versions = ch_versions.mix(QUAST_QC.out.versions)
        // Will keep this for now but this should be handled by another process that makes it into a nice table for the multiqc report
        ch_multiqc_files = ch_multiqc_files.mix(QUAST_QC.out.tsv)
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

