include { CHECKV_ENDTOEND                 } from '../../modules/nf-core/checkv/endtoend/main'
include { CAT_CAT as CAT_CAT_QC           } from '../../modules/nf-core/cat/cat/main'
include { QUAST  as QUAST_QC              } from '../../modules/nf-core/quast/main'
include { BLAST_BLASTN as BLAST_BLASTN_QC } from '../../modules/nf-core/blast/blastn/main'
include { MAFFT as MAFFT_QC               } from '../../modules/nf-core/mafft/main'

workflow CONSENSUS_QC  {

    take:
    ch_genome         // channel: [ val(meta), [ genome ] ]
    checkv_db         // channel: [ checkv_db ]
    blast_db          // channel: [ val(meta), [blast_db] ]
    skip_checkv       // boolean
    skip_quast        // boolean
    skip_blast_qc     // boolean
    skip_alignment_qc // boolean
    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    blast_txt        = Channel.empty()
    checkv_summary   = Channel.empty()
    quast_summary    = Channel.empty()

    if ( !skip_checkv || !skip_alignment_qc) {
        ch_genome
            .map{meta, genome -> [meta.subMap('id','cluster_id','sample'), genome]}
            .groupTuple()
            .set{ch_genome_grouped}

        CAT_CAT_QC(ch_genome_grouped)

        CAT_CAT_QC
            .out
            .file_out
            .set{ch_genome_collapsed}
        ch_versions = ch_versions.mix(CAT_CAT_QC.out.versions)
    }

    if ( !skip_checkv ) {
        // uses HMM and AA alignment to deterimine completeness
        CHECKV_ENDTOEND ( ch_genome_collapsed, checkv_db )
        checkv_summary = CHECKV_ENDTOEND.out.quality_summary
        ch_versions    = ch_versions.mix(CHECKV_ENDTOEND.out.versions)
    }

    if ( !skip_alignment_qc){
        // Align the different steps to each other to see how the sequences have changed
        // MAFFT doesn't like those that have only one sequence
        ch_genome_collapsed
            .branch { meta, scaffolds ->
                pass: scaffolds.countFasta() > 1
                fail: scaffolds.countFasta() == 1
            }
            .set{ch_genome_collapsed_branch}
        MAFFT_QC (
            ch_genome_collapsed_branch.pass,
            [],
        )
        ch_versions = ch_versions.mix(MAFFT_QC.out.versions)
    }

    if ( !skip_quast ) {
        // Basic summary statistics
        QUAST_QC (
            ch_genome,
            [[:],[]],
            [[:],[]]
        )
        ch_versions   = ch_versions.mix(QUAST_QC.out.versions)
        quast_summary = QUAST_QC.out.tsv
    }

    if ( !skip_blast_qc ){
        // Identify closest reference from the database
        BLAST_BLASTN_QC (
            ch_genome,
            blast_db
        )
        ch_versions = ch_versions.mix(BLAST_BLASTN_QC.out.versions)
        blast_txt   = BLAST_BLASTN_QC.out.txt
    }


    emit:
    blast_txt       = blast_txt         // channel: [ val(meta), [ txt ] ]
    checkv_summary  = checkv_summary    // channel: [ val(meta), [ tsv ] ]
    quast_summary   = quast_summary     // channel: [ val(meta), [ tsv ] ]
    mqc             = ch_multiqc_files  // channel: [ tsv ]
    versions        = ch_versions       // channel: [ versions.yml ]
}

