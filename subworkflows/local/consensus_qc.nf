include { CHECKV_DOWNLOADDATABASE           } from '../../modules/nf-core/checkv/downloaddatabase/main'
include { CHECKV_ENDTOEND                   } from '../../modules/nf-core/checkv/endtoend/main'
include { CAT_CAT as CAT_CAT_QC             } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as CAT_CAT_MMSEQS         } from '../../modules/nf-core/cat/cat/main'
include { QUAST  as QUAST_QC                } from '../../modules/nf-core/quast/main'
include { BLAST_BLASTN as BLASTN_QC         } from '../../modules/nf-core/blast/blastn/main'
include { MAFFT as MAFFT_ITERATIONS         } from '../../modules/nf-core/mafft/main'
include { MAFFT as MAFFT_QC                 } from '../../modules/nf-core/mafft/main'
include { MMSEQS_ANNOTATE                   } from './mmseqs_annotate.nf'

workflow CONSENSUS_QC  {

    take:
    ch_genome              // channel: [ val(meta), [ genome ] ]
    ch_aligned_raw_contigs // channel: [ val(meta), [ genome ] ]
    checkv_db              // channel: [ checkv_db ]
    refpool_db             // channel: [ val(meta), [refpool_db] ]
    annotation_db          // channel: [ val(meta), [annotation_db] ]s

    main:

    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()
    blast               = Channel.empty()
    checkv              = Channel.empty()
    quast               = Channel.empty()
    annotation          = Channel.empty()
    ch_genome_collapsed = Channel.empty()

    if ( !params.skip_checkv || !params.skip_alignment_qc) {
        ch_genome
            .map{meta, genome -> [meta.subMap('id','cluster_id','sample'), genome]}
            .groupTuple()
            .set{catcat_in}

        CAT_CAT_QC( catcat_in )

        CAT_CAT_QC
            .out
            .file_out
            .set{ch_genome_collapsed}
        ch_versions = ch_versions.mix(CAT_CAT_QC.out.versions)
    }

    if ( !params.skip_checkv ) {
        if ( !params.checkv_db ) {
            CHECKV_DOWNLOADDATABASE()
            checkv_db   = CHECKV_DOWNLOADDATABASE.out.checkv_db
            ch_versions = ch_versions.mix(CHECKV_DOWNLOADDATABASE.out.versions)
        }

        // uses HMM and AA alignment to deterimine completeness
        CHECKV_ENDTOEND ( ch_genome_collapsed, checkv_db)

        checkv      = CHECKV_ENDTOEND.out.quality_summary
        ch_versions = ch_versions.mix(CHECKV_ENDTOEND.out.versions)
    }

    // Align the different steps to each other to see how the sequences have changed
    if ( !params.skip_alignment_qc){

        // MAFFT doesn't like those that have only one sequence
        ch_genome_collapsed
            .branch { meta, scaffolds ->
                pass: scaffolds.countFasta() > 1
                fail: scaffolds.countFasta() == 1
            }
            .set{ch_genome_collapsed_branch}

        MAFFT_ITERATIONS ( ch_genome_collapsed_branch.pass, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]], false )

        ch_versions = ch_versions.mix(MAFFT_ITERATIONS.out.versions)
        contigs_mod = ch_aligned_raw_contigs.map{ meta, genome -> [meta.id, meta, genome] }

        // Make a channel that contains the alignment of the iterations with
        // the original contigs from the assemblers
        ch_genome_collapsed_branch
            .fail.mix(MAFFT_ITERATIONS.out.fas)                             // Combine alignments with single results
            .map{ meta, genome -> [meta.id, meta, genome] }                 // Set common delimiter
            .join( contigs_mod, by: 0)                                      // Combine with raw contigs
            .filter{id, meta_genome, scaffolds, meta_contigs, contigs  ->   // Make sure we have at least 2 sequences
                scaffolds.countFasta() + contigs.countFasta() > 1
            }
            .multiMap{ id, meta_genome, scaffolds, meta_contigs, contigs -> // Split in correct inputs
                scaffolds: [meta_genome, scaffolds]
                contigs: [meta_genome, contigs]
            }.set{mafftQC_in}

        MAFFT_QC ( mafftQC_in.scaffolds, mafftQC_in.contigs, [[:],[]], [[:],[]], [[:],[]], [[:],[]], false)

        ch_versions = ch_versions.mix(MAFFT_QC.out.versions)
    }

    // Contig summary statistics
    if ( !params.skip_quast ) {
        QUAST_QC ( ch_genome, [[:],[]], [[:],[]])

        ch_versions   = ch_versions.mix(QUAST_QC.out.versions)
        quast = QUAST_QC.out.tsv
    }

    // Identify closest reference from the reference pool database using blast
    if ( !params.skip_blast_qc ){
        BLASTN_QC ( ch_genome, refpool_db)

        blast       = BLASTN_QC.out.txt
        ch_versions = ch_versions.mix(BLASTN_QC.out.versions)
    }

    // use MMSEQS easy search to find best hits against annotation db
    if ( !params.skip_annotation){
        ch_genomes_collect = ch_genome.collect{it[1]}.map{files -> [[id:"all_genomes_annotation.hits"], files]}
        CAT_CAT_MMSEQS( ch_genomes_collect )
        ch_versions = ch_versions.mix(CAT_CAT_MMSEQS.out.versions)

        MMSEQS_ANNOTATE(CAT_CAT_MMSEQS.out.file_out,annotation_db)

        annotation  = MMSEQS_ANNOTATE.out.tsv
        ch_versions = ch_versions.mix(MMSEQS_ANNOTATE.out.versions)
    }

    emit:
    blast       = blast             // channel: [ val(meta), [ txt ] ]
    checkv      = checkv            // channel: [ val(meta), [ tsv ] ]
    quast       = quast             // channel: [ val(meta), [ tsv ] ]
    annotation  = annotation        // channel: [ val(meta), [ txt ] ]
    mqc         = ch_multiqc_files  // channel: [ tsv ]
    versions    = ch_versions       // channel: [ versions.yml ]
}

