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
    refpool_db             // channel: [ val(meta), [refpool_db], [refpool_seq] ]
    annotation_db          // channel: [ val(meta), [annotation_db] ]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    blast_txt        = Channel.empty()
    checkv_summary   = Channel.empty()
    quast_summary    = Channel.empty()
    annotation_txt  = Channel.empty()

    if ( !params.skip_checkv || !params.skip_alignment_qc) {
        ch_genome
            .map{meta, genome -> [meta.subMap('id','cluster_id','sample'), genome]}
            .groupTuple()
            .set{ch_genome_grouped}

        CAT_CAT_QC(
            ch_genome_grouped
        )

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
        CHECKV_ENDTOEND (
            ch_genome_collapsed,
            checkv_db
        )
        checkv_summary = CHECKV_ENDTOEND.out.quality_summary
        ch_versions    = ch_versions.mix(CHECKV_ENDTOEND.out.versions)
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

        MAFFT_ITERATIONS (
            ch_genome_collapsed_branch.pass,
            [[:],[]],
            [[:],[]],
            [[:],[]],
            [[:],[]],
            [[:],[]],
            false
        )
        ch_versions = ch_versions.mix(MAFFT_ITERATIONS.out.versions)

        // Mix single sequences with the (multiple) aligned ones
        ch_genome_collapsed_branch
            .fail
            .mix(MAFFT_ITERATIONS.out.fas)
            .map{ meta, genome -> [meta.id, meta, genome] }
            .set{ch_genome_collapsed_mod}

        ch_aligned_raw_contigs.map{ meta, genome -> [meta.id, meta, genome] }.set{ch_aligned_raw_contigs_mod}

        // Combine with the orignal contigs
        ch_genome_collapsed_mod
            .join( ch_aligned_raw_contigs_mod, by: 0)
            .map{ id, meta_genome, scaffolds, meta_contigs, contigs -> [meta_genome, scaffolds, contigs] }
            .filter{ meta, scaffolds, contigs ->
                scaffolds.countFasta() + contigs.countFasta() > 1
            }
            .set{ch_genome_collapsed_branch}

        ch_genome_collapsed_branch
            .map{ meta, scaffolds, contigs -> [meta, scaffolds] }
            .set{scaffolds}
        ch_genome_collapsed_branch
            .map{ meta, scaffolds, contigs -> [meta,contigs] }
            .set{addsequences}

        MAFFT_QC (
            scaffolds,
            addsequences,
            [[:],[]],
            [[:],[]],
            [[:],[]],
            [[:],[]],
            false
        )
        ch_versions = ch_versions.mix(MAFFT_QC.out.versions)
    }

    if ( !params.skip_quast ) {
        // Contig summary statistics
        QUAST_QC (
            ch_genome,
            [[:],[]],
            [[:],[]]
        )
        ch_versions   = ch_versions.mix(QUAST_QC.out.versions)
        quast_summary = QUAST_QC.out.tsv
    }

    if ( !params.skip_blast_qc ){

        // combine refpool_db based on specified samples in the reference_pools parameter
        ch_genome
            .combine(refpool_db)
            .filter{ meta_genome, genome, meta_db, blast_db, blast_seq ->
                meta_genome.sample == meta_db.sample || (meta_genome.sample != null && meta_db.sample == null)}
            .branch{ meta_genome, genome, meta_db, blast_db, blast_seq ->
                genome: [meta_genome, genome]
                db: [meta_db, blast_db]
            }
            .set{ch_blast_in}

        // Identify closest reference from the reference pool database using blast
        BLASTN_QC (
            ch_blast_in.genome,
            ch_blast_in.db
        )
        blast_txt   = BLASTN_QC.out.txt
        ch_versions = ch_versions.mix(BLASTN_QC.out.versions)
    }

    if ( !params.skip_annotation){
        ch_genomes_collect = ch_genome.collect{it[1]}.map{files -> [[id:"all_genomes_annotation.hits"], files]}
        CAT_CAT_MMSEQS(
            ch_genomes_collect
        )
        ch_versions = ch_versions.mix(CAT_CAT_MMSEQS.out.versions)
        // use MMSEQS easy search to find best hits against annotation db
        MMSEQS_ANNOTATE(
            CAT_CAT_MMSEQS.out.file_out,
            annotation_db
        )
        annotation_txt = MMSEQS_ANNOTATE.out.tsv
        ch_versions = ch_versions.mix(MMSEQS_ANNOTATE.out.versions)
    }


    emit:
    blast_txt       = blast_txt         // channel: [ val(meta), [ txt ] ]
    checkv_summary  = checkv_summary    // channel: [ val(meta), [ tsv ] ]
    quast_summary   = quast_summary     // channel: [ val(meta), [ tsv ] ]
    annotation_txt  = annotation_txt   // channel: [ val(meta), [ txt ] ]
    mqc             = ch_multiqc_files  // channel: [ tsv ]
    versions        = ch_versions       // channel: [ versions.yml ]
}

