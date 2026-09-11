include { CHECKV_DOWNLOADDATABASE         } from '../../../modules/nf-core/checkv/downloaddatabase/main'
include { CHECKV_ENDTOEND                 } from '../../../modules/nf-core/checkv/endtoend/main'
include { QUAST as QUAST_QC               } from '../../../modules/nf-core/quast/main'
include { BLAST_BLASTN as BLASTN_QC       } from '../../../modules/nf-core/blast/blastn/main'
include { MAFFT_ALIGN as MAFFT_ITERATIONS } from '../../../modules/nf-core/mafft/align/main'
include { MAFFT_ALIGN as MAFFT_QC         } from '../../../modules/nf-core/mafft/align/main'
include { PROKKA                          } from '../../../modules/nf-core/prokka/main'
include { MMSEQS_ANNOTATE                 } from '../mmseqs_annotate'

workflow CONSENSUS_QC {
    take:
    ch_genome              // channel: [ val(meta), [ genome ] ]
    ch_aligned_raw_contigs // channel: [ val(meta), [ genome ] ]
    ch_checkv_db           // channel: [ checkv_db ]
    ch_refpool_db          // channel: [ val(meta), [refpool_db] ]
    ch_annotation_db       // channel: [ val(meta), [annotation_db] ]
    ch_prokka_db           // channel: [ val(meta), [prokka_db] ]

    main:

    ch_blast           = channel.empty()
    ch_checkv          = channel.empty()
    ch_quast           = channel.empty()
    ch_annotation      = channel.empty()
    ch_genome_grouped  = channel.empty()

    // Combine all genomes into a single file
    ch_genomes_all = ch_genome
        .collectFile(name: "all_genomes.fa") {_meta, fasta -> fasta }
        .map { fasta -> [[id: "all_genomes"], fasta] }

    // combine the different iterations of a single consensus
    ch_genomes_mapped = ch_genome.multiMap { meta, fasta ->
        metadata: [meta.id, meta.subMap('id', 'cluster_id', 'sample')]
        fasta: [meta.id, fasta]
    }

    ch_genome_grouped = ch_genomes_mapped.fasta
        .collectFile { id, genome ->
            ["${id}.fa", genome]
        }
        .map { file -> [file.simpleName, file] }
        .join(ch_genomes_mapped.metadata.unique())
        .map { _id, genome, meta -> [meta, genome] }

    // Contig summary statistics
    if (!params.skip_quast) {
        QUAST_QC(ch_genome, [[:], []], [[:], []])
        ch_quast = QUAST_QC.out.tsv
    }

    // Identify closest reference from the reference pool database using blast
    if (!params.skip_blast_qc) {
        BLASTN_QC(ch_genomes_all, ch_refpool_db, [], [], [])
        ch_blast = BLASTN_QC.out.txt
    }

    // use MMSEQS easy search to find best hits against annotation db
    if (!params.skip_consensus_annotation) {
        MMSEQS_ANNOTATE(ch_genomes_all, ch_annotation_db)
        ch_annotation = MMSEQS_ANNOTATE.out.tsv
    }

    // Annotate proteins with prokka
    if (!params.skip_prokka) {
        // Run
        ch_genomes_final = ch_genome.filter { meta, _genome ->
            def isConstraint = meta.containsKey('isConstraint') ? meta.isConstraint : false
            isConstraint || meta.iteration == params.iterative_refinement_cycles
        }

        PROKKA(ch_genomes_final, ch_prokka_db, [])
    }


    // use Checkv to estimate Completeness and Contamination
    if (!params.skip_checkv) {
        if (!params.checkv_db) {
            CHECKV_DOWNLOADDATABASE()
            ch_checkv_db = CHECKV_DOWNLOADDATABASE.out.checkv_db
        }

        // uses HMM and AA alignment to deterimine completeness
        CHECKV_ENDTOEND(ch_genome_grouped, ch_checkv_db)
        ch_checkv = CHECKV_ENDTOEND.out.quality_summary
    }

    // Align the different steps to each other to see how the sequences have changed
    if (!params.skip_alignment_qc) {

        // MAFFT doesn't like those that have only one sequence
        ch_genome_grouped_branch = ch_genome_grouped.branch { _meta, scaffolds ->
            pass: scaffolds.countFasta() > 1
            fail: scaffolds.countFasta() == 1
        }

        MAFFT_ITERATIONS(ch_genome_grouped_branch.pass, [[:], []], [[:], []], [[:], []], [[:], []], [[:], []], false)

        ch_contigs_mod = ch_aligned_raw_contigs.map { meta, genome -> [meta.id, meta, genome] }

        // Make a channel that contains the alignment of the iterations with
        // the original contigs from the assemblers
        ch_mafftQC_in = ch_genome_grouped_branch.fail
            .mix(MAFFT_ITERATIONS.out.fas)
            .map { meta, genome -> [meta.id, meta, genome] }
            .join( ch_contigs_mod, by: 0)
            .filter { _id, _meta_genome, scaffolds, _meta_contigs, contigs ->
                // Make sure we have at least 2 sequences
                scaffolds.countFasta() + contigs.countFasta() > 1
            }
            .multiMap { _id, meta_genome, scaffolds, _meta_contigs, contigs ->
                scaffolds: [meta_genome, scaffolds]
                contigs: [meta_genome, contigs]
            }

        MAFFT_QC(ch_mafftQC_in.scaffolds, ch_mafftQC_in.contigs, [[:], []], [[:], []], [[:], []], [[:], []], false)
    }

    emit:
    blast      = ch_blast         // channel: [ val(meta), [ txt ] ]
    checkv     = ch_checkv        // channel: [ val(meta), [ tsv ] ]
    quast      = ch_quast         // channel: [ val(meta), [ tsv ] ]
    annotation = ch_annotation    // channel: [ val(meta), [ txt ] ]
}
