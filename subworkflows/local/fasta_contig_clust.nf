include { FASTA_BLAST_EXTRACT     } from '../../subworkflows/local/fasta_blast_extract'
include { FASTA_FASTQ_CLUST       } from '../../subworkflows/local/fasta_fastq_clust'
include { FASTA_CONTIG_PRECLUST   } from '../../subworkflows/local/fasta_contig_preclust'
include { CLUST_SEQ_EXTRACT       } from '../../subworkflows/local/clust_seq_extract'

workflow FASTA_CONTIG_CLUST {

    take:
    fasta_fastq     // channel: [ val(meta), [ fasta ],  [ fastq ] ]
    blast_db        // channel: [ val(meta), path(db) ]
    blast_db_fasta  // channel: [ val(meta), path(fasta) ]
    kraken2_db      // channel: [ val(meta), path(db) ]
    kaiju_db        // channel: [ val(meta), path(db) ]

    main:
    ch_versions = Channel.empty()
    fasta       = fasta_fastq.map{ meta, fasta, fastq -> [meta, fasta] }

    // Blast contigs to a reference database, to find a reference genome can be used for scaffolding
    FASTA_BLAST_EXTRACT (
        fasta,
        blast_db,
        blast_db_fasta
    )
    ch_versions   = ch_versions.mix(FASTA_BLAST_EXTRACT.out.versions)
    no_blast_hits = FASTA_BLAST_EXTRACT.out.no_blast_hits

    // precluster our reference hits and contigs using kraken & Kaiju to delineate contigs at a species level
    if (!params.skip_precluster) {
        FASTA_CONTIG_PRECLUST (
            fasta,
            kaiju_db,
            kraken2_db
        )
        ch_versions = ch_versions.mix(FASTA_CONTIG_PRECLUST.out.versions)
        fasta_fastq = FASTA_CONTIG_PRECLUST.out.classifications
    }
    fasta_ref_contigs = FASTA_BLAST_EXTRACT.out.fasta_ref_contigs

    // Combine with reads if vrhyme is used
    fasta_ref_contigs
        .join(fasta_fastq, by: [0])
        .map{meta, contigs_joined, contigs, reads -> [meta, contigs_joined, reads]}
        .set{ch_contigs_reads}

    // cluster our reference hits and contigs should make this a subworkflow
    FASTA_FASTQ_CLUST (
        ch_contigs_reads,
        params.cluster_method,
    )
    ch_versions = ch_versions.mix(FASTA_FASTQ_CLUST.out.versions)

    CLUST_SEQ_EXTRACT(
        FASTA_FASTQ_CLUST.out.clusters,
        params.cluster_method,
        fasta_ref_contigs
    )
    ch_versions = ch_versions.mix(CLUST_SEQ_EXTRACT.out.versions)
    ch_centroids_members = CLUST_SEQ_EXTRACT.out.seq_centroids_members

    emit:
    clusters              = FASTA_FASTQ_CLUST.out.clusters           // channel: [ [ meta ], [ clusters ] ]
    centroids_members     = ch_centroids_members                     // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    clusters_tsv          = CLUST_SEQ_EXTRACT.out.clusters_tsv       // channel: [ [ meta ], [ tsv ] ]
    clusters_summary      = CLUST_SEQ_EXTRACT.out.clusters_summary   // channel: [ [ meta ], [ tsv ] ]
    no_blast_hits_mqc     = no_blast_hits                            // channel: [ tsv ]
    versions              = ch_versions                              // channel: [ versions.yml ]

}

