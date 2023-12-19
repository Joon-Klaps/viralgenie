include { FASTA_BLAST_EXTRACT     } from '../../subworkflows/local/fasta_blast_extract'
include { FASTA_FASTQ_CLUST       } from '../../subworkflows/local/fasta_fastq_clust'
include { FASTA_CONTIG_PRECLUST   } from '../../subworkflows/local/fasta_contig_preclust'
include { CLUST_SEQ_EXTRACT       } from '../../subworkflows/local/clust_seq_extract'
include { EXTRACT_CLUSTER         } from '../../modules/local/extract_cluster/main'

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
    ch_versions       = ch_versions.mix(FASTA_BLAST_EXTRACT.out.versions)
    no_blast_hits     = FASTA_BLAST_EXTRACT.out.no_blast_hits
    fasta_ref_contigs = FASTA_BLAST_EXTRACT.out.fasta_ref_contigs

    // Combine with reads if vrhyme is used
    fasta_ref_contigs
        .join(fasta_fastq, by: [0])
        .map{meta, contigs_joined, contigs, reads -> [meta + [ntaxa: 1], contigs_joined, reads]} // ntaxa will use later
        .set{ch_contigs_reads}

    // precluster our reference hits and contigs using kraken & Kaiju to delineate contigs at a species level
    if (!params.skip_precluster) {
        FASTA_CONTIG_PRECLUST (
            ch_contigs_reads,
            kaiju_db,
            kraken2_db
        )
        ch_versions      = ch_versions.mix(FASTA_CONTIG_PRECLUST.out.versions)
        ch_contigs_reads = FASTA_CONTIG_PRECLUST.out.sequences_reads
    }

    // cluster our reference hits and contigs should make this a subworkflow
    FASTA_FASTQ_CLUST (
        ch_contigs_reads,
        params.cluster_method,
    )
    ch_versions = ch_versions.mix(FASTA_FASTQ_CLUST.out.versions)

    // Join cluster files with contigs & group based on number of preclusters (ntaxa)
    fasta_ref_contigs
        .map{ meta, fasta -> [meta.sample, meta, fasta] }                       // add sample for join
        .set{sample_fasta_ref_contigs}

    FASTA_FASTQ_CLUST
        .out
        .clusters
        .map{ meta, clusters ->
            tuple( groupKey(meta.sample, meta.ntaxa), meta, clusters )         // Set groupkey by sample and ntaxa
            }
        .groupTuple(remainder: true)
        .join(sample_fasta_ref_contigs, by: [0])                               // join with contigs
        .map{ sample, meta_clust, clusters, meta_contig, contigs ->
             [meta_contig, clusters, contigs]                                  // get rid of meta_clust & sample
        }
        .set{ch_clusters_contigs}

    ch_clusters_contigs.view{it -> "ch_clusters_contigs: $it"}

    EXTRACT_CLUSTER (
        ch_clusters_contigs,
        params.cluster_method
    )
    // CLUST_SEQ_EXTRACT(
    //     ch_clusters_contigs,
    //     params.cluster_method
    // )

    // ch_versions          = ch_versions.mix(CLUST_SEQ_EXTRACT.out.versions)
    // ch_centroids_members = CLUST_SEQ_EXTRACT.out.seq_centroids_members
    ch_centroids_members = ch_clusters_contigs
    clusters_summary = channel.empty()

    emit:
    clusters              = FASTA_FASTQ_CLUST.out.clusters           // channel: [ [ meta ], [ clusters ] ]
    centroids_members     = ch_centroids_members                     // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    // clusters_tsv          = CLUST_SEQ_EXTRACT.out.clusters_tsv       // channel: [ [ meta ], [ tsv ] ]
    // clusters_summary      = CLUST_SEQ_EXTRACT.out.clusters_summary   // channel: [ [ meta ], [ tsv ] ]
    clusters_summary      = clusters_summary   // channel: [ [ meta ], [ tsv ] ]
    no_blast_hits_mqc     = no_blast_hits                            // channel: [ tsv ]
    versions              = ch_versions                              // channel: [ versions.yml ]

}

