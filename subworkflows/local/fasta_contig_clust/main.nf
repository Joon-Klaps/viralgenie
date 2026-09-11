include { FASTA_BLAST_REFSEL      } from '../fasta_blast_refsel'
include { FASTA_FASTQ_CLUST       } from '../fasta_fastq_clust'
include { FASTA_CONTIG_PRECLUST   } from '../fasta_contig_preclust'
include { EXTRACT_CLUSTER         } from '../../../modules/local/extract_cluster/main'
include { getMapFromJson          } from '../utils_nfcore_viralmetagenome_pipeline'

workflow FASTA_CONTIG_CLUST {

    take:
    ch_fasta_fastq              // channel: [ val(meta), [ fasta ],  [ fastq ] ]
    ch_coverages                // channel: [ val(meta), [ idxstats* ] ]
    ch_blacklist                // channel: [ path(blacklist) ]
    ch_blast_db                 // channel: [ val(meta), path(db) ]
    ch_blast_db_fasta           // channel: [ val(meta), path(fasta) ]
    ch_kraken2_db               // channel: [ val(meta), path(db) ]
    ch_kaiju_db                 // channel: [ val(meta), path(db) ]
    contig_classifiers          // value ['kraken2','kaiju']
    cluster_method              // value ['vsearch','cdhitest','mmseqs-linclust','mmseqs-cluster','vrhyme', 'mash']
    identity_threshold          // value: 0.85
    skip_precluster             // boolean
    perc_reads_contig           // value: 5
    cluster_with_reference_pool // boolean: whether blast db is provided (if not, skip blast ref selection)

    main:
    ch_no_blast_hits     = channel.empty()
    ch_fasta             = ch_fasta_fastq.map{ meta, fasta, _fastq -> [meta, fasta] }
    ch_fasta_ref_contigs = ch_fasta

    if ( cluster_with_reference_pool ) {
        // Blast contigs to a reference database, to find a reference genome can be used for scaffolding
        FASTA_BLAST_REFSEL (
            ch_fasta,
            ch_blacklist,
            ch_blast_db,
            ch_blast_db_fasta
        )
        ch_no_blast_hits     = FASTA_BLAST_REFSEL.out.no_blast_hits
        ch_fasta_ref_contigs = FASTA_BLAST_REFSEL.out.fasta_ref_contigs
    }

    // Combine with reads if vrhyme is used
    ch_contigs_reads = ch_fasta_ref_contigs
        .join(ch_fasta_fastq, by: [0])
        .map{meta, contigs_joined, _contigs, reads -> [meta + [ntaxa: 1], contigs_joined, reads]} // ntaxa will use later

    // precluster our reference hits and contigs using kraken & Kaiju to delineate contigs at a species level
    if (!skip_precluster) {
        FASTA_CONTIG_PRECLUST (
            ch_contigs_reads,
            contig_classifiers,
            ch_kaiju_db,
            ch_kraken2_db
        )
        ch_contigs_reads = FASTA_CONTIG_PRECLUST.out.contigs_reads
    }

    // cluster our reference hits and contigs should make this a subworkflow
    FASTA_FASTQ_CLUST (
        ch_contigs_reads,
        cluster_method,
        identity_threshold
    )

    // if we have no coverage files, make the empty array else join with coverages
    if (perc_reads_contig == 0){
        sample_fasta_ref_contigs = ch_fasta_ref_contigs
            .map{ meta, fasta -> [meta.sample, meta, fasta,[]] }               // add sample for join
    } else {
        sample_coverages = ch_coverages
            .map{ meta, idxstats -> [meta.sample, meta, idxstats] }            // add sample for join

        sample_fasta_ref_contigs = ch_fasta_ref_contigs
            .map{ meta, fasta -> [meta.sample, meta, fasta] }                  // add sample for join
            .join(sample_coverages, by: [0])                                   // join with coverages
            .map{ sample, meta_fasta, fasta, _meta_coverages, coverages ->     // remove meta coverages
                [sample, meta_fasta, fasta, coverages]
                }
    }

    // Join cluster files with contigs & group based on number of preclusters (ntaxa)
    ch_clusters_contigs_coverages = FASTA_FASTQ_CLUST
        .out
        .clusters
        .map{ meta, clusters ->
            tuple( groupKey(meta.sample, meta.ntaxa), meta, clusters )         // Set groupkey by sample and ntaxa
            }
        .groupTuple(remainder: true)                                           // Has to be grouped to link different taxa preclusters to the same sample
        .combine(sample_fasta_ref_contigs)                                     // combine with contigs (regural join doesn't work)
        .filter{ sample, _meta_clust, _clusters, sample2, _meta_contig, _contigs, _coverages ->
            sample==sample2 }                                                  // filter for matching samples
        .map{ _sample, _meta_clust, clusters, _sample2, meta_contig, contigs, coverages ->
            [meta_contig, clusters, contigs, coverages]                        // get rid of meta_clust & sample
        }

    EXTRACT_CLUSTER (
        ch_clusters_contigs_coverages,
        cluster_method
    )

    ch_seq_centroids_members = EXTRACT_CLUSTER
        .out
        .members_centroids
        .transpose()                                                                   // wide to long
        .map { meta, seq_members, seq_centroids, json_file ->
            def json = getMapFromJson(json_file)                                       // deep-cloned HashMap; safe to drop into meta
            def map_json = [
                id                : "${meta.sample}_${json.cluster_id}",               // rename meta.id to include cluster number
                centroid          : json.centroid?.toString(),
                cluster_id        : json.cluster_id?.toString(),
                cluster_size      : (json.cluster_size as Integer),
                external_reference: (json.external_reference as Boolean),
                taxid             : json.taxid?.toString()
            ]
            return [meta + map_json, seq_centroids, seq_members]
        }

    emit:
    clusters              = FASTA_FASTQ_CLUST.out.clusters // channel: [ [ meta ], [ clusters ] ]
    centroids_members     = ch_seq_centroids_members       // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    clusters_tsv          = EXTRACT_CLUSTER.out.tsv        // channel: [ [ meta ], [ tsv ] ]
    clusters_summary      = EXTRACT_CLUSTER.out.summary    // channel: [ [ meta ], [ tsv ] ]
    no_blast_hits_mqc     = ch_no_blast_hits               // channel: [ tsv ]

}
