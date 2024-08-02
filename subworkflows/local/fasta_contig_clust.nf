include { FASTA_BLAST_REFSEL      } from '../../subworkflows/local/fasta_blast_refsel'
include { FASTA_FASTQ_CLUST       } from '../../subworkflows/local/fasta_fastq_clust'
include { FASTA_CONTIG_PRECLUST   } from '../../subworkflows/local/fasta_contig_preclust'
include { EXTRACT_CLUSTER         } from '../../modules/local/extract_cluster/main'

workflow FASTA_CONTIG_CLUST {

    take:
    fasta_fastq        // channel: [ val(meta), [ fasta ],  [ fastq ] ]
    blast_db           // channel: [ val(meta), [ db ], [ fasta ] ]
    coverages          // channel: [ val(meta), [ idxstats* ] ]
    contig_classifiers // value:   [ kaiju, kraken2 ]
    kraken2_db         // channel: [ val(meta), path(db) ]
    kaiju_db           // channel: [ val(meta), path(db) ]

    main:
    ch_versions = Channel.empty()
    fasta       = fasta_fastq.map{ meta, fasta, fastq -> [meta, fasta] }

    blast_db
    .map{ meta, db, fasta ->
        samples = meta.samples == null ? meta.samples: tuple(meta.samples.split(";"))  // Split up samples if meta.samples is not null
        [meta, samples, db, fasta]
    }
    .transpose(remainder:true)                                                         // unset
    .set{blast_db}

    // combine refpool_db based on specified samples in the reference_pools parameter
    fasta_fastq
        .combine(blast_db)
        .branch {
            meta_genome, genome, reads, meta_db, db_samples, blast_db, blast_seq ->
            specific: meta_genome.sample == db_samples
            generic: db_samples == null
        }.set{db}

    //combine specific and generic branches
    db.generic
        .mix(db.specific)
        .multiMap{ meta_genome, genome, reads, meta_db, db_samples, blast_db, blast_seq ->
            new_id = meta_genome.id + '-' +  meta_db.id
            contig: [meta_genome + [id: new_id, db_comb: new_id ], genome, reads]
            db: [meta_genome + [id: new_id, db_comb: new_id], blast_db]
            db_fasta: [meta_genome + [id: new_id, db_comb: new_id], blast_seq]
        }
        .set{ch_blastrefsel_in}

    // Blast contigs to a reference database, to find a reference genome can be used for scaffolding
    FASTA_BLAST_REFSEL (
        ch_blastrefsel_in.contig,
        ch_blastrefsel_in.db,
        ch_blastrefsel_in.db_fasta
    )
    ch_versions       = ch_versions.mix(FASTA_BLAST_REFSEL.out.versions)
    no_blast_hits     = FASTA_BLAST_REFSEL.out.no_blast_hits
    fasta_sel_fastq   = FASTA_BLAST_REFSEL.out.fasta_sel_fastq

    fasta_sel_fastq
        .map{meta, contigs_joined, contigs, reads -> [meta + [ntaxa: 1], contigs_joined, reads]} // ntaxa will use later
        .set{ch_contigs_reads}

    // precluster our reference hits and contigs using kraken & Kaiju to delineate contigs at a species level
    if (!params.skip_precluster) {
        FASTA_CONTIG_PRECLUST (
            ch_contigs_reads,
            contig_classifiers,
            kaiju_db,
            kraken2_db
        )
        ch_versions      = ch_versions.mix(FASTA_CONTIG_PRECLUST.out.versions)
        ch_contigs_reads = FASTA_CONTIG_PRECLUST.out.contigs_reads
    }

    // cluster our reference hits and contigs should make this a subworkflow
    FASTA_FASTQ_CLUST (
        ch_contigs_reads,
        params.cluster_method,
    )
    ch_versions = ch_versions.mix(FASTA_FASTQ_CLUST.out.versions)

    // if we have no coverage files, make the empty array else join with coverages
    if (params.perc_reads_contig == 0){
        sample_fasta_ref_contigs = ch_contigs_reads
            .map{ meta, fasta, reads -> [meta.db_comb, meta, fasta []] }       // add sample for join
    } else {
        sample_coverages = coverages
            .map{ meta, idxstats -> [meta.sample, meta, idxstats] }            // add sample for join

        sample_fasta_ref_contigs = ch_contigs_reads
            .map{ meta, fasta, reads -> [meta.sample, meta, fasta] }           // add sample for join
            .combine(sample_coverages)                                         // combine with coverages (need an carhesian product)
            .filter{it -> it[0]==it[3]}                                        // filter for matching samples
            .map{ sample, meta_fasta, fasta, sample_coverages, meta_coverages, coverages ->
                [meta_fasta.db_comb, meta_fasta, fasta, coverages]             // remove meta coverages
                }
    }

    sample_fasta_ref_contigs.dump(tag:"sample_fasta_ref_contigs", pretty:true)

    // Join cluster files with contigs & group based on number of preclusters (ntaxa)
    FASTA_FASTQ_CLUST
        .out
        .clusters
        .map{ meta, clusters ->
            tuple( groupKey(meta.db_comb, meta.ntaxa), meta, clusters )        // Set groupkey by sample-db and ntaxa
            }
        .groupTuple(remainder: true)                                           // Has to be grouped to link different taxa preclusters to the same sample
        .combine(sample_fasta_ref_contigs)                                     // combine with contigs (regural join doesn't work)
        .set{tmp}
        tmp.dump(tag:"tmp", pretty:true)
        tmp
        .filter{it -> it[0]==it[3]}                                            // filter for matching samples
        .map{ db_comb, meta_clust, clusters, sample2, meta_contig, contigs, coverages ->
            [meta_contig, clusters, contigs, coverages]                        // get rid of meta_clust & sample
        }
        .set{ch_clusters_contigs_coverages}

    EXTRACT_CLUSTER (
        ch_clusters_contigs_coverages,
        params.cluster_method
    )
    ch_versions = ch_versions.mix(EXTRACT_CLUSTER.out.versions.first())

    EXTRACT_CLUSTER
        .out
        .members_centroids
        .transpose()                                                                     // wide to long
        .map { meta, seq_members, seq_centroids, json_file ->
            json          = WorkflowCommons.getMapFromJson(json_file)                    // convert cluster metadata to Map
            new_meta      = meta + [ id: "${meta.db_comb}_${json.cluster_id}"] + json    // rename meta.id to include cluster number
            return [new_meta, seq_centroids, seq_members]
        }
        .set{seq_centroids_members}

    emit:
    clusters              = FASTA_FASTQ_CLUST.out.clusters            // channel: [ [ meta ], [ clusters ] ]
    centroids_members     = seq_centroids_members                     // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    clusters_tsv          = EXTRACT_CLUSTER.out.tsv                   // channel: [ [ meta ], [ tsv ] ]
    clusters_summary      = EXTRACT_CLUSTER.out.summary               // channel: [ [ meta ], [ tsv ] ]
    no_blast_hits_mqc     = no_blast_hits                            // channel: [ tsv ]
    versions              = ch_versions                              // channel: [ versions.yml ]

}

