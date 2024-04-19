
include { CDHIT_CDHITEST           } from '../../modules/nf-core/cdhit/cdhitest/main'
include { VSEARCH_CLUSTER          } from '../../modules/nf-core/vsearch/cluster/main'
include { MMSEQS_CREATEDB          } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_LINCLUST          } from '../../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_CLUSTER           } from '../../modules/nf-core/mmseqs/cluster/main'
include { MMSEQS_CREATETSV         } from '../../modules/nf-core/mmseqs/createtsv/main'
include { VRHYME_VRHYME            } from '../../modules/nf-core/vrhyme/vrhyme/main'
include { MASH_DIST                } from '../../modules/nf-core/mash/dist/main'
include { NETWORK_CLUSTER          } from '../../modules/local/network_cluster/main'

workflow FASTA_FASTQ_CLUST {

    take:
    fasta_fastq // channel: [ val(meta), [ fasta ], [ fastq ] ]
    cluster_method // string

    main:
    ch_versions = Channel.empty()
    ch_fasta    = fasta_fastq.map{ meta, fasta, fastq -> [meta, fasta] }
    ch_dist     = Channel.empty()

    // cluster our reference hits and contigs should make this a subworkflow
    if ( cluster_method == "vsearch" ){
        VSEARCH_CLUSTER (ch_fasta)

        ch_clusters = VSEARCH_CLUSTER.out.uc
        ch_versions = ch_versions.mix(VSEARCH_CLUSTER.out.versions.first())
    }
    else if (cluster_method == "cdhitest") {
        if (params.identity_threshold < 0.80 ) {
            log.warn "cdhitest identity threshold is set to ${params.identity_threshold}, which is below the minimum threshold of 0.80.\n AUTOFIX: Defaulting to 0.80"
        }
        CDHIT_CDHITEST(ch_fasta)

        ch_clusters = CDHIT_CDHITEST.out.clusters
        ch_versions = ch_versions.mix(CDHIT_CDHITEST.out.versions.first())
    }
    else if (cluster_method == "mmseqs-linclust" || cluster_method == "mmseqs-cluster") {
        MMSEQS_CREATEDB ( ch_fasta )
        ch_versions = ch_versions.mix(MMSEQS_CREATEDB.out.versions.first())

        if (cluster_method == "mmseqs-linclust") {
            MMSEQS_LINCLUST ( MMSEQS_CREATEDB.out.db )
            ch_versions = ch_versions.mix(MMSEQS_LINCLUST.out.versions.first())
            MMSEQS_LINCLUST
                .out
                .db_cluster
                .join(MMSEQS_CREATEDB.out.db, by: [0])
                .set{ mmseqs_cluster } // channel: [ [ meta ], [ db_cluster ], [ db ] ]
        }
        else {
            MMSEQS_CLUSTER ( MMSEQS_CREATEDB.out.db )
            ch_versions = ch_versions.mix(MMSEQS_CLUSTER.out.versions.first())
            MMSEQS_CLUSTER
                .out
                .db_cluster
                .join(MMSEQS_CREATEDB.out.db, by: [0])
                .set{ mmseqs_cluster } // channel: [ [ meta ], [ db_cluster ], [ db ] ]
        }

        MMSEQS_CREATETSV ( mmseqs_cluster )
        ch_versions = ch_versions.mix(MMSEQS_CREATETSV.out.versions.first())

        ch_clusters = MMSEQS_CREATETSV.out.tsv
    }
    else if (cluster_method == "vrhyme") {
        VRHYME_VRHYME (fasta_fastq)
        ch_clusters = VRHYME_VRHYME.out.membership
        ch_versions = ch_versions.mix(VRHYME_VRHYME.out.versions.first())
    }
    else if (cluster_method == "mash") {

        // Calculate distances
        MASH_DIST(ch_fasta)
        ch_versions = ch_versions.mix(MASH_DIST.out.versions.first())
        ch_dist = MASH_DIST.out.dist
    }
    else {
        error "Unknown cluster method: ${cluster_method}"
    }

    // Calculate clusters for distance based methods (mash)
    if (cluster_method in ["mash"] ) {
        // Calculate clusters using leidenalg
        NETWORK_CLUSTER(ch_dist, cluster_method, params.network_clustering)
        ch_clusters = NETWORK_CLUSTER.out.clusters
        ch_versions = ch_versions.mix(NETWORK_CLUSTER.out.versions.first())
    }

    emit:
    clusters =  ch_clusters // channel: [ [ meta ], [ clusters ] ]

    versions = ch_versions  // channel: [ versions.yml ]
}

