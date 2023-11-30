
include { CDHIT_CDHITEST           } from '../../modules/nf-core/cdhit/cdhitest/main'
include { VSEARCH_CLUSTER          } from '../../modules/nf-core/vsearch/cluster/main'
include { MMSEQS_CREATEDB          } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_LINCLUST          } from '../../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_CREATETSV         } from '../../modules/nf-core/mmseqs/createtsv/main'
include { VRHYME_VRHYME            } from '../../modules/nf-core/vrhyme/vrhyme/main'
include { GUNZIP as GUNZIP_CONTIGS } from '../../modules/nf-core/gunzip/main'

workflow FASTA_FASTQ_CLUST {

    take:
    fasta_fastq // channel: [ val(meta), [ fasta ], [ fastq ] ]
    cluster_method // string

    main:
    ch_versions = Channel.empty()
    ch_fasta    = fasta_fastq.map{ meta, fasta, fastq -> [meta, fasta] }   

    // cluster our reference hits and contigs should make this a subworkflow
    if ( cluster_method == "vsearch" ){
        VSEARCH_CLUSTER (ch_fasta)

        ch_clusters = VSEARCH_CLUSTER.out.uc
        ch_versions = ch_versions.mix(VSEARCH_CLUSTER.out.versions.first())
    }
    else if (cluster_method == "cdhitest") {
        CDHIT_CDHITEST(ch_fasta)

        ch_clusters = CDHIT_CDHITEST.out.clusters
        ch_versions = ch_versions.mix(CDHIT_CDHITEST.out.versions.first())
    }
    else if (cluster_method == "mmseqs") {
        MMSEQS_CREATEDB ( ch_fasta )
        ch_versions = ch_versions.mix(MMSEQS_CREATEDB.out.versions.first())
        MMSEQS_LINCLUST ( MMSEQS_CREATEDB.out.db )
        ch_versions = ch_versions.mix(MMSEQS_LINCLUST.out.versions.first())

        mmseqs = MMSEQS_LINCLUST
            .out
            .db_cluster
            .join(MMSEQS_CREATEDB.out.db, by: [0])
        
        mmseqs_cluster = mmseqs.map{ meta, db_cluster, db -> [meta, db_cluster] }
        mmseqs_db      = mmseqs.map{ meta, db_cluster, db -> [meta, db] }

        MMSEQS_CREATETSV (
            mmseqs_cluster,
            mmseqs_db,
            [[:],[]]
        )
        ch_versions = ch_versions.mix(MMSEQS_CREATETSV.out.versions.first())

        ch_clusters = MMSEQS_CREATETSV.out.tsv
    }
    else if (cluster_method == "vrhyme") { 

        GUNZIP_CONTIGS (ch_fasta)
        ch_versions = ch_versions.mix(GUNZIP_CONTIGS.out.versions.first())

        fasta_fastq
            .join(GUNZIP_CONTIGS.out.gunzip, by: [0])
            .map{ meta, fasta, fastq, gunzip -> [meta, gunzip, fastq] }
            .set{ch_gunzip_fastq}
        VRHYME_VRHYME (ch_gunzip_fastq)
        ch_clusters = VRHYME_VRHYME.out.membership
        ch_versions = ch_versions.mix(VRHYME_VRHYME.out.versions.first())
    }
    else {
        error "Unknown cluster method: ${cluster_method}"
    }
    emit:
    clusters =  ch_clusters // channel: [ [ meta ], [ clusters ] ]

    versions = ch_versions  // channel: [ versions.yml ]
}

