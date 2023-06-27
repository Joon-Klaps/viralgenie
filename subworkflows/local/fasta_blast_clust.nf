include { GUNZIP            } from '../../modules/nf-core/gunzip/main'
include { BLAST_BLASTN      } from '../../modules/local/blast_blastn'
include { SEQKIT_GREP       } from '../../modules/nf-core/seqkit/grep/main'
include { CDHIT_CDHITEST      } from '../../modules/nf-core/cdhit/cdhitest/main'
include { VSEARCH_CLUSTER   } from '../../modules/nf-core/vsearch/cluster/main'
include { CAT_CAT           } from '../../modules/nf-core/cat/cat/main'
include { CLUST_SEQ_EXTRACT } from '../../subworkflows/local/clust_seq_extract'

workflow FASTA_BLAST_CLUST {

    take:
    fasta          // channel: [ val(meta), [ fasta.gz ] ]
    blast_db       // channel: [ path(db) ]
    blast_db_fasta // channel: [path(db) ]
    cluster_method // string

    main:
    ch_versions = Channel.empty()

    GUNZIP(fasta)
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())
    ch_fasta_gunzip = GUNZIP.out.gunzip

    // Blast results, to a reference database, to find a complete genome that's already assembled
    BLAST_BLASTN (
        ch_fasta_gunzip,
        blast_db
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    // give the references a meta again so it can be used in seqkit
    blast_db_fasta
        .combine(BLAST_BLASTN.out.txt)
        .map{
            it -> [it[1],it[0]]
        }
        .set{ch_blast_db_fasta}

    ch_hits = BLAST_BLASTN.out.txt.map{it -> it[1]}
    // isolate the hits from the database to a fasta file
    SEQKIT_GREP (ch_blast_db_fasta, ch_hits)
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    // put reference hits and contigs together
    fasta.join(
            SEQKIT_GREP.out.filter,
            by: [0],
            remainder: true
            )
    .map {
            it ->
                [it[0],it[1..-1]]
            }
        .set {ch_reference_contigs_comb}
    CAT_CAT(ch_reference_contigs_comb)

    // cluster our reference hits and contigs
    if (cluster_method == "vsearch"){
        VSEARCH_CLUSTER (CAT_CAT.out.file_out)
        ch_clusters = VSEARCH_CLUSTER.out.uc
        ch_versions = ch_versions.mix(VSEARCH_CLUSTER.out.versions.first())
    }
    else if (cluster_method == "cdhitest") {
        CDHIT_CDHITEST(CAT_CAT.out.file_out)
        ch_clusters = CDHIT_CDHITEST.out.clusters
        ch_versions = ch_versions.mix(CDHIT_CDHITEST.out.versions.first())
    }

    CLUST_SEQ_EXTRACT(
        ch_clusters,
        cluster_method,
        CAT_CAT.out.file_out
    )

    emit:
    clusters    = ch_clusters                           // channel: [ [ meta ], [ clusters ] ]
    members     = CLUST_SEQ_EXTRACT.out.seq_members     // channel: [ [ meta ], [ seq_members.fa] ]
    centroids   = CLUST_SEQ_EXTRACT.out.seq_centroids   // channel: [ [ meta ], [ seq_centroids.fa ] ]
    versions    = ch_versions                           // channel: [ versions.yml ]

}
