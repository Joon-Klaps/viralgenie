include { BLAST_BLASTN      } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_FILTER      } from '../../modules/local/blast_filter'
include { GUNZIP            } from '../../modules/nf-core/gunzip/main'
include { SEQKIT_GREP       } from '../../modules/nf-core/seqkit/grep/main'
include { CDHIT_CDHITEST    } from '../../modules/nf-core/cdhit/cdhitest/main'
include { VSEARCH_CLUSTER   } from '../../modules/nf-core/vsearch/cluster/main'
include { CAT_CAT           } from '../../modules/nf-core/cat/cat/main'
include { CLUST_SEQ_EXTRACT } from '../../subworkflows/local/clust_seq_extract'

workflow FASTA_BLAST_CLUST {

    take:
    fasta          // channel: [ val(meta), [ fasta.gz ] ]
    blast_db       // channel: [ val(meta), path(db) ]
    blast_db_fasta // channel: [ val(meta), path(fasta) ]
    cluster_method // string

    main:
    ch_versions = Channel.empty()

    // Blast results, to a reference database, to find a complete genome that's already assembled
    BLAST_BLASTN (
        fasta,
        blast_db
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    BLAST_BLASTN
        .out
        .txt
        .branch { meta, txt ->
            no_hits : txt.countLines() == 0
            hits    : txt.countLines() > 0
        }
        .set { ch_blast_txt }

    //TODO throw an warning if no hits are found

    BLAST_FILTER (
        ch_blast_txt.hits
    )
    ch_versions = ch_versions.mix(BLAST_FILTER.out.versions.first())

    // give the references a meta again so it can be used in seqkit
    blast_db_fasta
        .combine(BLAST_FILTER.out.hits)
        .view()

    blast_db_fasta
        .combine(BLAST_FILTER.out.hits)
        .map{ meta_blast, db_fasta, meta_filter, filter -> [meta_filter, db_fasta]
        }
        .set{ch_blast_db_fasta}

    ch_hits = BLAST_FILTER.out.hits.map{it -> it[1]}

    // isolate the hits from the database to a fasta file
    SEQKIT_GREP (ch_blast_db_fasta, ch_hits)
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    GUNZIP(SEQKIT_GREP.out.filter)
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())

    // put reference hits and contigs together
    fasta
        .mix(GUNZIP.out.gunzip)
        .groupTuple(size :2, remainder: false) // group and throw away those that don't fit
        .set {ch_reference_contigs_comb}

    CAT_CAT(ch_reference_contigs_comb)

    // cluster our reference hits and contigs
    if ( cluster_method == "vsearch" ){
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
    ch_centroids_members = CLUST_SEQ_EXTRACT.out.seq_centroids_members

    emit:
    clusters              = ch_clusters           // channel: [ [ meta ], [ clusters ] ]
    centroids_members     = ch_centroids_members  // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    versions              = ch_versions           // channel: [ versions.yml ]

}

