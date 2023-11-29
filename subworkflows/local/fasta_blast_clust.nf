include { BLAST_BLASTN      } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_FILTER      } from '../../modules/local/blast_filter'
include { GUNZIP            } from '../../modules/nf-core/gunzip/main'
include { SEQKIT_GREP       } from '../../modules/nf-core/seqkit/grep/main'
include { CDHIT_CDHITEST    } from '../../modules/nf-core/cdhit/cdhitest/main'
include { VSEARCH_CLUSTER   } from '../../modules/nf-core/vsearch/cluster/main'
include { CAT_CAT           } from '../../modules/nf-core/cat/cat/main'
include { MMSEQS_CREATEDB   } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_LINCLUST   } from '../../modules/nf-core/mmseqs/linclust/main'
include { MMSEQS_CREATETSV  } from '../../modules/nf-core/mmseqs/createtsv/main'
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

    // Throw an warning if no hits are found
    ch_blast_txt
        .hits
        .collect()
        .ifEmpty{ log.warn "No blast hits were found in any samples of the given BLAST database. Consider updating the search parameters or the database: \n ${params.reference_pool} "}

    // if (ch_blast_txt.hits.count() ==0 ){
    //     println(" IWAS HERRREREREZNJKGNDSJKGNJKDS VKSD VKLDS FKLDS L")
    //     Nextflow.warn("No blast hits were found in any samples of the given BLAST database. Consider updating the search parameters or the database: \n ${params.reference_pool} ")
    // }

    ch_blast_txt
        .no_hits
        .join(fasta)
        .map { meta, txt, fasta ->
            def n_fasta = fasta.countFasta()
            ["$meta.sample\t$n_fasta"]}
        .collect()
        .map {
            tsv_data ->
                def comments = [
                    "id: 'samples_without_blast_hits'",
                    "anchor: 'WARNING: Filtered samples'",
                    "section_name: 'Samples without blast hits'",
                    "format: 'tsv'",
                    "description: 'Samples that did not have any blast hits for their contigs (using ${params.assemblers}) were not included in further assembly polishing'",
                    "plot_type: 'table'"
                ]
                def header = ['Sample', "Number of contigs"]
                return WorkflowCommons.multiqcTsvFromList(tsv_data, header, comments) // make it compatible with the other mqc files
        }
        .collectFile(name:'samples_no_blast_hits_mqc.tsv')
        .set { no_blast_hits }

    BLAST_FILTER (
        ch_blast_txt.hits
    )
    ch_versions = ch_versions.mix(BLAST_FILTER.out.versions.first())

    // give the references a meta again so it can be used in seqkit
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

    // cluster our reference hits and contigs should make this a subworkflow
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
    else if (cluster_method == "mmseqs") {
        MMSEQS_CREATEDB ( CAT_CAT.out.file_out )
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

    CLUST_SEQ_EXTRACT(
        ch_clusters,
        cluster_method,
        CAT_CAT.out.file_out
    )
    ch_centroids_members = CLUST_SEQ_EXTRACT.out.seq_centroids_members

    emit:
    clusters              = ch_clusters                              // channel: [ [ meta ], [ clusters ] ]
    centroids_members     = ch_centroids_members                     // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    clusters_tsv          = CLUST_SEQ_EXTRACT.out.clusters_tsv       // channel: [ [ meta ], [ tsv ] ]
    clusters_summary      = CLUST_SEQ_EXTRACT.out.clusters_summary   // channel: [ [ meta ], [ tsv ] ]
    no_blast_hits_mqc     = no_blast_hits                            // channel: [ tsv ]
    versions              = ch_versions                              // channel: [ versions.yml ]

}

