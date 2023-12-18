// modules
include { EXTRACT_CLUSTER                             } from '../../modules/local/extract_cluster/main'
include { SEQKIT_GREP_MULTI as SEQKIT_GREP_MEMBERS    } from '../../modules/local/seqkit_grep_multi'
include { SEQKIT_GREP_MULTI as SEQKIT_GREP_CENTROIDS  } from '../../modules/local/seqkit_grep_multi'

workflow CLUST_SEQ_EXTRACT {

    take:
    ch_clusters     // channel: [ [ meta ], [cluster1, cluster2, ... ], [ fasta ] ]
    cluster_method  // value: cdhitest| vsearch | mmseqs | ...

    main:
    ch_versions         = Channel.empty()

    // extract members and centroids from clusters
    EXTRACT_CLUSTER(ch_clusters, cluster_method)
    ch_versions =  ch_versions.mix(EXTRACT_CLUSTER.out.versions.first())

    db_seq_reads
        .map{ meta, seq, reads -> [meta, seq] }
        .set{db_seq}

    EXTRACT_CLUSTER
        .out
        .members_centroids
        .join(db_seq, by: [0])
        .set{members_centroids_seq}
    // [sample_meta, [members], [centroids] ]

    ch_centroids  = members_centroids_seq.map {meta, members, centroids,seq -> [meta, centroids,seq] }
    ch_members    = members_centroids_seq.map {meta, members, centroids,seq -> [meta, members,seq] }

    // extract members and centroids from db
    SEQKIT_GREP_MEMBERS (ch_members )
    ch_versions =  ch_versions.mix(SEQKIT_GREP_MEMBERS.out.versions.first())

    SEQKIT_GREP_CENTROIDS( ch_centroids)
    ch_versions =  ch_versions.mix(SEQKIT_GREP_CENTROIDS.out.versions.first())

    //
    // NEXT STEPS are to make [ meta, [ centroids, members ] ]
    // Here the meta contains metadata on the sample and cluster, so instead of $it represinting samples, $it represents a single cluster
    // We will use the json file from EXTRACT_CLUSTER to get the cluster id and the necessary metadata
    // > [ cluster_id, centroid, members, cluster_size]
    SEQKIT_GREP_CENTROIDS.out.filter
        .join(SEQKIT_GREP_MEMBERS.out.filter, remainder: true)
        .join(EXTRACT_CLUSTER.out.json, remainder: true)
        .transpose() //wide to long
        .map { meta, seq_centroids, seq_members, json_file ->
            id            = seq_centroids.baseName.replace('_centroid','')
            json          = WorkflowCommons.getMapFromJson(json_file)
            new_meta      = meta + [ id: String.valueOf(id)] + json
            return [new_meta, seq_centroids, seq_members]
        }.set{seq_centroids_members}

    emit:
    seq_centroids_members     = seq_centroids_members       // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    clusters_tsv              = EXTRACT_CLUSTER.out.tsv     // channel: [ [ meta ], [ tsv ] ]
    clusters_summary          = EXTRACT_CLUSTER.out.summary // channel: [ [ meta ], [ tsv ] ]
    versions                  = ch_versions
}
