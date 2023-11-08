// modules
include { CLUSTER_EXTRACT                       } from '../../modules/local/cluster_extract'
include { SEQKIT_GREP_MULTI as SEQKIT_GREP_MEMBERS    } from '../../modules/local/seqkit_grep_multi'
include { SEQKIT_GREP_MULTI as SEQKIT_GREP_CENTROIDS  } from '../../modules/local/seqkit_grep_multi'

workflow CLUST_SEQ_EXTRACT {

    take:
    ch_clusters     // channel: [ [ meta ], [ clusters ] ]
    cluster_method  // value: cdhitest| vsearch
    db_seq          // channel: [ [ meta ], [ db ] ]


    main:
    ch_versions         = Channel.empty()

    // extract members and centroids from clusters
    CLUSTER_EXTRACT(ch_clusters, cluster_method)
    ch_versions =  ch_versions.mix(CLUSTER_EXTRACT.out.versions.first())

    CLUSTER_EXTRACT
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
    // We will use the json file from CLUSTER_EXTRACT to get the cluster id and the necessary metadata
    // > [ cluster_id, centroid, members, cluster_size]
    SEQKIT_GREP_CENTROIDS.out.filter
        .join(SEQKIT_GREP_MEMBERS.out.filter, remainder: true)
        .join(CLUSTER_EXTRACT.out.json, remainder: true)
        .transpose() //wide to long
        .map { meta, seq_centroids, seq_members, json_file ->
            id            = seq_centroids.baseName.replace('_centroid','')
            json          = getClustersFromJson(json_file)
            new_meta      = meta + [ id: String.valueOf(id)] + json
            return [new_meta, seq_centroids, seq_members]
        }.set{seq_centroids_members}

    emit:
    seq_centroids_members     = seq_centroids_members       // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    versions                  = ch_versions
}

//
// Function that parses fastp json output file to get total number of reads after trimming
//
import groovy.json.JsonSlurper

def getClustersFromJson(json_file) {
    def Map json = (Map) new JsonSlurper().parse(json_file)
    return json
}

//
// Function to get list of [ meta, [ centroids, members ] ]
// Here the meta contains metadata on the sample and cluster, so instead of $it represinting samples, $it represents a single cluster
//
def create_member_ref_channel(row) {
    meta         = row[0]
    centroids    = row[1]
    members      = row[2]

    sample       = String.valueOf(meta.id) // just to make sure we don't pass by reference
    id           = centroids.baseName.replaceAll("_centroids.fa", "")

    new_meta = meta + [ id: id]

    def result = [ new_meta, centroids, members]
    return result
}

