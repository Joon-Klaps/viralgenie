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

    SEQKIT_GREP_CENTROIDS.out.filter
        .join(SEQKIT_GREP_MEMBERS.out.filter, remainder: true)
        .transpose() //wide to long
        .join(CLUSTER_EXTRACT.out.json)
        .map { create_member_ref_channel(it) }
        .set { seq_centroids_members }

    seq_centroids_members.view()

    emit:
    seq_centroids_members     = seq_centroids_members        // channel: [ [ meta ], [ seq_centroids.fa], [ seq_members.fa] ]
    versions                  = ch_versions
}

//
// Function that parses fastp json output file to get total number of reads after trimming
//
import groovy.json.JsonSlurper

def getClustersFromJson(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text)
    return json
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_member_ref_channel(ArrayList row) {
    meta         = row[0]
    centroids    = row[1]
    members      = row[2]
    json         = getClustersFromJson(row[3])

    sample       = String.valueOf(meta.id) // just to make sure we don't pass by reference
    cluster      = json['id']
    id           = "${sample}_${cluster}"

    new_meta = meta + [ id: id, cluster: cluster, cluster_size: json['size'], sample: sample, centroid: json['centroid'], members: json['members'] ]

    def result = [ new_meta, centroids, members]
    return result
}

