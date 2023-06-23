// modules
include { CLUSTER_EXTRACT                           } from '../../modules/local/cluster_extract'
include { SEQKIT_GREP as SEQKIT_GREP_MEMBERS        } from '../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_CENTROIDS     } from '../../modules/nf-core/seqkit/grep/main'

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

    // join with sequence data based on meta - transpose to a long format so we can split up the channels
    CLUSTER_EXTRACT.out.members_centroids
        .join(db_seq)
        .transpose() //wide to long
        .set { members_centroids_transposed }

    // Update the meta data to include the cluster ID
    ch_members_centroids = members_centroids_transposed.map { create_member_ref_channel(it) }

    // split up members, centroids and db for downstream processing
    ch_members_centroids
        .map { [ it[1] ] }
        .set { ch_members }
    ch_members_centroids
        .map { [ it[2] ] }
        .set { ch_centroids }
    ch_members_centroids
        .map { [ it[0], it[3] ] }
        .set { ch_db_seq_anno }

    SEQKIT_GREP_MEMBERS(ch_db_seq_anno,ch_members )
    ch_seq_members   =  SEQKIT_GREP_MEMBERS.out.filter

    SEQKIT_GREP_CENTROIDS(ch_db_seq_anno,ch_centroids )
    ch_seq_centroids = SEQKIT_GREP_CENTROIDS.out.filter

    ch_versions =  ch_versions.mix(SEQKIT_GREP_CENTROIDS.out.versions.first())

    emit:
    // seq_members     = ch_seq_members        // channel: [ [ meta ], [ seq_members.fa] ]
    // seq_centroids  = ch_seq_centroids     // channel: [ [ meta ], [ seq_centroids.fa ] ]
    versions        = ch_versions
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_member_ref_channel(ArrayList row) {
    new_meta = row[0].clone()
    members = row[1]
    centroids = row[2]
    sequence = row[3]

    new_meta.sample     = String.valueOf(new_meta.id) // just to make sure we don't pass by reference
    regex       = (members =~ /.*\/.*_([0-9]+)_n([0-9]+)_members.txt/) // try and extract the correct cluster ID and size associated to the sample
    new_meta.cluster    = regex[0][1]
    new_meta.size       = regex[0][2]
    new_meta.id = "${new_meta.sample}_${new_meta.cluster}"

    def result = [ new_meta, members, centroids, sequence]
    return result
}

