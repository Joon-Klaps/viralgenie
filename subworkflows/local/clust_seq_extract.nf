// modules
include { CLUSTER_EXTRACT                           } from '../../modules/nf-core/local/cluster_extract'
include { SEQKIT_GREP as SEQKIT_GREP_MEMBERS        } from '../../modules/nf-core/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_REFERENCES     } from '../../modules/nf-core/nf-core/seqkit/grep/main'

workflow CLUST_SEQ_EXTRACT {

    take:
    ch_clusters     // channel: [ [ meta ], [ clusters ] ]
    cluster_method  // value: cdhit | vsearch
    ch_db           // channel: [ [ meta ], [ db ] ]


    main:
    ch_versions         = Channel.empty()

    // extract members and references from clusters
    CLUSTER_EXTRACT(ch_clusters, cluster_method)
        .out.members_references
        .map { create_member_ref_channel(it) }
        .set { ch_members_references }

    ch_versions =  ch_versions.mix(CLUSTER_EXTRACT.out.versions.first())

    // extract reads from members and references
    ch_members_references
        .map { [ it[0], it[1] ] }
        .set { ch_members }
    ch_members_references
        .map { [ it[0], it[2] ] }
        .set { ch_references }

    // TODO: ch_db need to be combined with ch_members and ch_members must lose it meta data
    SEQKIT_GREP_MEMBERS(ch_members, ch_db)
        .out
        .set { ch_seq_members }

    SEQKIT_GREP_REFERENCES(ch_references, ch_db)
        .out
        .set { ch_seq_references }

    ch_versions =  ch_versions.mix(SEQKIT_GREP_REFERENCES.out.versions.first())

    emit:
    seq_members     = ch_seq_members        // channel: [ [ meta ], [ seq_members.fa] ]
    seq_references  = ch_seq_references     // channel: [ [ meta ], [ seq_references.fa ] ]
    versions        = ch_versions
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_member_ref_channel(LinkedHashMap row) {
    new_meta = row[0].copy()
    members = row[1]
    references = row[2]

    new_meta.sample     = new_meta.id.copy()
    new_meta.cluster    = (members =~ /.*\/.*([0-9]+)_members.txt/)[0][1] // try and extract the correct cluster ID associated to the sample
    new_meta.id = "${new_meta.sample}_${new_meta.cluster}"

    def result = [ new_meta, members, references ]
    return result
}

