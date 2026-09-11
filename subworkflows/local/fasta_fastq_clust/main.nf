include { CDHIT_CDHITEST       } from '../../../modules/nf-core/cdhit/cdhitest/main'
include { VSEARCH_CLUSTER      } from '../../../modules/nf-core/vsearch/cluster/main'
include { MMSEQS_FASTA_CLUSTER } from '../../../subworkflows/nf-core/mmseqs_fasta_cluster/main'
include { VRHYME_VRHYME        } from '../../../modules/nf-core/vrhyme/vrhyme/main'
include { MASH_DIST            } from '../../../modules/nf-core/mash/dist/main'
include { CLUSTY               } from '../../../modules/nf-core/clusty/main'

workflow FASTA_FASTQ_CLUST {
    take:
    ch_fasta_fastq     // channel: [ val(meta), [ fasta ], [ fastq ] ]
    cluster_method     // string
    identity_threshold // value: 0.85

    main:
    ch_fasta = ch_fasta_fastq.map { meta, fasta, _fastq -> [meta, fasta] }

    // cluster our reference hits and contigs should make this a subworkflow
    if (cluster_method == "vsearch") {
        VSEARCH_CLUSTER(ch_fasta)
        ch_clusters = VSEARCH_CLUSTER.out.uc
    }
    else if (cluster_method == "cdhitest") {
        if (identity_threshold < 0.80) {
            log.warn("cdhitest identity threshold is set to ${identity_threshold}, which is below the minimum threshold of 0.80.\n AUTOFIX: Defaulting to 0.80")
        }
        CDHIT_CDHITEST(ch_fasta)

        ch_clusters = CDHIT_CDHITEST.out.clusters
    }
    else if (cluster_method == "mmseqs-linclust" || cluster_method == "mmseqs-cluster") {
        MMSEQS_FASTA_CLUSTER(ch_fasta, cluster_method == "mmseqs-linclust" ? "linclust" : "cluster")
        ch_clusters = MMSEQS_FASTA_CLUSTER.out.clusters
    }
    else if (cluster_method == "vrhyme") {
        vhryme_in = ch_fasta_fastq.branch { meta, fasta, reads ->
            fasta: [meta, fasta]
            reads: [meta, reads]
        }
        VRHYME_VRHYME(vhryme_in.reads, vhryme_in.fasta)
        ch_clusters = VRHYME_VRHYME.out.membership
    }
    else if (cluster_method == "mash") {

        // Calculate distances
        MASH_DIST(ch_fasta)

        // Fix bug with clusty not accepting singletons
        ch_dist = MASH_DIST.out.dist
            .branch{ _meta, dist ->
                def line_count = dist.withInputStream { stream -> stream.readLines().take(3).size() }
                multiple: line_count > 2
                single: line_count <= 2 // if only two lines, then singleton cluster, no need to cluster
            }

        // Determine clusters from distances
        CLUSTY(ch_dist.multiple, [[:], []])
        ch_clusters = ch_dist.single.mix(CLUSTY.out.assignments)
    }

    emit:
    clusters = ch_clusters // channel: [ [ meta ], [ clusters ] ]
}
