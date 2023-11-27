//
// Run quast and remove contigs based on the results
//

include { QUAST as QUAST_FILTER  } from '../../modules/nf-core/quast/main'

// Function to extract the contig size & N' per 100 kbp from the QUAST report
def getQuastStats(report_file) {
    def contig_size = 0
    def n_100 = 100000 // assume the worst
    report_file.eachLine { line ->
        def contig_size_match = line =~ /Largest contig\s([\d]+)/
        def n_100_match       = line =~ /N's per 100 kbp\s([\d\.]+)/
        if (contig_size_match) contig_size = contig_size_match[0][1].toFloat()
        if (n_100_match) n_100 = n_100_match[0][1].toFloat()
    }
    return [contig_size : contig_size, n_100 : n_100]
}

workflow FASTA_CONTIG_FILTERING {

    take:
    ch_contigs // channel: [ val(meta), [ fasta ] ]
    min_len    // integer: min_length
    n_100      // integer: n_100

    main:

    ch_versions = Channel.empty()


    QUAST_FILTER ( ch_contigs, [[:],[]], [[:],[]] )
    ch_versions = ch_versions.mix(QUAST_FILTER.out.versions.first())

    ch_contigs
        .join(QUAST_FILTER.out.tsv, by:[0])
        .map { meta, fasta, report -> [ meta, fasta, getQuastStats( report ) ] }
        .branch { meta, fasta, stats ->
            pass: stats.contig_size >= min_len.toInteger() && stats.n_100 <= n_100.toInteger()
                return [ meta, fasta ]
            fail: stats.contig_size < min_len.toInteger() || stats.n_100 > n_100.toInteger()
                return [ meta, fasta, stats ]}
        .set { ch_contigs_filtered }

    ch_contigs_filtered
        .fail
        .map { meta, fasta, stats ->
            ["$meta.id\t$meta.sample\t$meta.cluster_id\t$meta.step\t$stats.contig_size\t$stats.n_100"]
        }
        .collect()
        .map {
            tsv_data ->
                def comments = [
                    "id: 'failed_contig_quality'",
                    "anchor: 'WARNING: Filtered contigs'",
                    "section_name: 'Failed contig quality'",
                    "format: 'tsv'",
                    "description: 'Contigs that are not of minimum size ${min_len} or have more then ${n_100} ambigous bases per 100 kbp were filtered out'",
                    "plot_type: 'table'"
                ]
                def header = ['Id','Sample', 'Cluster','Step','Contig size', 'N\'s per 100 kbp']
                return WorkflowCommons.multiqcTsvFromList(tsv_data, header, comments) // make it compatible with other mqc files
        }
        .collectFile(name:'failed_contig_quality_mqc.tsv')
        .set { ch_contig_qc_fail_mqc }



    emit:
    contigs            = ch_contigs_filtered.pass      // channel: [ val(meta), [ fasta ] ]
    contig_qc_fail_mqc = ch_contig_qc_fail_mqc         // channel: [ tsv ]
    versions           = ch_versions                   // channel: [ versions.yml ]
}

