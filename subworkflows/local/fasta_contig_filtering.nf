//
// Run quast and remove contigs based on the results
//

include { QUAST  } from '../../modules/nf-core/quast/main'

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


    QUAST ( ch_contigs, [[:],[]], [[:],[]] )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    ch_contigs
        .join(QUAST.out.tsv, by:[0])
        .map { meta, fasta, report -> [ meta, fasta, getQuastStats( report ) ] }
        .filter { meta, fasta, stats -> stats.contig_size >= min_len.toInteger() && stats.n_100 <= n_100.toInteger() }
        .set { ch_contigs_filtered_stats }

    ch_contigs_filtered_stats
        .map { meta, fasta, stats -> [ meta, fasta ] }
        .set { ch_contigs_filtered }

    emit:
    contigs       = ch_contigs_filtered           // channel: [ val(meta), [ fasta ] ]
    contigs_stats = ch_contigs_filtered_stats     // channel: [ val(meta), [ fasta] , [ stats ] ]
    versions      = ch_versions                   // channel: [ versions.yml ]
    }

