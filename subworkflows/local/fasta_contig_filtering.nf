workflow FASTA_CONTIG_FILTERING {

    take:
    ch_contigs // channel: [ val(meta), [ fasta ] ]
    min_len    // integer: min_length
    n_100      // integer: n_100

    main:

    ch_versions = Channel.empty()

    ch_contigs
        .map { meta, fasta -> [ meta, fasta, WorkflowCommons.getLengthAndAmbigous( fasta ) ] }
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
                def header = ['Id','Sample', 'Cluster','Step','Contig size', 'N\'s %']
                return WorkflowCommons.multiqcTsvFromList(tsv_data, header, comments) // make it compatible with other mqc files
        }
        .collectFile(name:'failed_contig_quality_mqc.tsv')
        .set { ch_contig_qc_fail_mqc }

    ch_contigs_filtered.pass.view()

    emit:
    contigs            = ch_contigs_filtered.pass      // channel: [ val(meta), [ fasta ] ]
    contig_qc_fail_mqc = ch_contig_qc_fail_mqc         // channel: [ tsv ]
    versions           = ch_versions                   // channel: [ versions.yml ]
}

