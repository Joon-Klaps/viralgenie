def filterContigs(contig, min_len, n_100) {
    contig
    // //     .branch { meta, fasta, stats ->
    // //         pass: return [ meta, fasta ]
    // //         fail: return [ meta, fasta, stats ]}
    // contig
    //     .map { meta, fasta -> [ meta, fasta, WorkflowCommons.getLengthAndAmbigous( fasta ) ] }
    //     .view { it }
    //     .branch { meta, fasta, stats ->
    //         pass: stats.contig_size >= min_len.toInteger() && stats.n_100 <= n_100.toInteger()
    //             return [ meta, fasta ]
    //         fail: stats.contig_size < min_len.toInteger() || stats.n_100 > n_100.toInteger()
    //             return [ meta, fasta, stats ]}
}

def failedContigsToMultiQC(tsv_data, min_len, n_100) {
    WorkflowCommons.multiqcTsvFromList(
        tsv_data,
        [
            "id: 'failed_contig_quality'",
            "anchor: 'WARNING: Filtered contigs'",
            "section_name: 'Failed contig quality'",
            "format: 'tsv'",
            "description: 'Contigs that are not of minimum size ${min_len} or have more then ${n_100} ambigous bases per 100 kbp were filtered out'",
            "plot_type: 'table'"
        ],
        ['Id','Sample', 'Cluster','Step','Contig size', 'N\'s %']
    )
}

def failedMappedReadsToMultiQC(tsv_data, min_mapped_reads) {
    WorkflowCommons.multiqcTsvFromList(tsv_data,
        [
            "id: 'failed_mapped'",
            "anchor: 'WARNING: Filtered contigs'",
            "section_name: 'Minimum mapped reads'",
            "format: 'tsv'",
            "description: 'Contigs that did not have more then ${min_mapped_reads} mapped reads were filtered out'",
            "plot_type: 'table'"
        ],
        ['Id','Sample', 'Cluster','Step','Mapped reads']
    )
}

def noBlastHitsToMultiQC(tsv_data, assemblers) {
    WorkflowCommons.multiqcTsvFromList(tsv_data,
        [
            "id: 'samples_without_blast_hits'",
            "anchor: 'WARNING: Filtered samples'",
            "section_name: 'Samples without blast hits'",
            "format: 'tsv'",
            "description: 'Samples that did not have any blast hits for their contigs (using ${assemblers}) were not included in further assembly polishing'",
            "plot_type: 'table'"
        ],
        ['Sample', "Number of contigs"]
    )
}
