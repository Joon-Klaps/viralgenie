def filterContigs(contig, min_len, n_100) {
    contig
        .map { meta, fasta -> [ meta, fasta, WorkflowCommons.getLengthAndAmbigous( fasta ) ] }
        .branch { meta, fasta, stats ->
            pass: stats.contig_size >= min_len.toInteger() && stats.n_100 <= n_100.toInteger()
                return [ meta, fasta ]
            fail: stats.contig_size < min_len.toInteger() || stats.n_100 > n_100.toInteger()
                return [ meta, fasta, stats ]}
}

def failedContigsToMultiQC(tsv_data, min_len, n_100) {
    tsv_data
        .map { meta, fasta, stats -> ["$meta.id\t$meta.sample\t$meta.cluster_id\t$meta.previous_step\t$stats.contig_size\t$stats.n_100"] }
        .collect()
        .map { tsv ->
            WorkflowCommons.multiqcTsvFromList(
                tsv,
                ['Id','sample name', 'cluster','step','contig size', 'N\'s %'],
                [
                    "id: 'failed_contig_quality'",
                    "anchor: 'WARNING: Filtered contigs'",
                    "section_name: 'Failed contig quality'",
                    "format: 'tsv'",
                    "description: 'Contigs that are not of minimum size ${min_len} or have more then ${n_100} ambigous bases per 100 kbp were filtered out'",
                    "plot_type: 'table'"
                ]
            )
        }
}

def failedMappedReadsToMultiQC(tsv_data, min_mapped_reads) {
    tsv_data
        .map { meta, bam, mapped_reads ->
            ["$meta.id\t$meta.sample\t$meta.cluster_id\t$meta.previous_step\t$mapped_reads"]
            }
        .collect()
        .map { tsv ->
            WorkflowCommons.multiqcTsvFromList(tsv,
                ['Id','sample name', 'cluster','step','mapped reads'],
                [
                    "id: 'failed_mapped'",
                    "anchor: 'WARNING: Filtered contigs'",
                    "section_name: 'Minimum mapped reads'",
                    "format: 'tsv'",
                    "description: 'Contigs that did not have more then ${min_mapped_reads} mapped reads were filtered out'",
                    "plot_type: 'table'"
                ]
            )
        }
}

def noBlastHitsToMultiQC(tsv_data, assemblers) {
    tsv_data
        .map { meta, txt, fasta ->
            def n_fasta = fasta.countFasta()
            ["$meta.sample\t$n_fasta"]}
        .collect()
        .map { tsv ->
            WorkflowCommons.multiqcTsvFromList(tsv,
                ['sample name', "number of contigs"],
                [
                    "id: 'samples_without_blast_hits'",
                    "anchor: 'WARNING: Filtered samples'",
                    "section_name: 'Samples without blast hits'",
                    "format: 'tsv'",
                    "description: 'Samples that did not have any blast hits for their contigs (using ${assemblers}) were not included in further assembly polishing'",
                    "plot_type: 'table'"
                ]
            )
        }
}
