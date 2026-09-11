include { KRAKEN2_KRAKEN2 as KRAKEN2_HOST_REMOVE } from '../../../modules/nf-core/kraken2/kraken2/main'
include { FASTQC as FASTQC_HOST                  } from '../../../modules/nf-core/fastqc/main'

workflow FASTQ_KRAKEN_HOST_REMOVE {
    take:
    ch_reads
    ch_kraken2_host_db
    skip_fastqc
    min_reads

    main:
    ch_multiqc_files = channel.empty()

    // remove host reads & keep unclassified reads [true, true]
    KRAKEN2_HOST_REMOVE(
        ch_reads,
        ch_kraken2_host_db,
        true,
        true,
    )

    ch_multiqc_files     = ch_multiqc_files.mix(KRAKEN2_HOST_REMOVE.out.report)
    ch_reads_hostremoved = KRAKEN2_HOST_REMOVE.out.unclassified_reads_fastq
        .join(KRAKEN2_HOST_REMOVE.out.report)
        .map { meta, fastq, tsv -> [meta, fastq, getReadsAfterHostRemove(tsv)] }
        .branch { meta, fastq, n_reads ->
            pass: n_reads > min_reads
            return [meta, fastq]
            fail: n_reads <= min_reads
            return [meta, n_reads]
        }

    // fastqc
    if (!skip_fastqc) {
        FASTQC_HOST(
            ch_reads_hostremoved.pass
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HOST.out.html)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_HOST.out.zip)
    }

    emit:
    reads_hostremoved      = ch_reads_hostremoved.pass // channel: [ [ meta ], [ fastq ] ]
    reads_hostremoved_fail = ch_reads_hostremoved.fail // channel: [ [ meta ], [ n_reads ] ]
    mqc                    = ch_multiqc_files          // channel: [ multiqc_files ]
}

def getReadsAfterHostRemove(tsv) {
    def n_reads = 0
    def lines = tsv.readLines()
    if (!lines) return 0
    def firstLine = lines.first()
    if (firstLine =~ /(unclassified)/) {
        def valuePart = firstLine.split('\t')[1]
        // Fetching the value from the second column
        n_reads += valuePart.toLong()
    }
    return n_reads.toLong()
}
