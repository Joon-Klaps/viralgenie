// modules

include { lowReadSamplesToMultiQC            } from '../../modules/local/functions'
include { PRINSEQPLUSPLUS as PRINSEQ_READS   } from '../../modules/nf-core/prinseqplusplus/main'
include { HUMID                              } from '../../modules/nf-core/humid/main'
include { BBMAP_BBDUK                        } from '../../modules/nf-core/bbmap/bbduk/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC  } from './fastq_fastqc_umitools_trimmomatic'
include { FASTQ_FASTQC_UMITOOLS_FASTP        } from '../nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_KRAKEN_HOST_REMOVE           } from './fastq_kraken_host_remove'

workflow PREPROCESSING_ILLUMINA {

    take:
    ch_reads                   // channel: [ [ meta ], [ ch_reads ] ]
    ch_kraken2_host_db         // channel: [ path(kraken2_host_db) ]
    ch_adapter_fasta           // channel: [ path(adapter_fasta) ]
    ch_contaminants            // channel: [ path(contaminants_fasta) ]

    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()
    trim_read_count     = Channel.empty()

    // QC & UMI & Trimming with fastp or trimmomatic
    if (params.trim_tool == 'trimmomatic') {
        FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC (
            ch_reads,
            params.skip_fastqc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            params.save_trimmed_fail,
            params.save_merged,
            params.min_trimmed_reads
            )
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_raw_zip)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_trim_html)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.trim_log)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.umi_log)

        trim_read_count = trim_read_count.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.trim_read_count)

        ch_reads_trim = FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.reads
    }
    else if (params.trim_tool == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            ch_reads,
            params.skip_fastqc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            ch_adapter_fasta,
            params.save_trimmed_fail,
            params.save_merged,
            params.min_trimmed_reads
            )

        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.umi_log)

        trim_read_count = trim_read_count.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count)

        ch_reads_trim = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
    }
    else {
        throw new Exception("Unknown trim tool: ${trim_tool}")
    }

    // Keeping track of failed reads for reporting
    trim_read_count
        .filter{meta,num_reads -> num_reads < params.min_trimmed_reads.toLong() }
        .set { failed_reads }

    // deduplicate UMI's with HUMID
    if (params.with_umi && ['read', 'both'].contains(params.umi_deduplicate) && params.deduplicate ) {
        HUMID (
            ch_reads_trim,
            [[:],[]]
        )
        ch_reads_dedup   = HUMID.out.dedup
        ch_multiqc_files = ch_multiqc_files.mix(HUMID.out.stats)
        ch_versions      = ch_versions.mix(HUMID.out.versions)
    }
    else {
        ch_reads_dedup = ch_reads_trim
    }


    // Decomplexification with BBDuk
    if (!params.skip_complexity_filtering) {
        if {params.decomplexifier == 'bbduk'} {
            BBMAP_BBDUK (
                ch_reads_dedup,
                ch_contaminants,
                params.decomplexifier
            )
            ch_reads_decomplexified = BBMAP_BBDUK.out.reads
            ch_multiqc_files        = ch_multiqc_files.mix(BBMAP_BBDUK.out.log)
            ch_versions             = ch_versions.mix(BBMAP_BBDUK.out.versions)
        } else if {params.decomplexifier == 'prinseq'} {
            prinseq_in = ch_reads_dedup.map { meta, reads -> [meta, reads, fasta] }
            PRINSEQ_READS (
                ch_reads_dedup,
                params.decomplexifier
            )
            ch_reads_decomplexified = PRINSEQ_READS.out.good_reads
            ch_multiqc_files        = ch_multiqc_files.mix(PRINSEQ_READS.out.log)
            ch_versions             = ch_versions.mix(PRINSEQ_READS.out.versions)
        }
    } else {
        ch_reads_decomplexified = ch_reads_dedup
    }

    // Host removal with kraken2
    if (!params.skip_hostremoval){
        FASTQ_KRAKEN_HOST_REMOVE (
            ch_reads_decomplexified,
            ch_kraken2_host_db,
            params.host_k2_library,
            params.skip_host_fastqc,
            params.min_trimmed_reads,
        )

        ch_reads_hostremoved   = FASTQ_KRAKEN_HOST_REMOVE.out.reads_hostremoved
        failed_reads           = failed_reads.mix(FASTQ_KRAKEN_HOST_REMOVE.out.reads_hostremoved_fail)
        ch_multiqc_files       = ch_multiqc_files.mix( FASTQ_KRAKEN_HOST_REMOVE.out.mqc )
        ch_versions            = ch_versions.mix( FASTQ_KRAKEN_HOST_REMOVE.out.versions )

    } else {
        ch_reads_hostremoved = ch_reads_decomplexified
    }

    //
    // Create a section that reports failed samples and their read counts
    //
    lowReadSamplesToMultiQC(failed_reads, params.min_trimmed_reads)
        .collectFile(name:'samples_low_reads_mqc.tsv')
        .set{low_reads_mqc}

    emit:
    reads                   = ch_reads_hostremoved            // channel: [ [ meta ], [ reads ] ]
    reads_decomplexified    = ch_reads_decomplexified         // channel: [ [ meta ], [ reads ] ]
    reads_trimmed           = ch_reads_trim                   // channel: [ [ meta ], [ reads ] ]
    mqc                     = ch_multiqc_files                // channel: [ [ meta ], [ mqc ] ]
    low_reads_mqc           = low_reads_mqc                   // channel: [ mqc ]
    versions                = ch_versions                     // channel: [ versions.yml ]
}

