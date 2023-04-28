// modules
include { BBMAP_BBDUK } from '../../modules/nf-core/bbmap/bbduk/main'

// Subworkflows
// > local
include { FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC } from './fastq_fastqc_umitools_trimmomatic'
include { FASTQ_BOWTIE2_SAMTOOLS            } from './fastq_bowtie2_samtools'
// > nf-core
include { FASTQ_FASTQC_UMITOOLS_FASTP       } from '../nf-core/fastq_fastqc_umitools_fastp/main'

workflow PREPROCESSING_ILLUMINA {

    take:
    ch_reads                   // channel: [ [ meta ], [ ch_reads ] ]
    ch_host                    // channel: [ path(host_fasta) ]
    ch_index                   // channel: [ path(index) ]
    ch_adapter_fasta           // channel: [ path(adapter_fasta) ]
    ch_contaminants            // channel: [ path(contaminants_fasta) ]

    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

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
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.versions.first())

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_raw_zip)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_trim_html)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.trim_log)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.umi_log)

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

        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions.first())

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_html)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.umi_log)

        ch_reads_trim = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
    }
    else {
        throw new Exception("Unknown trim tool: ${trim_tool}")
    }

    // Decomplexification with BBDuk
    if (!params.skip_complexity_filtering) {
        BBMAP_BBDUK ( ch_reads_trim, ch_contaminants )
        ch_reads_decomplexified = BBMAP_BBDUK.out.reads
        ch_multiqc_files = BBMAP_BBDUK.out.log
    } else {
        ch_reads_decomplexified = ch_reads_trim
    }

    // Host removal with Bowtie2
    if (!params.skip_hostremoval){
        FASTQ_BOWTIE2_SAMTOOLS ( ch_reads_decomplexified, ch_host, ch_index )
        ch_reads_hostremoved   = FASTQ_BOWTIE2_SAMTOOLS.out.reads

        ch_multiqc_files       = ch_multiqc_files.mix( FASTQ_BOWTIE2_SAMTOOLS.out.mqc )
        ch_versions            = ch_versions.mix(FASTQ_BOWTIE2_SAMTOOLS.out.versions)

    } else {
        ch_reads_hostremoved = ch_reads_decomplexified
    }

    emit:
    reads                   = ch_reads_hostremoved            // channel: [ [ meta ], [ reads ] ]
    reads_decomplexified    = ch_reads_decomplexified         // channel: [ [ meta ], [ reads ] ]
    reads_trimmed           = ch_reads_trim                   // channel: [ [ meta ], [ reads ] ]
    mqc                     = ch_multiqc_files
    versions                = ch_versions                     // channel: [ versions.yml ]
}
