// modules
include { BMAP_BBDUK } from '../modules/nf-core/bbmap/bbduk/main'

// Subworkflows
include { FASTQ_FASTQC_UMITOOLS_FASTP       } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC } from '../subworkflows/local/fastq_fastqc_umitools_trimmomatic/main'
include { FASTQ_BOWTIE2_SAMTOOLS            } from '../subworkflows/local/fastq_bowtie2_samtools/main'

workflow PREPROCESSING_ILLUMINA {

    take:
    reads                   // channel: [ [ meta ], [ reads ] ]
    reference               //    file: /path/to.fasta
    index                   //    file: /path/to.index
    adapterlist             //    file: /path/to.adapterlist

    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    // QC & UMI & Trimming with fastp or trimmomatic
    if (trim_tool == 'trimmomatic') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            reads,
            params.skip_fastqc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            params.skip_trimmed_fail,
            params.save_merged,
            params.min_trimmed_reads
            )
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.versions.first())

        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_raw_zip)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.fastqc_trim_html)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.trim_json)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.umi_log)

        ch_reads_trim = FASTQ_FASTQC_UMITOOLS_TRIMMOMATIC.out.reads
    }
    else if (trim_tool == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            reads,
            params.skip_fastqc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            adapterlist,
            params.skip_trimmed_fail,
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
        exit 1
    }

    // Decomplexification with BBDuk
    if (params.skip_complexity_filtering) {
        BMAP_BBDUK ( ch_reads_trim, params.contaminants )
        ch_reads_decomplexified = BMAP_BBDUK.out.reads
        ch_multiqc_files = BMAP_BBDUK.out.log
    } else {
        ch_reads_decomplexified = ch_reads_trim
    }

    // Host removal with Bowtie2
    if (params.skip_hostremoval){
        FASTQ_BOWTIE2_SAMTOOLS ( ch_reads_decomplexified, reference, index )
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

