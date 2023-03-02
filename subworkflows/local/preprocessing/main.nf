//
// Perform QC - read trimming - QC - complexity filtering - QC - host removal - QC.
//

//modules
include { QC as RAW_QC                 } from '../qc/main'
include { QC as TRIM_QC                } from '../qc/main'
include { QC as COMPLEXITYFILTER_QC    } from '../qc/main'
include { QC as HOSTFILTER_QC          } from '../qc/main'

//subworkflows
include { TRIMMING                 } from '../trimming/main'
include { COMPLEXITYFILTERING      } from '../complexity_filtering/main'
include { HOSTREMOVAL              } from '../hostremoval/main'

// TODO: TEST this subworkflow
workflow PREPROCESSING_ILLUMINA {

    take:
    reads                   // channel: [ [ meta ], [ reads ] ]
    reference               // reference: /path/to.fasta
    index                   // /path/to.index
    qc_tool                 // value 'falco' or 'fastqc'
    trim_tool               // value 'trimmomatic' or 'fastp'
    complexityfilter_tool   // value 'bbduk' or 'prinseqplusplus' or 'fastp'
    adapterlist             // file

    main:

    ch_versions         = Channel.empty()

    ch_multiqc_files    = Channel.empty()

    // QC pre trimming & complexity filtering
    RAW_QC ( reads, qc_tool)
    ch_versions      = ch_versions.mix( RAW_QC.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( RAW_QC.out.mqc )


    // TODO: Decide whether to include fastp complexity filtering as an option as well.
    // Trimming step
    if (params.skip_trimming) {
        TRIMMING( reads, params.trim_tool,adapterlist)
        ch_reads_trimmed   = TRIMMING.out.reads

        ch_multiqc_files   = ch_multiqc_files.mix(TRIMMING.out.mqc )
        ch_versions        = ch_versions.mix(TRIMMING.out.versions)

        // QC subworkflow step post trimming
        TRIM_QC ( trimmed_reads, qc_tool)

        ch_multiqc_files   = ch_multiqc_files.mix( TRIM_QC.out.mqc )
        ch_versions        = ch_versions.mix( TRIM_QC.out.versions )

    } else {
        ch_reads_trimmed = reads
    }

    // Complexity filtering subworkflow step
    if (params.skip_complexityfilter) {
        COMPLEXITYFILTERING ( ch_reads_trimmed , complexityfilter_tool)
        ch_reads_decomplexified    = COMPLEXITYFILTERING.out.reads

        ch_multiqc_files           = ch_multiqc_files.mix(COMPLEXITYFILTERING.out.mqc)
        ch_versions                = ch_versions.mix( COMPLEXITYFILTERING.out.versions )

        // QC subworkflow step post complexity filtering
        COMPLEXITYFILTER_QC ( reads, qc_tool)

        ch_multiqc_files           = ch_multiqc_files.mix( COMPLEXITYFILTER_QC.out.mqc )
        ch_versions                = ch_versions.mix( COMPLEXITYFILTER_QC.out.versions )

    } else {
        ch_reads_decomplexified = ch_reads_trimmed
    }

    // add host removal filtering step post complexity filtering
    if (params.skip_hostremoval){
        HOSTREMOVAL_BOWTIE2_SAMTOOLS ( ch_reads_decomplexified, reference, index )
        ch_reads_hostremoved   = HOSTREMOVAL_BOWTIE2_SAMTOOLS.out.reads

        ch_multiqc_files       = ch_multiqc_files.mix( HOSTREMOVAL_BOWTIE2_SAMTOOLS.out.mqc )
        ch_versions            = ch_versions.mix(HOSTREMOVAL_BOWTIE2_SAMTOOLS.out.versions)

        // QC subworkflow step post host removal filtering
        HOSTFILTER_QC ( reads, qc_tool)

        ch_multiqc_files       = ch_multiqc_files.mix( HOSTFILTER_QC.out.mqc )
        ch_versions            = ch_versions.mix( HOSTFILTER_QC.out.versions )


    } else {
        ch_reads_hostremoved = ch_reads_decomplexified
    }

    emit:
    reads                   = ch_reads_hostremoved            // channel: [ [ meta ], [ reads ] ]
    reads_decomplexified    = ch_reads_decomplexified         // channel: [ [ meta ], [ reads ] ]
    reads_trimmed           = ch_reads_trimmed                // channel: [ [ meta ], [ reads ] ]
    mqc                     = ch_multiqc_files
    versions                = ch_versions                     // channel: [ versions.yml ]
}

