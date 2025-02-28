include { BWAMEM2_MEM       } from '../../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX     } from '../../modules/nf-core/bwamem2/index/main'
include { BOWTIE2_ALIGN     } from '../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD     } from '../../modules/nf-core/bowtie2/build/main'
include { BWA_MEM           } from '../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX         } from '../../modules/nf-core/bwa/index/main'

workflow MAP_READS  {

    take:
    reference_reads     // channel: [ val(meta), [ fasta ], [ reads ] ]
    mapper              // val: 'bwamem2' or 'bowtie2' or 'bwa'

    main:

    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    ch_reads    = reference_reads.map{meta, fasta,fastq -> [ meta, fastq ]}
    reference   = reference_reads.map{meta, fasta,fastq -> [ meta, fasta ]}

    if ( mapper == 'bwamem2' ) {
        BWAMEM2_INDEX ( reference )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

        bwamem2_input = reference_reads
            .join(BWAMEM2_INDEX.out.index, by: [0])
            .multiMap{meta, fasta, fastq, index ->
                reads: [ meta, fastq]
                index: [ meta, index ]
                fasta: [ meta, fasta ]
            }

        BWAMEM2_MEM ( bwamem2_input.reads, bwamem2_input.index, bwamem2_input.fasta, true )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
        //no mqc for bwamem2

        ch_bam      = BWAMEM2_MEM.out.bam
    }
    else if ( mapper == 'bowtie2' ) {
        BOWTIE2_BUILD ( reference )
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions.first())

        bowtie2_input = reference_reads
            .join(BOWTIE2_BUILD.out.index, by: [0])
            .multiMap{meta, fasta, fastq, index ->
                reads: [ meta, fastq]
                index: [ meta, index ]
                fasta: [ meta, fasta ]
            }

        BOWTIE2_ALIGN ( bowtie2_input.reads, bowtie2_input.index, bowtie2_input.fasta, false, true)
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())

        ch_bam      = BOWTIE2_ALIGN.out.bam
        ch_multiqc  = ch_multiqc.mix(BOWTIE2_ALIGN.out.log)
    } else if ( mapper == "bwa") {
        BWA_INDEX ( reference )
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())
        ch_reads
            .join(BWA_INDEX.out.index, by: [0])
            .join(reference, by: [0])
            .multiMap{meta, reads, index, fasta ->
                reads: [ meta, reads ]
                index: [ meta, index ]
                fasta: [ meta, fasta ]
            }
            .set { bwamem_in }

        BWA_MEM ( bwamem_in.reads, bwamem_in.index, bwamem_in.fasta, true )
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
        //no mqc for bwa

        ch_bam      = BWA_MEM.out.bam

    } else {
        error ("Unknown mapper: ${mapper}")
    }

    emit:
    bam      = ch_bam      // channel: [ val(meta), [ bam ] ]
    ref      = reference   // channel: [ val(meta), [ fasta ] ]
    mqc      = ch_multiqc  // channel: [ multiqc ]

    versions = ch_versions // channel: [ versions.yml ]
}

