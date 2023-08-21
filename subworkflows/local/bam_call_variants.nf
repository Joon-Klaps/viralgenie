include { TABIX_TABIX    } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../../modules/nf-core/bcftools/stats/main'

workflow  {

    take:
    bam             // channel: [ val(meta), [ bam ] ]
    fasta           // channel: [val (meta), [ fasta] ]
    variant_caller  // value: [ bcftools | ivar ]


    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

