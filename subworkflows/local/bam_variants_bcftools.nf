include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM    } from '../../modules/nf-core/bcftools/norm/main'

workflow  {

    take:
    // TODO nf-core: edit input (take) channels
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:
 
    ch_versions = Channel.empty()

    //
    // Call variants
    //
    BCFTOOLS_MPILEUP (
        bam.map{ meta, bam_file -> [ meta, bam_file, [] ] },
        fasta,
        params.save_mpileup
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

    // Filter out samples with 0 variants
    BCFTOOLS_MPILEUP
        .out
        .vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .join(BCFTOOLS_MPILEUP.out.stats)
        .filter { meta, vcf, tbi, stats -> WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0 }
        .set { ch_vcf_tbi_stats }

    ch_vcf_tbi_stats
        .map { meta, vcf, tbi, stats -> [ meta, vcf ] }
        .set { ch_vcf }

    ch_vcf_tbi_stats
        .map { meta, vcf, tbi, stats -> [ meta, tbi ] }
        .set { ch_tbi }

    ch_vcf_tbi_stats
        .map { meta, vcf, tbi, stats -> [ meta, stats ] }
        .set { ch_stats }

    //
    // Split multi-allelic positions
    //
    BCFTOOLS_NORM (
        ch_vcf.join(ch_tbi, by: [0]),
        fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())


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

