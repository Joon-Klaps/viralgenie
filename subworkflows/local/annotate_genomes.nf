
include { BLAST_BLASTN      as BLASTN_ANNOTATION      } from '../../../modules/nf-core/blast/blastn/main'

workflow ANNOTATE_GENOMES {

    take:
    genomes,
    annotation_db

    main:

    BLASTN_ANNOTATION( genomes, annotation_db )


    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

