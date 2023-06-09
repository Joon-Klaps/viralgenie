include { BLAST_MAKEBLASTDB } from '../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN      } from '../../modules/nf-core/blast/blastn/main'
include { SEQKIT_GREP       } from '../../modules/nf-core/seqkit/grep/main'
include { CDHIT_CDHIT       } from '../../modules/nf-core/cdhit/cdhit/main'
include { UNTAR             } from '../../modules/nf-core/untar/main'
include { CAT_CAT           } from '../../modules/nf-core/cat/cat/main'

workflow  {

    take:
    fasta    // channel: [ val(meta), [ fasta ] ]
    db       // channel: [ path(db_fasta) ]

    main:
    ch_versions = Channel.empty()


    if (db.endsWith('.gz')) {
                UNTAR (
                    [ [:], db ]
                )
                ch_db = UNTAR.out.untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR.out.versions)
            } else {
                        ch_db = Channel.value(file(db))
            }

    BLAST_MAKEBLASTDB ( ch_db )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions.first())

    BLAST_BLASTN (
        fasta,
        BLAST_MAKEBLASTDB.out.db
    )

    SEQKIT_GREP (ch_db, BLAST_BLASTN.out.txt)
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    fasta.join(
            SEQKIT_GREP.out.filter,
            by: [0],
            remainder: true
            )
        .map {
            it ->
                [it[0],it[1..-1]]
            }
        .set {ch_reference_contigs_comb}

    CDHIT_CDHIT (ch_reference_contigs_comb)
    ch_versions = ch_versions.mix(CDHIT_CDHIT.out.versions.first())


    //TODO: Make a script that extracts the members of the clusters from the cdhit output

    //TODO: Reextract the members from the original fasta file

    // TODO: New subworkflow - Map to cluster representative & create frakenstein

    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())



    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

