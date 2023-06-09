include { GUNZIP            } from '../../modules/nf-core/gunzip/main'
include { BLAST_BLASTN      } from '../../modules/nf-core/blast/blastn/main'
include { SEQKIT_GREP       } from '../../modules/nf-core/seqkit/grep/main'
include { CDHIT_CDHIT       } from '../../modules/nf-core/cdhit/cdhit/main'
include { CAT_CAT           } from '../../modules/nf-core/cat/cat/main'

workflow FASTA_BLAST_CDHIT {

    take:
    fasta          // channel: [ val(meta), [ fasta.gz ] ]
    blast_db       // channel: [ path(db) ]
    blast_db_fasta // channel: [path(db) ]

    main:
    ch_versions = Channel.empty()

    GUNZIP(fasta)
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())
    ch_fasta_gunzip = GUNZIP.out.gunzip

    // Blast results, to a reference database, to find a complete genome that's already assembled
    BLAST_BLASTN (
        ch_fasta_gunzip,
        blast_db
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    // give the result meta again so we can use it again
    blast_db_fasta
        .combine(BLAST_BLASTN.out.txt)
        .map{
            it -> [it[1],it[0]]
        }
        .set{ch_blast_db_fasta}

    ch_hits = BLAST_BLASTN.out.txt.map{it -> it[1]}
    // isolate the hits from the database to a fasta file
    SEQKIT_GREP (ch_blast_db_fasta, ch_hits)
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    // put refernece hits and contigs together
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

    ch_reference_contigs_comb.view()

    CAT_CAT(ch_reference_contigs_comb)

    // cluster our reference hits and contigs
    CDHIT_CDHIT (CAT_CAT.out.file_out)
    ch_versions = ch_versions.mix(CDHIT_CDHIT.out.versions.first())


    //TODO: Make a script that extracts the members of the clusters from the cdhit output

    //TODO: Reextract the members from the original fasta file

    // TODO: New subworkflow - Map to cluster representative & create frakenstein



    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}

