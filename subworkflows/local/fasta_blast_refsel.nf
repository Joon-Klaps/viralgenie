include { noBlastHitsToMultiQC  } from '../../modules/local/functions'
include { BLAST_BLASTN          } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_FILTER          } from '../../modules/local/blast_filter'

workflow FASTA_BLAST_REFSEL {

    take:
    fasta_fastq    // channel: [ val(meta), path(fasta), path(fastq)]
    blast_db       // channel: [ val(meta), path(db) ]
    blast_db_fasta // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()
    fasta = fasta_fastq.map{ meta, fasta, fastq -> [meta, fasta] }

    // Blast results, to a reference database, to find a complete genome that's already assembled
    BLAST_BLASTN (
        fasta,
        blast_db
    )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())

    BLAST_BLASTN
        .out
        .txt
        .branch { meta, txt ->
            no_hits : txt.countLines() == 0
            hits    : txt.countLines() > 0
        }
        .set { ch_blast_txt }

    // Make a table of samples that did not have any blast hits
    no_blast_hits_mqc = Channel.empty()
    ch_blast_txt
        .no_hits
        .join(fasta)
        .set{no_blast_hits}

    no_blast_hits_mqc = noBlastHitsToMultiQC(no_blast_hits,params.assemblers).collectFile(name:'samples_no_blast_hits_mqc.tsv')

    // Filter out false positve hits that based on query length, alignment length, identity, e-score & bit-score
    ch_blast_txt
        .hits
        .join(fasta, by:[0], remainder:true)
        .join(blast_db_fasta, by:[0], remainder:true)
        .multiMap{
            meta, txt, fasta, db_fasta ->
            hits : [meta, txt ? txt : []]
            contigs : [meta, fasta]
            db_fasta : [meta, db_fasta]
        }
        .set{input_blast_filter}

    BLAST_FILTER (
        input_blast_filter.hits,
        input_blast_filter.contigs,
        input_blast_filter.db_fasta
    )
    ch_versions = ch_versions.mix(BLAST_FILTER.out.versions.first())

    ch_fasta_sel_fastq = BLAST_FILTER.out.sequence.join(fasta_fastq, by:[0])

    emit:
    fasta_sel_fastq   = ch_fasta_sel_fastq  // channel: [ val(meta), [ fasta_all ], [contigs], [ fastq ] ]
    no_blast_hits     = no_blast_hits_mqc   // channel: [ val(meta), [ mqc ] ]
    versions          = ch_versions         // channel: [ versions.yml ]
}

