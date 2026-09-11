include { noBlastHitsToMultiQC  } from '../utils_nfcore_viralmetagenome_pipeline'
include { BLAST_BLASTN          } from '../../../modules/nf-core/blast/blastn/main'
include { BLAST_FILTER          } from '../../../modules/local/blast_filter'

workflow FASTA_BLAST_REFSEL {
    take:
    ch_fasta          // channel: [ val(meta), path(fasta)]
    ch_blacklist      // channel: [ path(blacklist) ]
    ch_blast_db       // channel: [ val(meta), path(db) ]
    ch_blast_db_fasta // channel: [ val(meta), path(fasta) ]

    main:

    // Find a complete genome that's already assembled
    BLAST_BLASTN(
        ch_fasta,
        ch_blast_db,
        [], // taxidlist
        [], // taxids
        []  // negative_tax
    )

    ch_blast_txt = BLAST_BLASTN.out.txt.branch { _meta, txt ->
        no_hits: txt.countLines() == 0
        hits: txt.countLines() > 0
    }

    // Make a table of samples that did not have any blast hits
    ch_no_blast_hits = channel.empty()
    ch_no_blast_hits = ch_blast_txt.no_hits.join(ch_fasta)

    ch_no_blast_hits_mqc = noBlastHitsToMultiQC(ch_no_blast_hits,params.assemblers).collectFile(name:'samples_no_blast_hits_mqc.tsv')

    // Filter out false positve hits that based on query length, alignment length, identity, e-score & bit-score
    ch_input_blast_filter = ch_blast_txt.hits
        .join(ch_fasta, by: [0], remainder: true)
        .multiMap { meta, txt, fasta ->
            hits: [meta, txt ? txt : []]
            contigs: [meta, fasta]
        }

    BLAST_FILTER(
        ch_input_blast_filter.hits,
        ch_input_blast_filter.contigs,
        ch_blacklist,
        ch_blast_db_fasta,
    )

    emit:
    fasta_ref_contigs = BLAST_FILTER.out.sequence // channel: [ val(meta), [ fasta ] ]
    no_blast_hits     = ch_no_blast_hits_mqc      // channel: [ val(meta), [ mqc ] ]
}
