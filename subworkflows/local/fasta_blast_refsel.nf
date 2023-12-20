include { BLAST_BLASTN          } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_FILTER          } from '../../modules/local/blast_filter'

workflow FASTA_BLAST_REFSEL {

    take:
    fasta          // channel: [ val(meta), path(fasta)]
    blast_db       // channel: [ val(meta), path(db) ]
    blast_db_fasta // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()
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

    // Throw an warning if no hits are found
    ch_blast_txt
        .hits
        .collect()
        .ifEmpty{ log.warn "No blast hits were found in any samples of the given BLAST database. Consider updating the search parameters or the database: \n ${params.reference_pool} "}

    // Make a table of samples that did not have any blast hits
    ch_no_blast_hits = Channel.empty()
    ch_blast_txt
        .no_hits
        .join(fasta)
        .map { meta, txt, fasta ->
            def n_fasta = fasta.countFasta()
            ["$meta.sample\t$n_fasta"]}
        .collect()
        .map {
            tsv_data ->
                def comments = [
                    "id: 'samples_without_blast_hits'",
                    "anchor: 'WARNING: Filtered samples'",
                    "section_name: 'Samples without blast hits'",
                    "format: 'tsv'",
                    "description: 'Samples that did not have any blast hits for their contigs (using ${params.assemblers}) were not included in further assembly polishing'",
                    "plot_type: 'table'"
                ]
                def header = ['Sample', "Number of contigs"]
                return WorkflowCommons.multiqcTsvFromList(tsv_data, header, comments) // make it compatible with the other mqc files
        }
        .collectFile(name:'samples_no_blast_hits_mqc.tsv')
        .set { ch_no_blast_hits }

    // Filter out false positve hits that based on query length, alignment length, identity, e-score & bit-score
    ch_blast_txt
        .hits
        .join(fasta, by:[0], remainder:false)
        .set{ hits_contigs }

    BLAST_FILTER (
        hits_contigs, 
        blast_db_fasta
    )
    ch_versions = ch_versions.mix(BLAST_FILTER.out.versions.first())

    emit:
    fasta_ref_contigs = BLAST_FILTER.out.sequence  // channel: [ val(meta), [ fasta ] ]
    no_blast_hits     = ch_no_blast_hits           // channel: [ val(meta), [ mqc ] ]
    versions          = ch_versions                // channel: [ versions.yml ]
}

