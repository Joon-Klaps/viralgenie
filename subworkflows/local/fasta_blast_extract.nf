include { BLAST_BLASTN      } from '../../modules/nf-core/blast/blastn/main'
include { BLAST_FILTER      } from '../../modules/local/blast_filter'
include { GUNZIP            } from '../../modules/nf-core/gunzip/main'
include { SEQKIT_GREP       } from '../../modules/nf-core/seqkit/grep/main'
include { CAT_CAT           } from '../../modules/nf-core/cat/cat/main'

workflow FASTA_BLAST_EXTRACT {

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
    BLAST_FILTER (
        ch_blast_txt.hits
    )
    ch_versions = ch_versions.mix(BLAST_FILTER.out.versions.first())

    // give the references a meta again so it can be used in seqkit
    blast_db_fasta
        .combine(BLAST_FILTER.out.hits)
        .map{ meta_blast, db_fasta, meta_filter, filter -> [meta_filter, db_fasta]
        }
        .set{ch_blast_db_fasta}

    ch_hits = BLAST_FILTER.out.hits.map{it -> it[1]}

    // isolate the hits from the database to a fasta file
    SEQKIT_GREP (ch_blast_db_fasta, ch_hits)
    ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions.first())

    GUNZIP(SEQKIT_GREP.out.filter)
    ch_versions = ch_versions.mix(GUNZIP.out.versions.first())

    // put reference hits and contigs together
    fasta
        .mix(GUNZIP.out.gunzip)
        .groupTuple(size :2, remainder: false) // group and throw away those that don't fit
        .set {ch_reference_contigs_comb}

    CAT_CAT(ch_reference_contigs_comb)

    emit:
    fasta_ref_contigs = CAT_CAT.out.file_out     // channel: [ val(meta), [ bam ] ]
    no_blast_hits     = ch_no_blast_hits         // channel: [ val(meta), [ bai ] ]

    versions          = ch_versions              // channel: [ versions.yml ]
}

