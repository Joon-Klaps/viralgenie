/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    spades_modes     : ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio']
]

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

def createFileChannel(param) {
    return param ? Channel.fromPath(param, checkIfExists: true) : Channel.empty()
}

def createChannel(dbPath, dbName, skipFlag) {
    return dbPath && skipFlag ? Channel.fromPath(dbPath, checkIfExists: true).map { db -> [[id: dbName], db] } : Channel.empty()
}

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config, params.adapter_fasta, params.contaminants,
    params.spades_yml,params.spades_hmm
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input            ) { ch_input = file(params.input)                                      } else { exit 1, 'Input samplesheet not specified!'                              }

// Optional parameters
ch_adapter_fasta = createFileChannel(params.adapter_fasta)
ch_metadata      = createFileChannel(params.metadata)
ch_contaminants  = createFileChannel(params.contaminants)
ch_spades_yml    = createFileChannel(params.spades_yml)
ch_spades_hmm    = createFileChannel(params.spades_hmm)

def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL & NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Preprocessing
include { PREPROCESSING_ILLUMINA          } from '../subworkflows/local/preprocessing_illumina'
include { UNPACK_DB as UNPACK_DB_KRAKEN2_HOST  } from '../subworkflows/local/unpack_db'

// metagenomic diversity
include { FASTQ_KRAKEN_KAIJU              } from '../subworkflows/local/fastq_kraken_kaiju'
include { UNPACK_DB as UNPACK_DB_KRAKEN   } from '../subworkflows/local/unpack_db'
include { UNPACK_DB as UNPACK_DB_BRACKEN  } from '../subworkflows/local/unpack_db'
include { UNPACK_DB as UNPACK_DB_KAIJU    } from '../subworkflows/local/unpack_db'

// Assembly
include { FASTQ_SPADES_TRINITY_MEGAHIT    } from '../subworkflows/local/fastq_spades_trinity_megahit'

// Consensus polishing of genome
include { FASTA_CONTIG_CLUST              } from '../subworkflows/local/fasta_contig_clust'
include { BLAST_MAKEBLASTDB               } from '../modules/nf-core/blast/makeblastdb/main'
include { ALIGN_COLLAPSE_CONTIGS          } from '../subworkflows/local/align_collapse_contigs'
include { UNPACK_DB as UNPACK_DB_BLAST    } from '../subworkflows/local/unpack_db'
include { FASTQ_FASTA_ITERATIVE_CONSENSUS } from '../subworkflows/local/fastq_fasta_iterative_consensus'
include { SINGLETON_FILTERING             } from '../subworkflows/local/singleton_filtering'

// Variant calling
include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_CONSTRAIN } from '../modules/local/rename_fasta_header'
include { FASTQ_FASTA_MAP_CONSENSUS                            } from '../subworkflows/local/fastq_fasta_map_consensus'

// QC consensus
include { UNPACK_DB as UNPACK_DB_CHECKV   } from '../subworkflows/local/unpack_db'
include { CHECKV_DOWNLOADDATABASE         } from '../modules/nf-core/checkv/downloaddatabase/main'
include { CONSENSUS_QC                    } from '../subworkflows/local/consensus_qc'

// Report generation
include { CREATE_MULTIQC_TABLES           } from '../modules/local/create_multiqc_tables'
include { MULTIQC as MULTIQC_DATAPREP     } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_REPORT       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VIRALGENIE {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Importing samplesheet
    ch_reads = Channel.fromSamplesheet(
        'input'
        ).map{
            meta, read1, read2 ->
            [meta , [read1, read2]]
            }

    // unpack host & contamination DB
    ch_k2_host = Channel.empty()
    if (!params.skip_hostremoval && params.host_k2_db){
        UNPACK_DB_KRAKEN2_HOST(params.host_k2_db)
        ch_k2_host = UNPACK_DB_KRAKEN2_HOST.out.db
    }

    // preprocessing illumina reads (trim - decomplexify - host remove)
    PREPROCESSING_ILLUMINA (
        ch_reads,
        ch_k2_host,
        ch_adapter_fasta,
        ch_contaminants)
    ch_host_trim_reads      = PREPROCESSING_ILLUMINA.out.reads
    ch_decomplex_trim_reads = PREPROCESSING_ILLUMINA.out.reads_decomplexified
    ch_multiqc_files        = ch_multiqc_files.mix(PREPROCESSING_ILLUMINA.out.mqc.collect{it[1]}.ifEmpty([]))
    ch_versions             = ch_versions.mix(PREPROCESSING_ILLUMINA.out.versions)


    // Prepare database for read & contig annotation
    ch_kraken2_db = Channel.empty()
    ch_kaiju_db   = Channel.empty()
    if(!params.skip_metagenomic_diversity || !params.skip_precluster){
        UNPACK_DB_KRAKEN(params.kraken2_db)
        ch_kraken2_db = UNPACK_DB_KRAKEN.out.db
        UNPACK_DB_KAIJU(params.kaiju_db)
        ch_kaiju_db   = UNPACK_DB_KAIJU.out.db
    }

    // Determining metagenomic diversity
    if (!params.skip_metagenomic_diversity) {
        ch_bracken_db   = Channel.empty()
        if (!params.skip_bracken){
            ch_bracken_db = UNPACK_DB_BRACKEN(params.bracken_db).db
        }
        FASTQ_KRAKEN_KAIJU(
            ch_host_trim_reads,
            ch_kraken2_db,
            ch_bracken_db,
            ch_kaiju_db
            )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_KRAKEN_KAIJU.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_versions      = ch_versions.mix(FASTQ_KRAKEN_KAIJU.out.versions)
    }

    // Assembly
    ch_unaligned_raw_contigs   = Channel.empty()
    // Channel for consensus sequences that have been generated across different iteration
    ch_consensus               = Channel.empty()
    // Channel for consensus sequences that have been generated at the LAST iteration
    ch_consensus_results_reads = Channel.empty()
    // Channel for summary table of cluseters to include in mqc report
    ch_clusters_summary        = Channel.empty()

    if (!params.skip_assembly) {
        // run different assemblers and combine contigs
        FASTQ_SPADES_TRINITY_MEGAHIT(
            ch_host_trim_reads,
            assemblers,
            ch_spades_yml,
            ch_spades_hmm)

        ch_versions      = ch_versions.mix(FASTQ_SPADES_TRINITY_MEGAHIT.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_SPADES_TRINITY_MEGAHIT.out.mqc.collect{it[1]}.ifEmpty([]))

        // Filter out empty scaffolds
        FASTQ_SPADES_TRINITY_MEGAHIT
            .out
            .scaffolds
            .branch { meta, scaffolds ->
                pass: scaffolds.countFasta() > 0
                fail: scaffolds.countFasta() == 0
            }
            .set{ch_contigs}

        ch_contigs
            .fail
            .map{ meta, scaffolds ->
                def n_fasta = scaffolds.countFasta()
                ["$meta.sample\t$n_fasta"]
            }
            .collect()
            .map {
                tsv_data ->
                    def comments = [
                        "id: 'samples_without_contigs'",
                        "anchor: 'WARNING: Filtered samples'",
                        "section_name: 'Samples without contigs'",
                        "format: 'tsv'",
                        "description: 'Samples that did not have any contigs (using ${params.assemblers}) were not included in further assembly polishing'",
                        "plot_type: 'table'"
                    ]
                    def header = ['Sample', "Number of contigs"]
                    return WorkflowCommons.multiqcTsvFromList(tsv_data, header, comments) // make it compatible with the other mqc files
            }
            .collectFile(name:'samples_no_contigs_mqc.tsv')
            .set{no_contigs}

        ch_multiqc_files = ch_multiqc_files.mix(no_contigs.ifEmpty([]))

        if (!params.skip_polishing || !params.skip_consensus_qc) {
        UNPACK_DB_BLAST (params.reference_pool)
        UNPACK_DB_BLAST
            .out
            .db
            .map{db -> [[id:"blast_db"], db]}
            .set{ch_ref_pool}
        ch_versions         = ch_versions.mix(UNPACK_DB_BLAST.out.versions)

        BLAST_MAKEBLASTDB ( ch_ref_pool )
        ch_blast_db  = BLAST_MAKEBLASTDB.out.db
        ch_versions  = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
        }

        if (!params.skip_polishing){
            // blast contigs against reference & identify clusters of (contigs & references)
            ch_contigs
                .pass
                .join(ch_host_trim_reads, by: [0], remainder: false)
                .set{ch_contigs_reads}

            FASTA_CONTIG_CLUST (
                ch_contigs_reads,
                ch_blast_db,
                ch_ref_pool,
                ch_kraken2_db,
                ch_kaiju_db
                )
            ch_versions = ch_versions.mix(FASTA_CONTIG_CLUST.out.versions)

            // Split up clusters into singletons and clusters of multiple contigs
            FASTA_CONTIG_CLUST
                .out
                .centroids_members
                .map { meta, centroids, members ->
                    [ meta, centroids, members ]
                }
                .branch { meta, centroids, members ->
                    singletons: meta.cluster_size == 0
                        return [ meta + [step:"singleton"], centroids ]
                    multiple: meta.cluster_size >   0
                        return [ meta + [step:"consensus"], centroids, members ]
                }
                .set{ch_centroids_members}

            ch_clusters_summary    = FASTA_CONTIG_CLUST.out.clusters_summary.collect{it[1]}.ifEmpty([])

            ch_multiqc_files       =  ch_multiqc_files.mix(FASTA_CONTIG_CLUST.out.no_blast_hits_mqc.ifEmpty([]))

            // Align clustered contigs & collapse into a single consensus per cluster
            ALIGN_COLLAPSE_CONTIGS (
                ch_centroids_members.multiple,
                params.contig_align_method
                )
            ch_versions = ch_versions.mix(ALIGN_COLLAPSE_CONTIGS.out.versions)


            SINGLETON_FILTERING (
                ch_centroids_members.singletons,
                params.min_contig_size,
                params.max_n_1OOkbp
                )
            ch_versions = ch_versions.mix(SINGLETON_FILTERING.out.versions)

            if ( !params.skip_singleton_filtering ) {
                ch_singletons = SINGLETON_FILTERING.out.filtered
            } else {
                ch_singletons = SINGLETON_FILTERING.out.renamed
            }

            ALIGN_COLLAPSE_CONTIGS
                .out
                .consensus
                .mix( ch_singletons )
                .set{ ch_consensus }

            ch_unaligned_raw_contigs = ALIGN_COLLAPSE_CONTIGS.out.unaligned_fasta
            ch_unaligned_raw_contigs = ch_unaligned_raw_contigs.mix( SINGLETON_FILTERING.out.renamed )

            // We want the meta from the reference channel to be used downstream as this is our varying factor
            // To do this we combine the channels based on sample
            // Extract the reference meta's and reads
            // Make cartesian product of identified references & reads so all references will be mapped against again.
                ch_decomplex_trim_reads
                    .map{meta, fastq -> [meta.sample,meta, fastq]}
                    .set{ch_reads_tmp}

                ch_consensus
                    .map{meta, fasta -> [meta.sample,meta, fasta]}
                    .set{ch_consensus_tmp}

                ch_consensus_tmp
                    .combine(ch_reads_tmp, by: [0])
                    .map{
                        sample, meta_ref, fasta, meta_reads, fastq -> [meta_ref, fasta, fastq]
                    }
                    .set{ch_consensus_results_reads_intermediate}

            if (!params.skip_iterative_refinement) {
                FASTQ_FASTA_ITERATIVE_CONSENSUS (
                    ch_consensus_results_reads_intermediate,
                    params.iterative_refinement_cycles,
                    params.intermediate_mapper,
                    params.with_umi,
                    params.deduplicate,
                    params.call_intermediate_variants,
                    params.intermediate_variant_caller,
                    params.intermediate_consensus_caller,
                    params.get_intermediate_stats,
                    params.min_mapped_reads,
                    params.min_contig_size,
                    params.max_n_1OOkbp
                )
                ch_consensus               = ch_consensus.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.consensus_allsteps)
                ch_consensus_results_reads = FASTQ_FASTA_ITERATIVE_CONSENSUS.out.consensus_reads
                ch_versions                = ch_versions.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.versions)
                ch_multiqc_files           = ch_multiqc_files.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.mqc.ifEmpty([])) //collect already done in subworkflow
            } else {
                ch_consensus_results_reads = ch_consensus_results_reads_intermediate
            }
        }
    }

    // add last step to it
    ch_consensus_results_reads
        .map{ meta, fasta, fastq ->
            [meta + [step: "variant-calling", iteration:'variant-calling', previous_step: meta.step], fasta, fastq]
            }
        .set{ch_consensus_results_reads}

    if (params.mapping_sequence && !params.skip_variant_calling ) {
        ch_mapping_sequence = Channel.fromPath(params.mapping_sequence)

        //get header names of sequences
        ch_mapping_sequence
            .splitFasta(record: [id: true])
            .set {seq_names}

        //split up multifasta
        ch_mapping_sequence
            .splitFasta(file : true, by:1)
            .merge(seq_names)
            .set{ch_map_seq_anno}

        //combine with all input samples & rename the meta.id
        // TODO: consider adding another column to the samplesheet with the mapping sequence
        ch_host_trim_reads
            .combine( ch_map_seq_anno )
            .map
                {
                    meta, reads, seq, name ->
                    sample = meta.id
                    id = "${sample}_${name.id}"
                    new_meta = meta + [
                        id: id,
                        sample: sample,
                        cluster_id: "${name.id}",
                        step: "constrain",
                        constrain: true,
                        reads: reads,
                        iteration: 'variant-calling',
                        previous_step: 'constrain'
                        ]
                    return [new_meta, seq]
                }
            .set{ch_map_seq_anno_combined}

        // For QC we keep original sequence to compare to
        ch_unaligned_raw_contigs = ch_unaligned_raw_contigs.mix(ch_map_seq_anno_combined)

        //rename fasta headers
        RENAME_FASTA_HEADER_CONSTRAIN (ch_map_seq_anno_combined,[])
        ch_versions = ch_versions.mix(RENAME_FASTA_HEADER_CONSTRAIN.out.versions)

        RENAME_FASTA_HEADER_CONSTRAIN
            .out
            .fasta
            .map{ meta, fasta -> [meta, fasta, meta.reads] }
            .set{constrain_consensus_reads}

        //Add to the consensus channel, the mapping sequences will now always be mapped against
        ch_consensus_results_reads = ch_consensus_results_reads.mix(constrain_consensus_reads)
    }

    // After consensus sequences have been made, we still have to map against it and call variants
    if ( !params.skip_variant_calling ) {

        FASTQ_FASTA_MAP_CONSENSUS(
            ch_consensus_results_reads,
            params.mapper,
            params.with_umi,
            params.deduplicate,
            true,
            params.variant_caller,
            true,
            params.consensus_caller,
            params.get_stats,
            params.min_mapped_reads,
            params.min_contig_size,
            params.max_n_1OOkbp
        )
        ch_consensus     = ch_consensus.mix(FASTQ_FASTA_MAP_CONSENSUS.out.consensus_all)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTA_MAP_CONSENSUS.out.mqc.ifEmpty([])) // collect already done in subworkflow

    }

    ch_checkv_summary = Channel.empty()
    ch_quast_summary  = Channel.empty()
    ch_blast_summary  = Channel.empty()

    if ( !params.skip_consensus_qc || (params.skip_assembly && params.skip_variant_calling) ) {
        ch_checkv_db = Channel.empty()
        if (!params.skip_checkv ) {
            if (params.checkv_db){
                ch_checkv_db = UNPACK_DB_CHECKV(params.checkv_db)
            } else {
                ch_checkv_db = CHECKV_DOWNLOADDATABASE().checkv_db
            }
        }

        CONSENSUS_QC(
            ch_consensus,
            ch_unaligned_raw_contigs,
            ch_checkv_db,
            ch_blast_db,
            params.skip_checkv,
            params.skip_quast,
            params.skip_blast_qc,
            params.skip_alignment_qc
            )
        ch_versions       = ch_versions.mix(CONSENSUS_QC.out.versions)
        ch_multiqc_files  = ch_multiqc_files.mix(CONSENSUS_QC.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_checkv_summary = CONSENSUS_QC.out.checkv_summary.collect{it[1]}.ifEmpty([])
        ch_quast_summary  = CONSENSUS_QC.out.quast_summary.collect{it[1]}.ifEmpty([])
        ch_blast_summary  = CONSENSUS_QC.out.blast_txt.collect{it[1]}.ifEmpty([])

    }

    MULTIQC_DATAPREP (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
    )

    multiqc_data = MULTIQC_DATAPREP.out.data.ifEmpty([])

    // Prepare MULTIQC custom tables
    CREATE_MULTIQC_TABLES (
            ch_clusters_summary,
            ch_metadata,
            ch_checkv_summary,
            ch_quast_summary,
            ch_blast_summary,
            multiqc_data
            )
    ch_multiqc_files = ch_multiqc_files.mix(CREATE_MULTIQC_TABLES.out.summary_clusters_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CREATE_MULTIQC_TABLES.out.sample_metadata_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CREATE_MULTIQC_TABLES.out.contigs_overview_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CREATE_MULTIQC_TABLES.out.summary_checkv_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CREATE_MULTIQC_TABLES.out.summary_blast_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CREATE_MULTIQC_TABLES.out.summary_quast_mqc.ifEmpty([]))
    ch_versions      = ch_versions.mix(CREATE_MULTIQC_TABLES.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowViralgenie.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowViralgenie.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)


    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))

    MULTIQC_REPORT (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC_REPORT.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
