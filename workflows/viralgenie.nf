/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    trim_tool        : ['fastp', 'trimmomatic'],
    assemblers       : ['spades', 'trinity', 'megahit'],
    spades_modes     : ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio'],
    cluster_method   : ['cdhitest', 'vsearch']
]

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Validate input parameters
WorkflowViralgenie.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config, params.adapter_fasta,
    params.host_index,params.host_genome,params.contaminants,
    params.spades_yml,params.spades_hmm
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input            ) { ch_input = file(params.input)                                      } else { exit 1, 'Input samplesheet not specified!'                              }
if (params.adapter_fasta    ) { ch_adapter_fasta = file(params.adapter_fasta)                      } else { ch_adapter_fasta  = []                                                  }
if (params.host_genome      ) { ch_host_genome = file(params.host_genome)                          } else { ch_host_genome = file(WorkflowMain.getGenomeAttribute(params, 'fasta')) }
if (params.host_index       ) { ch_host_index = Channel.fromPath(params.host_index).map{[[], it]}  } else { ch_host_index = []                                                      }
if (params.contaminants     ) { ch_contaminants = file(params.contaminants)                        } else { ch_contaminants = []                                                    }
if (params.spades_yml       ) { ch_spades_yml = file(params.spades_yml)                            } else { ch_spades_yml = []                                                      }
if (params.spades_hmm       ) { ch_spades_hmm = file(params.spades_hmm)                            } else { ch_spades_hmm = []                                                      }

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
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
include { PREPROCESSING_ILLUMINA          } from '../subworkflows/local/preprocessing_illumina'
include { FASTQ_KRAKEN_KAIJU              } from '../subworkflows/local/fastq_kraken_kaiju'
include { FASTQ_SPADES_TRINITY_MEGAHIT    } from '../subworkflows/local/fastq_spades_trinity_megahit'

//  Add consensus reconstruction of genome
include { FASTA_BLAST_CLUST               } from '../subworkflows/local/fasta_blast_clust'
include { ALIGN_COLLAPSE_CONTIGS          } from '../subworkflows/local/align_collapse_contigs'
include { UNPACK_DB as UNPACK_DB_BLAST    } from '../subworkflows/local/unpack_db'

// variant analysis
include { FASTQ_FASTA_ITERATIVE_CONSENSUS } from '../subworkflows/local/fastq_fasta_iterative_consensus'

// QC consensus

include { RENAME_FASTA_HEADER as RENAME_FASTA_HEADER_CONSTRAIN } from '../modules/local/rename_fasta_header'
include { UNPACK_DB as UNPACK_DB_CHECKV   } from '../subworkflows/local/unpack_db'
include { CONSENSUS_QC                    } from '../subworkflows/local/consensus_qc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BLAST_MAKEBLASTDB           } from '../modules/nf-core/blast/makeblastdb/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CHECKV_DOWNLOADDATABASE     } from '../modules/nf-core/checkv/downloaddatabase/main'


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
    ch_samplesheet = Channel.fromSamplesheet(
        'input'
        ).map{
            meta, read1, read2 ->
            [meta , [read1, read2]]
            }

    // preprocessing illumina reads
    PREPROCESSING_ILLUMINA (
        ch_samplesheet,
        ch_host_genome,
        ch_host_index,
        ch_adapter_fasta,
        ch_contaminants)
    ch_host_trim_reads      = PREPROCESSING_ILLUMINA.out.reads
    ch_decomplex_trim_reads = PREPROCESSING_ILLUMINA.out.reads_decomplexified
    ch_multiqc_files        = ch_multiqc_files.mix(PREPROCESSING_ILLUMINA.out.mqc.collect{it[1]}.ifEmpty([]))
    ch_versions             = ch_versions.mix(PREPROCESSING_ILLUMINA.out.versions)

    // Determining metagenomic diversity
    if (!params.skip_metagenomic_diversity) {
        FASTQ_KRAKEN_KAIJU(
            ch_host_trim_reads,
            params.kraken2_db,
            params.bracken_db,
            params.kaiju_db )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_KRAKEN_KAIJU.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_versions      = ch_versions.mix(FASTQ_KRAKEN_KAIJU.out.versions)
    }

    // Assembly
    ch_consensus = Channel.empty()
    if (!params.skip_assembly) {
        FASTQ_SPADES_TRINITY_MEGAHIT(
            ch_host_trim_reads,
            assemblers,
            ch_spades_yml,
            ch_spades_hmm)

        ch_versions      = ch_versions.mix(FASTQ_SPADES_TRINITY_MEGAHIT.out.versions)
        ch_contigs       = FASTQ_SPADES_TRINITY_MEGAHIT.out.scaffolds

        if (!params.skip_polishing){
            unpacked_references = UNPACK_DB_BLAST (params.reference_fasta).db
            ch_versions         = ch_versions.mix(UNPACK_DB_BLAST.out.versions)

            BLAST_MAKEBLASTDB ( unpacked_references )
            blast_clust_db      = BLAST_MAKEBLASTDB.out.db
            ch_versions         = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)


            // blast contigs against reference & identify clusters of (contigs & references)
            FASTA_BLAST_CLUST (
                ch_contigs,
                blast_clust_db,
                unpacked_references,
                params.cluster_method
                )
            ch_versions = ch_versions.mix(FASTA_BLAST_CLUST.out.versions)

            // Split up clusters into singletons and clusters of multiple contigs
            FASTA_BLAST_CLUST
                .out
                .centroids_members
                .branch{ meta, centroids, members ->
                    singletons: meta.cluster_size == 0
                    multiple: meta.cluster_size > 0
                }
                .set{ch_centroids_members}

            ch_centroids_members
                .singletons
                .map{ meta, centroids, members -> [ meta, centroids ] }
                .set{ ch_single_centroids }

            // Align clustered contigs & collapse into a single consensus per cluster
            ALIGN_COLLAPSE_CONTIGS (
                ch_centroids_members.multiple,
                params.contig_align_method
                )
            ch_versions = ch_versions.mix(ALIGN_COLLAPSE_CONTIGS.out.versions)


            SINGLETON_FILTERING (
                ch_single_centroids,
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
                    .set{ch_reference_reads}

            if (!params.skip_iterative_refinement) {
                FASTQ_FASTA_ITERATIVE_CONSENSUS (
                    ch_reference_reads,
                    params.iterative_repeats,
                    params.intermediate_mapper,
                    params.mapper,
                    params.with_umi,
                    params.deduplicate,
                    params.intermediate_variant_caller,
                    params.variant_caller,
                    params.intermediate_consensus_caller,
                    params.consensus_caller,
                    params.get_intermediate_stats,
                    params.get_stats,
                    params.min_contig_size,
                    params.max_n_1OOkbp
                )
            ch_consensus     = ch_consensus.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.consensus)
            ch_versions      = ch_versions.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTA_ITERATIVE_CONSENSUS.out.mqc.collect{it[1]}.ifEmpty([]))
            }
        }
    }

    // TODO: After consensus sequences have been made, we still have to map against it and call variants

    if (params.mapping_sequence ) {
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

        // TODO: map against constrain sequences

        // TODO: update consensus against constrain sequences

        //combine with all input samples & rename the meta.id
        // TODO: consider adding another column to the samplesheet with the mapping sequence
        ch_samplesheet
            .combine(ch_map_seq_anno)
            .map {meta, reads, seq, name ->
                sample = meta.id
                id = "${sample}_${name.id}"
                new_meta = meta + [id: id, sample: sample]
                return [new_meta, seq]
            }.set{ch_map_seq_anno_combined}

        //rename fasta headers
        RENAME_FASTA_HEADER_CONSTRAIN (ch_map_seq_anno_combined)
        ch_versions = ch_versions.mix(RENAME_FASTA_HEADER_CONSTRAIN.out.versions)

        //Add to the consensus channel, the mapping sequences will now always be mapped against
        ch_consensus = ch_consensus.mix(RENAME_FASTA_HEADER_CONSTRAIN.out.fasta)
        }

        if (!params.skip_consensus_qc) {
            ch_checkv_db = Channel.empty()
            if (!params.skip_checkv) {
                checkv_db = params.checkv_db
                if (!checkv_db) {
                    ch_checkv_db = CHECKV_DOWNLOADDATABASE().checkv_db
                } else {
                    ch_checkv_db = UNPACK_DB_CHECKV(checkv_db).db
                }
            }

            CONSENSUS_QC(
                ch_consensus,
                ch_checkv_db,
                params.skip_checkv,
                params.skip_quast
                )
            ch_versions      = ch_versions.mix(CONSENSUS_QC.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(CONSENSUS_QC.out.mqc.collect{it[1]}.ifEmpty([]))
        }

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



    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
