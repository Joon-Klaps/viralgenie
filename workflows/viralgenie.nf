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
if (params.input         ) { ch_input = file(params.input)                                      } else { exit 1, 'Input samplesheet not specified!' }
if (params.adapter_fasta ) { ch_adapter_fasta = file(params.adapter_fasta)                      } else { ch_adapter_fasta  = []                     }
if (params.host_genome   ) { ch_host_genome = file(params.host_genome)                          } else { ch_host_genome = []                        }
if (params.host_index    ) { ch_host_index = Channel.fromPath(params.host_index).map{[[], it]}  } else { ch_host_index = []                         }
if (params.contaminants  ) { ch_contaminants = file(params.contaminants)                        } else { ch_contaminants = []                       }
if (params.spades_yml    ) { ch_spades_yml = file(params.spades_yml)                            } else { ch_spades_yml = []                         }
if (params.spades_hmm    ) { ch_spades_hmm = file(params.spades_hmm)                            } else { ch_spades_hmm = []                         }

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
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                   } from '../subworkflows/local/input_check'
include { PREPROCESSING_ILLUMINA        } from '../subworkflows/local/preprocessing_illumina'
include { FASTQ_KRAKEN_KAIJU            } from '../subworkflows/local/fastq_kraken_kaiju'
include { FASTQ_SPADES_TRINITY_MEGAHIT  } from '../subworkflows/local/fastq_spades_trinity_megahit'
//  Add consensus reconstruction of genome
include { FASTA_BLAST_CLUST             } from '../subworkflows/local/fasta_blast_clust'
include { ALIGN_COLLAPSE_CONTIGS        } from '../subworkflows/local/align_collapse_contigs'
include { CONSENSUS_QC                  } from '../subworkflows/local/consensus_qc'
include { UNPACK_DB as UNPACK_DB_BLAST  } from '../subworkflows/local/unpack_db'
include { UNPACK_DB as UNPACK_DB_CHECKV } from '../subworkflows/local/unpack_db'

// TODO: Add identification intrahost variability

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

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    // Taken from viralrecon
    //
    INPUT_CHECK(ch_input)

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // preprocessing illumina reads
    PREPROCESSING_ILLUMINA (
        INPUT_CHECK.out.reads,
        ch_host_genome,
        ch_host_index,
        ch_adapter_fasta,
        ch_contaminants)
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING_ILLUMINA.out.mqc.collect{it[1]}.ifEmpty([]))
    ch_versions      = ch_versions.mix(PREPROCESSING_ILLUMINA.out.versions)

    // Determining metagenomic diversity
    if (!params.skip_metagenomic_diversity) {
        FASTQ_KRAKEN_KAIJU(
            PREPROCESSING_ILLUMINA.out.reads,
            params.kraken2_db,
            params.bracken_db,
            params.kaiju_db )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_KRAKEN_KAIJU.out.mqc.collect{it[1]}.ifEmpty([]))
        ch_versions      = ch_versions.mix(FASTQ_KRAKEN_KAIJU.out.versions)
    }

    if (!params.skip_assembly) {
        FASTQ_SPADES_TRINITY_MEGAHIT(
            PREPROCESSING_ILLUMINA.out.reads,
            assemblers,
            ch_spades_yml,
            ch_spades_hmm)
        ch_versions      = ch_versions.mix(FASTQ_SPADES_TRINITY_MEGAHIT.out.versions)

        if (!params.skip_polishing){
            ch_db_blast = UNPACK_DB_BLAST (params.reference_fasta).db
            ch_versions = ch_versions.mix(UNPACK_DB_BLAST.out.versions)
            BLAST_MAKEBLASTDB ( ch_db_blast )
            ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

            // blast contigs against reference & identify clustering contigs & references
            FASTA_BLAST_CLUST (
                FASTQ_SPADES_TRINITY_MEGAHIT.out.scaffolds,
                BLAST_MAKEBLASTDB.out.db,
                ch_db_blast,
                params.cluster_method
                )

            //TODO: branch those with a single or no members
            ch_members = FASTA_BLAST_CLUST.out.members
            ch_centroids  = FASTA_BLAST_CLUST.out.centroids
            ALIGN_COLLAPSE_CONTIGS(
                ch_members,
                ch_centroids
                params.align_method
                )

            ch_consensus = ALIGN_COLLAPSE_CONTIGS.out.consensus

            //TODO: subworkflow for iterative refinement, contains another subworkflow if we just give a single reference

            }
        } else {
            ch_reference = params.mapping_sequence

        }

        if (!params.skip_consensus_qc) {
            if (!skip_checkv) {
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
                params.skip_quast,
                )
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

    ch_multiqc_files = Channel.empty()
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
