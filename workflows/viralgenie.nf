/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    trim_tool        : ['fastp', 'trimmomatic'],
    assemblers       : ['spades', 'trinity', 'megahit'],
    spades_modes     : ['rnaviral', 'corona', 'metaviral', 'meta', 'metaplasmid', 'plasmid', 'isolate', 'rna', 'bio'],
    cluster_method   : ['cdhit', 'vsearch']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowViralgenie.initialise(params, log,valid_params)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config, params.adapter_fasta,
    params.host_index,params.host_reference,params.contaminants,
    params.spades_yml,params.spades_hmm
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input         ) { ch_input = file(params.input)                                      } else { exit 1, 'Input samplesheet not specified!' }
if (params.adapter_fasta ) { ch_adapter_fasta = file(params.adapter_fasta)                      } else { ch_adapter_fasta  = []                     }
if (params.host_reference) { ch_host_reference = file(params.host_reference)                    } else { ch_host_reference = []                     }
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
include { INPUT_CHECK                  } from '../subworkflows/local/input_check'
include { PREPROCESSING_ILLUMINA       } from '../subworkflows/local/preprocessing_illumina'
include { FASTQ_KRAKEN_KAIJU           } from '../subworkflows/local/fastq_kraken_kaiju'
include { FASTQ_SPADES_TRINITY_MEGAHIT } from '../subworkflows/local/fastq_spades_trinity_megahit'
//  Add consensus reconstruction of genome
//include { FASTA_FASTQ_BOWTIE2_METABAT2 } from '../subworkflows/local/fasta_fastq_bowtie2_metabat2'
include { FASTA_BLAST_CLUST            } from '../subworkflows/local/fasta_blast_clust'
// TODO: Add identification intrahost variability

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { UNTAR as UNTAR_BLAST_DB     } from '../modules/nf-core/untar/main'
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
    INPUT_CHECK (
        ch_input,
    )

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // preprocessing illumina reads
    PREPROCESSING_ILLUMINA (
        INPUT_CHECK.out.reads,
        ch_host_reference,
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

        //TODO: Binning of assemblies
        if (!params.skip_polishing){
            // FASTA_FASTQ_BOWTIE2_METABAT2(
            //     FASTQ_SPADES_TRINITY_MEGAHIT.out.scaffolds,
            //     PREPROCESSING_ILLUMINA.out.reads
            // )
            db_blast= params.reference_fasta
            if (db_blast.endsWith('.gz')) {
                UNTAR_BLAST_DB (
                    [ [:], db_blast ]
                )
                ch_db_blast = UNTAR_BLAST_DB.out.untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_BLAST_DB.out.versions)
            } else {
                        ch_db_blast = Channel.value(file(db_blast))
            }
            BLAST_MAKEBLASTDB ( ch_db_blast )
            ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions.first())

            FASTA_BLAST_CLUST (
                FASTQ_SPADES_TRINITY_MEGAHIT.out.scaffolds,
                BLAST_MAKEBLASTDB.out.db,
                ch_db_blast,
                params.cluster_method
                )

            //TODO: Filter bins further down if necessary

            //TODO: reference Identification


            //TODO: Scaffolding & consensus reconstruction of genome

            }
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

    methods_description    = WorkflowViralgenie.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
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
