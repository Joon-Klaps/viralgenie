//
// This file holds several functions specific to the main.nf workflow in the Joon-Klaps/viralgenie pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Define a global prefix
    //
    public static String getGlobalPrefix(workflow,params) {
        def date_stamp = new java.util.Date().format( 'yyyyMMdd')
        if (params.prefix) {
            return "${params.prefix}_${date_stamp}_${workflow.manifest.version}_${workflow.runName}"
        }
        return null
    }

    //
    // Print warning if genome fasta has more than one sequence
    //
    public static void isMultiFasta(fasta_file, log) {
        def count = 0
        def line  = null
        fasta_file.withReader { reader ->
            while (line = reader.readLine()) {
                if (line.contains('>')) {
                    count++
                    if (count > 1) {
                        log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            "  Multi-fasta genome files are not well supported by bowtie2 and bwamem2\n\n" +
                            "            Consider rerunning the pipeline with '--mapper bwa' \n" +
                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                        break
                    }
                }
            }
        }
    }

    //
    // Print warning if genome fasta has more than one sequence
    //
def calculateFastaLengthAndAmbiguousCount(File fastaFile) {
    def length = 0
    def ambiguousCount = 0

    fastaFile.eachLine { line ->
        if (line.startsWith(">")) {
            // Ignore header lines starting with ">"
            return
        } else {
            // Count the length of the sequence
            length += line.trim().length()
            // Count the occurrences of 'N' and ambiguous bases (e.g., 'R', 'Y', 'M', etc.)
            ambiguousCount += line.trim().findAll { it == 'N' }.size()
        }
    }

    return [length, ambiguousCount]
}


    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log, args) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that the profile doesn't contain spaces and doesn't end with a trailing comma
        checkProfile(workflow.profile, args, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            Nextflow.error("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
        }
    }

    //
    // Exit pipeline if --profile contains spaces
    //
    private static void checkProfile(profile, args, log) {
        if (profile.endsWith(',')) {
            Nextflow.error "Profile cannot end with a trailing comma. Please remove the comma from the end of the profile string.\nHint: A common mistake is to provide multiple values to `-profile` separated by spaces. Please use commas to separate profiles instead,e.g., `-profile docker,test`."
        }
        if (args[0]) {
            log.warn "nf-core pipelines do not accept positional arguments. The positional argument `${args[0]}` has been detected.\n      Hint: A common mistake is to provide multiple values to `-profile` separated by spaces. Please use commas to separate profiles instead,e.g., `-profile docker,test`."
        }
    }
}
