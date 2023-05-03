//
// This file holds several functions specific to the workflow/viralgenie.nf in the nf-core/viralgenie pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowViralgenie {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {

        if (!valid_params['trim_tool'].contains(params.trim_tool)) {
            Nextflow.error("Please specify a valid trimming tool: 'fastp' or 'trimmomatic' not ${params.trim_tool}.")
        }
        //check if all values of assembler are in valid_params.assemblers
        if (params.assemblers) {
            for (assembler in params.assemblers.split(',').collect{ it.trim().toLowerCase() }) {
                if (!(assembler in valid_params['assemblers'])) {
                    Nextflow.error("${assembler} is not a valid assembler. Please choose from ${valid_params['assemblers'].join(', ')}")
                }
            }
        }
        if (!valid_params['spades_modes'].contains(params.spades_mode)) {
            Nextflow.error("${params.spades_modes} is not a valid spades mode. Please choose from ${valid_params['spades_mode'].join(', ')}")
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }
}
