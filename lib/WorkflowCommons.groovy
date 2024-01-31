//
// This file holds several functions common to the multiple workflows in the nf-core/viralrecon pipeline
//
import groovy.json.JsonSlurper

class WorkflowCommons {

    //
    // Create MultiQC tsv custom content from a list of values
    //
    public static String multiqcTsvFromList(tsv_data, header, comments) {
        def tsv_string = ""
        if (tsv_data.size() > 0) {
            if (comments) tsv_string += "# ${comments.join('\n# ')}\n"
            tsv_string += "${header.join('\t')}\n"
            tsv_string += tsv_data.join('\n')
        }
        return tsv_string
    }

    //
    // Function to get column entries from a file
    //
    public static ArrayList getColFromFile(input_file, col=0, uniqify=false, sep='\t') {
        def vals = []
        input_file.eachLine { line ->
            def val = line.split(sep)[col]
            if (uniqify) {
                if (!vals.contains(val)) {
                    vals << val
                }
            } else {
                vals << val
            }
        }
        return vals
    }

    //
    // Function that returns the number of lines in a file
    //
    public static Integer getNumLinesInFile(input_file) {
        def num_lines = 0
        input_file.eachLine { line ->
            num_lines ++
        }
        return num_lines
    }

    //
    // Function to get number of variants reported in BCFTools stats file
    //
    public static Integer getNumVariantsFromBCFToolsStats(bcftools_stats) {
        def num_vars = 0
        bcftools_stats.eachLine { line ->
            def matcher = line =~ /SN\s*0\s*number\sof\srecords:\s*([\d]+)/
            if (matcher) num_vars = matcher[0][1].toInteger()
        }
        return num_vars
    }

    //
    // Function to get a Map from a JSON file
    //
    public static Map getMapFromJson(json_file) {
        def Map json = (Map) new JsonSlurper().parse(json_file)
        return json
    }
}
