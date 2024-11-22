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

    public static Map getLengthAndAmbigous(fastaFile) {
        def length = 0
        def ambiguousCount = 0
        def ambiguousPerc = 0

        fastaFile.eachLine { line ->
            if (line.startsWith(">")) {
                // Ignore header lines starting with ">"
                return
            } else {
                // Count the length of the sequence
                length += line.trim().length()
                // Count the occurrences of 'N'
                ambiguousCount += line.count("N")
            }
        }
        if (length.toInteger() > 0) {
            ambiguousPerc = (ambiguousCount / length) * 100
        }

        return [contig_size: length.toInteger(), n_100 :ambiguousPerc.toInteger()]
    }

    //
    // Function that parses and returns the number of mapped reasds from flagstat files
    //
    public static Integer getFlagstatMappedReads(flagstat_file) {
        def mapped_reads = 0
        flagstat_file.eachLine { line ->
            if (line.contains(' mapped (')) {
                mapped_reads = line.tokenize().first().toInteger()
            }
        }
        return mapped_reads
}

}
