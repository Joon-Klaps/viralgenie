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
    // Function that parses and returns the number of mapped reasds from stats files
    //
    public static Integer getStatsMappedReads(statsFile) {
        def n_reads = 0
        statsFile.eachLine { line ->
            if (line =~ /SN\treads mapped:\s+(\d+)/) {
                n_reads = line.split('\t')[2].toInteger()
            }}
        return n_reads
    }

}
