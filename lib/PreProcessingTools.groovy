//
// Functions that are needed within the subworkflow of preprocessing
//
import groovy.json.JsonSlurper

class PreProcessingTools {
    //
    // Function that parses fastp json output file to get total number of reads before trimming
    //
    public static Integer getFastpReadsBeforeFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['before_filtering']['total_reads'].toInteger()
    }

    public static Integer getFastpReadsAfterFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['after_filtering']['total_reads'].toInteger()
    }

    //
    // Function that parses Trimmomatic summary output file to get the total number of surviving reads
    //
    public static Integer getTrimmomaticReadsBeforeFiltering(file) {
        String fileContents = new File(file).text
        def pattern = /Input Read Pairs: (\d+\.?\d*)/
        def match = fileContents.find(pattern)
        return match.tokenize(': ').last().toInteger()
    }

    //
    // Function that parses Trimmomatic summary output file to get percentage of dropped reads after trimming
    //
    public static Integer getTrimmomaticReadsAfterFiltering(file) {
        String fileContents = new File(file).text
        def pattern = /Dropped Reads: (\d+\.?\d*)/
        def match = fileContents.find(pattern)
        return (getTrimmomaticReadsBeforeFiltering(file) - match.tokenize(': ').last().toInteger())
    }

    //
    // Functions that parses the tool functions in one to make it easier to understand
    //
    public static Integer getReadsBeforeStep(file, step){
        switch(tool) {
            case 'trimmomatic'  : return getTrimmomaticReadsBeforeFiltering(file)
            case 'fastp'        : return getFastpReadsBeforeFiltering(file)
            default             : break;
        }
    }

    public static Integer getReadsAfterStep(file, step){
        switch(tool) {
            case 'trimmomatic'  : return getTrimmomaticReadsAfterFiltering(file)
            case 'fastp'        : return getFastpReadsAfterFiltering(file)
            default             : break;
        }
    }
}
