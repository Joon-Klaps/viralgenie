// Helper workflows for getMapFromJson tests.
//
// PARSE      — wraps the production helper so the test can snapshot the
//              runtime class + Serializable status on our side of nf-test's
//              JSON conversion.
// See https://github.com/nf-core/viralmetagenome/issues/290.

include { getMapFromJson } from '../main.nf'

workflow PARSE {
    take:
    json_file
    lazy

    main:
    def parsed_map = !lazy ? getMapFromJson(json_file) : new groovy.json.JsonSlurper().parseText(json_file.text)

    emit:
    parsed          = Channel.value(parsed_map)
    class_name      = Channel.value(parsed_map.getClass().name)
    is_serializable = Channel.value(parsed_map instanceof java.io.Serializable)
}
