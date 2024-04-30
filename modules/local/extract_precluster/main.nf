process EXTRACT_PRECLUSTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78':
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(kaiju_classifications)
    tuple val(meta2), path(kraken_classifications), path(kraken_report)
    tuple val(meta3), path(sequence)
    path(kaiju_db)

    output:
    tuple val(meta), path("*.fa"), path("*.json") , emit: sequences, optional: true
    tuple val(meta), path("*.resolved.txt")       , emit: resolved
    path "versions.yml"                           , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kaiju = kaiju_classifications ? "--kaiju-classifications  <(sort -k2,2 ${kaiju_classifications})" : ''
    def kaiju_db = kaiju_db ? "--database ${kaiju_db}" : ''
    def kraken = kraken_classifications ? "--kraken-classifications <(sort -k2,2  ${kraken_classifications})" : ''
    def kraken_report = kraken_report ? "--kraken-report ${kraken_report}" : ''
    """
    extract_preclust.py \\
        $args \\
        ${kaiju} \\
        ${kraken} \\
        ${kraken_report} \\
        ${kaiju_db} \\
        --sequences ${sequence} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_taxid29070.fa
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
