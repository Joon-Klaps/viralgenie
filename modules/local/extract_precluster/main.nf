process EXTRACT_PRECLUSTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78':
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(classifications), path(sequence)

    output:
    tuple val(meta), path("*fa")  , emit: sequences
    path "versions.yml"           , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_preclust \\
        ${classifications} \\
        ${sequence} \\
        ${prefix} \\
        $args \\

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
