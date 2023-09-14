process BLAST_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda:bioframe==0.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioframe:0.4.1--pyhdfd78af_0':
        'biocontainers/bioframe:0.4.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(blast)

    output:
    tuple val(meta), path("*.hits.txt")  , emit: hits
    tuple val(meta), path("*.filter.tsv"), emit: filter
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    blast_filter.py \\
        $blast \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filter.txt
    touch ${prefix}.hits.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
