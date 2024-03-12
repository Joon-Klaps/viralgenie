process LOWCOV_TO_REFERENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0':
        'biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0' }"

    input:
    tuple val(meta), path(reference), path(consensus), path(mpileup)

    output:
    tuple val(meta), path("*.fasta") , emit: sequence
    tuple val(meta), path("*.txt")   , emit: txt, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lowcov_to_reference.py \\
        $args \\
        --reference ${reference} \\
        --consensus ${consensus} \\
        --mpileup ${mpileup} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        numpy: \$(pip show numpy | grep Version | sed 's/Version: //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        numpy: \$(pip show numpy | grep Version | sed 's/Version: //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
