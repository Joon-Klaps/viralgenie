process LOWCOV_TO_REFERENCE {
    tag "$meta.id"
    label 'process_single'
    errorStrategy { task.exitStatus == 4 ? 'ignore' : 'retry' } // can fail if mpileup empty

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://docker.io/jklaps/viralgenie:4be7b5d404272a48':
        'docker.io/jklaps/viralgenie:4be7b5d404272a48' }"

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
        pymuscle5: \$
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
