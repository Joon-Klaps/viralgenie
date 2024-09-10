process LOWCOV_TO_REFERENCE {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy { task.exitStatus == 4 ? 'ignore' : 'retry' } // can fail if mpileup empty

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/biopython_matplotlib:307f47953b7d175d':
        'community.wave.seqera.io/library/biopython_matplotlib:1991a6cc932d6beb' }"

    input:
    tuple val(meta), path(reference), path(consensus), path(mpileup)

    output:
    tuple val(meta), path("*.fa") , emit: sequence
    tuple val(meta), path("*.txt"), emit: txt, optional: true
    path "versions.yml"           , emit: versions

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
        matplotlib: \$(pip show matplotlib | grep Version | sed 's/Version: //g')
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
        matplotlib: \$(pip show matplotlib | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
