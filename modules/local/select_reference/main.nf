process SELECT_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/16fd0599cbc5e52a5ac51f8668ed2c6988b4f44d461606e37953afcd581cd52d/data'
        : 'community.wave.seqera.io/library/biopython_pandas_python:671653bb7f9c4d5b'}"

    input:
    tuple val(meta), path(screen), path(reference), path(reads)

    output:
    tuple val(meta), path("*.json"), path("*_reference.fa"), path(reads), emit: fasta_reads
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_reference.py \\
        ${args} \\
        --mash ${screen} \\
        --reference ${reference} \\
        --prefix ${prefix}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    touch ${prefix}_reference.fa
    touch ${prefix}.json
    """
}
