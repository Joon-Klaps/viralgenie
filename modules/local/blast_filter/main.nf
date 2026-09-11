process BLAST_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/16/16fd0599cbc5e52a5ac51f8668ed2c6988b4f44d461606e37953afcd581cd52d/data'
        : 'community.wave.seqera.io/library/biopython_pandas_python:671653bb7f9c4d5b'}"

    input:
    tuple val(meta), path(blast)
    tuple val(meta2), path(contigs)
    path(blacklist)
    tuple val(meta3), path(db)

    output:
    tuple val(meta), path("*.hits.txt"), emit: hits, optional: true
    tuple val(meta), path("*.fa"), emit: sequence
    tuple val(meta), path("*.filter.tsv"), emit: filter, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def blast_command = blast ? "-i ${blast}" : ""
    def blacklist_arg = blacklist ? "-k ${blacklist}" : ""
    """
    blast_filter.py \\
        ${args} \\
        ${blast_command} \\
        ${blacklist_arg} \\
        -c ${contigs} \\
        -r ${db} \\
        -p ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filter.tsv
    touch ${prefix}.filter.hits.txt
    touch ${prefix}_withref.fa
    """
}
