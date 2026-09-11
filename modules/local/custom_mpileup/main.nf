process CUSTOM_MPILEUP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysamstats:1.1.2--py39he47c912_12'
        : 'quay.io/biocontainers/pysamstats:1.1.2--py311h0152c62_12'}"

    input:
    tuple val(meta), path(bam), path(ref)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions
    tuple val("${task.process}"), val('numpy'), eval("python -c 'import numpy; print(numpy.__version__)'"), topic: versions
    tuple val("${task.process}"), val('pysamstats'), eval("python -c 'import pysamstats; print(pysamstats.__version__)'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    custom_mpileup.py \\
        ${args} \\
        --alignment ${bam} \\
        --reference ${ref} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}
