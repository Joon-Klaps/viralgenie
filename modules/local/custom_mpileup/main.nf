process CUSTOM_MPILEUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysamstats:1.1.2--py39he47c912_12':
        'quay.io/biocontainers/pysamstats:1.1.2--py311h0152c62_12' }"

    input:
    tuple val(meta), path(bam), path(ref)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    custom_mpileup.py \\
        $args \\
        --alignment ${bam} \\
        --reference ${ref} \\
        --prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version: | sed 's/Version: //g')
        pysam: \$(pip show pysam | grep Version: | sed 's/Version: //g')
        pysamstats: \$(pip show pysamstats | grep Version: | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version: | sed 's/Version: //g')
        pysam: \$(pip show pysam | grep Version: | sed 's/Version: //g')
        pysamstats: \$(pip show pysamstats | grep Version: | sed 's/Version: //g')
    """
}
