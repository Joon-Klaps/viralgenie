process BLAST_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0':
        'community.wave.seqera.io/library/pandas_pip_biopython:465ab8e47f7a0510' }"

    input:
    tuple val(meta), path(blast)
    tuple val(meta2), path(contigs)
    tuple val(meta3), path(db)

    output:
    tuple val(meta), path("*.hits.txt")  , emit: hits, optional: true
    tuple val(meta), path("*.fa")        , emit: sequence
    tuple val(meta), path("*.filter.tsv"), emit: filter, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def blast_command = blast ? "-i ${blast}" : ""
    """
    blast_filter.py \\
        $args \\
        ${blast_command} \\
        -c ${contigs} \\
        -r ${db} \\
        -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filter.tsv
    touch ${prefix}.filter.hits.txt
    touch ${prefix}_withref.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
