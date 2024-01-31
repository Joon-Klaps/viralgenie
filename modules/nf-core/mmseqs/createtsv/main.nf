
process MMSEQS_CREATETSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:15.6f452--pl5321h6a68c12_0':
        'biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_0' }"

    input:
    tuple val(meta), path(db_result), path(db_query)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args ?: "*.dbtype"
    def args3 = task.ext.args ?: "*.dbtype"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract files with specified args based suffix | remove suffix | isolate longest common substring of files
    DB_RESULT_PATH_NAME=\$(find -L "$db_result/" -maxdepth 1 -name "$args2" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )
    DB_QUERY_PATH_NAME=\$(find -L "$db_query/" -maxdepth 1 -name "$args3" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        createtsv \\
        \$DB_QUERY_PATH_NAME \\
        \$DB_QUERY_PATH_NAME \\
        \$DB_RESULT_PATH_NAME \\
        ${prefix}.tsv \\
        $args \\
        --threads ${task.cpus} \\
        --compressed 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """
}
