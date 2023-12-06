process NETWORK_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://docker.io/jklaps/viralgenie:igraph_leidenalg_matplotlib_pycairo_pandas--c1a94e30d4ecf531':
        'docker.io/jklaps/viralgenie:igraph_leidenalg_matplotlib_pycairo_pandas--413697ff28400e7c' }"

    input:
    tuple val(meta), path(dist)
    val(cluster_method)

    output:
    tuple val(meta), path("*.tsv"), emit: clusters
    tuple val(meta), path("*.png"), emit: png, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    network_cluster.py \\
        $args \\
        --method $cluster_method \\
        --prefix $prefix \\
        $dist \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
        matplotlib: \$(pip show matplotlib | grep Version | sed 's/Version: //g')
        igraph: \$(pip show igraph | grep Version | sed 's/Version: //g')
        leidenalg: \$(pip show leidenalg | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
        matplotlib: \$(pip show matplotlib | grep Version | sed 's/Version: //g')
        igraph: \$(pip show igraph | grep Version | sed 's/Version: //g')
        leidenalg: \$(pip show leidenalg | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
