process CLUSTER_EXTRACT {
    tag "meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta) , path(cluster)
    val(module)

    output:
    tuple val(meta), path('*_members.txt'), path('*_centroid.txt')   , emit: members_centroids
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def module = task.ext.args ?: "${module}"
    """
    extract_clust.py \\
        $module \\
        $cluster \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
