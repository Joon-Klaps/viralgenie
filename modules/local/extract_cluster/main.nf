process EXTRACT_CLUSTER {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78':
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(clusters), path(seq)
    val(module)

    output:
    tuple val(meta), path('*_members.fa'), path('*_centroid.fa')  , emit: members_centroids
    tuple val(meta), path("*.json")                               , emit: json
    tuple val(meta), path("*.clusters.tsv")                       , emit: tsv
    tuple val(meta), path("*.summary_mqc.tsv")                    , emit: summary
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    extract_clust.py \\
        $args \\
        -m $module \\
        -c ${clusters.join(' ')} \\
        -s $seq \\
        -p $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    touch ${prefix}_cl0_members.fa
    touch ${prefix}_cl0_centroid.fa
    touch ${prefix}_summary_mqc.tsv
    touch ${prefix}_clusters.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        biopython: \$(pip show biopython | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

}
