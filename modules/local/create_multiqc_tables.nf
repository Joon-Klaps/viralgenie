process CREATE_MULTIQC_TABLES {
    label 'process_single'

    conda "bioconda:bioframe==0.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioframe:0.4.1--pyhdfd78af_0':
        'biocontainers/bioframe:0.4.1--pyhdfd78af_0' }"

    input:
    path clusters_summary_files

    output:
    path("summary_clusters_mqc.tsv")    , emit: summary_clusters_mqc, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def clusters_summary_files = clusters_summary_files ? "--clusters_summary ${clusters_summary_files.join(' ')}" : ''

    """
    create_multiqc_custom_tables.py\\
        $args \\
        $clusters_summary_files \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch cluster_summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
