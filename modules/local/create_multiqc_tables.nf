process CREATE_MULTIQC_TABLES {
    label 'process_single'

    conda "bioconda:bioframe==0.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioframe:0.4.1--pyhdfd78af_0':
        'biocontainers/bioframe:0.4.1--pyhdfd78af_0' }"

    input:
    path clusters_summary_files
    path sample_metadata
    path checkv_files
    path blast_files

    output:
    path("summary_clusters_mqc.tsv")    , emit: summary_clusters_mqc, optional: true
    path("sample_metadata_mqc.tsv")     , emit: sample_metadata_mqc , optional: true
    path("summary_checkv_mqc.tsv")      , emit: summary_checkv_mqc  , optional: true
    path("summary_blast_mqc.tsv")       , emit: summary_blast_mqc   , optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def clusters_summary_files = clusters_summary_files ? "--clusters_summary ${clusters_summary_files.join(' ')}" : ''
    def sample_metadata        = sample_metadata        ? "--sample_metadata ${sample_metadata}"                   : ''
    def checkv_files           = checkv_files           ? "--checkv_files ${checkv_files.join(' ')}"               : ''
    def blast_files            = blast_files            ? "--blast_files ${blast_files.join(' ')}"                 : ''

    """
    create_multiqc_custom_tables.py\\
        $args \\
        $clusters_summary_files \\
        $sample_metadata \\
        $checkv_files \\
        $blast_files

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
    touch sample_metadata_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
