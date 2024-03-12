process CUSTOM_MULTIQC_TABLES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioframe:0.5.1--pyhdfd78af_0':
        'biocontainers/bioframe:0.5.1--pyhdfd78af_0' }"

    input:
    path clusters_summary_files
    path sample_metadata
    path checkv_files, stageAs: "?/*"
    path quast_files
    path blast_files
    path mapping_constrains
    path multiqc_dir
    path comment_headers
    path custom_table_headers

    output:
    path("summary_clusters_mqc.tsv")          , emit: summary_clusters_mqc  , optional: true
    path("sample_metadata_mqc.tsv")           , emit: sample_metadata_mqc   , optional: true
    path("contigs_overview_mqc.tsv")          , emit: contigs_overview_mqc  , optional: true
    path("summary_checkv_mqc.tsv")            , emit: summary_checkv_mqc    , optional: true
    path("summary_quast_mqc.tsv")             , emit: summary_quast_mqc     , optional: true
    path("summary_blast_mqc.tsv")             , emit: summary_blast_mqc     , optional: true
    path("mapping_constrains_mqc.tsv")        , emit: mapping_constrains_mqc, optional: true
    path("mapping_constrains_summary_mqc.tsv"), emit: constrains_summary_mqc, optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def clusters_summary_files = clusters_summary_files ? "--clusters_summary ${clusters_summary_files.join(' ')}" : ''
    def sample_metadata        = sample_metadata        ? "--sample_metadata ${sample_metadata}"                   : ''
    def checkv_files           = checkv_files           ? "--checkv_files ${checkv_files.join(' ')}"               : ''
    def quast_files            = quast_files            ? "--quast_files ${quast_files.join(' ')}"                 : ''
    def blast_files            = blast_files            ? "--blast_files ${blast_files.join(' ')}"                 : ''
    def mapping_constrains     = mapping_constrains     ? "--mapping_constrains ${mapping_constrains}"             : ''
    def multiqc_dir            = multiqc_dir            ? "--multiqc_dir ${multiqc_dir}"                           : ''
    def comment_headers        = comment_headers        ? "--comment_dir ${comment_headers}"                       : ''
    def custom_table_headers   = custom_table_headers   ? "--table_headers ${custom_table_headers}"                : ''

    """
    create_custom_multiqc_tables.py\\
        $args \\
        $clusters_summary_files \\
        $sample_metadata \\
        $checkv_files \\
        $quast_files \\
        $blast_files \\
        $mapping_constrains \\
        $comment_headers \\
        $custom_table_headers \\
        $multiqc_dir


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
        yaml: \$(pip show pyyaml | grep Version | sed 's/Version: //g')
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
        yaml: \$(pip show pyyaml | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
