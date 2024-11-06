process CUSTOM_MULTIQC_TABLES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73a8b9ad157a08e5aa8a83cdb1975f9d4ff46f112351ebb9427b27dc32763c76/data':
        'community.wave.seqera.io/library/pip_multiqc_pandas:d84d38acb8d47ed1' }"

    input:
    path multiqc_files, stageAs: "multiqc_files/?/*"
    path multiqc_config
    path clusters_summary_files
    path sample_metadata
    path checkv_files, stageAs: "?/checkv/*"
    path quast_files, stageAs: "?/quast/*"
    path blast_files, stageAs: "?/blast/*"
    path mapping_constrains
    path anno_files, stageAs: "?/annotation/*"
    path clusters_tsv, stageAs: "?/clusters/*"
    path screen_files, stageAs: "?/screen/*"
    path comment_headers
    path custom_table_headers

    output:

    path("summary_clusters_mqc.tsv")          , emit: summary_clusters_mqc  , optional: true
    path("sample_metadata_mqc.tsv")           , emit: sample_metadata_mqc   , optional: true
    path("contigs_overview_mqc.tsv")          , emit: contigs_overview_mqc  , optional: true
    path("summary_checkv_mqc.tsv")            , emit: summary_checkv_mqc    , optional: true
    path("summary_quast_mqc.tsv")             , emit: summary_quast_mqc     , optional: true
    path("summary_blast_mqc.tsv")             , emit: summary_blast_mqc     , optional: true
    path("summary_anno_mqc.tsv")              , emit: summary_anno_mqc      , optional: true
    path("clusters_barchart.tsv")             , emit: clusters_barchart_mqc , optional: true
    path("mapping_constrains_summary_mqc.tsv"), emit: constrains_summary_mqc, optional: true
    path "*multiqc_report.html"               , emit: report                , optional: true
    path "*_data"                             , emit: data                  , optional: true
    path "*_plots"                            , emit: plots                 , optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def multiqc_files_command          = multiqc_files          ? "--multiqc_files multiqc_files"                : '' // Just refer to the dir for now.
    def multiqc_config_command         = multiqc_config         ? "--multiqc_config ${multiqc_config}"           : ''
    def clusters_summary_files_command = clusters_summary_files ? "--clusters_summary ${clusters_summary_files}" : ''
    def sample_metadata_command        = sample_metadata        ? "--sample_metadata ${sample_metadata}"         : ''
    def checkv_files_command           = checkv_files           ? "--checkv_files ${checkv_files}"               : ''
    def quast_files_command            = quast_files            ? "--quast_files ${quast_files}"                 : ''
    def blast_files_command            = blast_files            ? "--blast_files ${blast_files}"                 : ''
    def annotation_files               = anno_files             ? "--annotation_files ${anno_files}"             : ''
    def clusters_files                 = clusters_tsv           ? "--clusters_files ${clusters_tsv}"             : ''
    def mapping_constrains_command     = mapping_constrains     ? "--mapping_constrains ${mapping_constrains}"   : ''
    def screen_files_command           = screen_files           ? "--screen_files ${screen_files}"               : ''
    def comment_headers_command        = comment_headers        ? "--comment_dir ${comment_headers}"             : ''
    def custom_table_headers_command   = custom_table_headers   ? "--table_headers ${custom_table_headers}"      : ''

    """
    custom_multiqc_tables.py \\
        $args \\
        $multiqc_files_command \\
        $multiqc_config_command \\
        $clusters_summary_files_command \\
        $sample_metadata_command \\
        $checkv_files_command \\
        $quast_files_command \\
        $blast_files_command \\
        $annotation_files \\
        $clusters_files \\
        $mapping_constrains_command \\
        $screen_files_command \\
        $comment_headers_command \\
        $custom_table_headers_command \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version: | sed 's/Version: //g')
        yaml: \$(pip show pyyaml | grep Version: | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch cluster_summary_mqc.tsv
    touch sample_metadata_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version: | sed 's/Version: //g')
        yaml: \$(pip show pyyaml | grep Version: | sed 's/Version: //g')
    END_VERSIONS
    """
}
