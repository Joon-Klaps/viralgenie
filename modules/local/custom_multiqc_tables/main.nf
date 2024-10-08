process CUSTOM_MULTIQC_TABLES {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_dash_multiqc_pandas:9472df82fedd8448':
        'community.wave.seqera.io/library/pip_dash_multiqc_pandas:1046792731cd7619' }"

    input:
    path clusters_summary_files
    path sample_metadata
    path checkv_files, stageAs: "?/checkv/*"
    path quast_files, stageAs: "?/quast/*"
    path blast_files, stageAs: "?/blast/*"
    path bed_files, stageAs: "?/bed/*"
    path mapping_constrains
    path anno_files, stageAs: "?/annotation/*"
    path clusters_tsv, stageAs: "?/clusters/*"
    path screen_files, stageAs: "?/screen/*"
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
    path("summary_anno_mqc.tsv")              , emit: summary_anno_mqc      , optional: true
    path("clusters_barchart.tsv")             , emit: clusters_barchart_mqc , optional: true
    path("contig_custom_table_mqc.html")      , emit: contig_html           , optional: true
    path("constrain_custom_table_mqc.html")   , emit: mapping_constrains_mqc, optional: true
    path("mapping_constrains_mqc.tsv")        , emit: mapping_constrains_mqc , optional: true
    path("mapping_constrains_summary_mqc.tsv"), emit: constrains_summary_mqc , optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def clusters_summary_files = clusters_summary_files ? "--clusters_summary ${clusters_summary_files}" : ''
    def sample_metadata        = sample_metadata        ? "--sample_metadata ${sample_metadata}"         : ''
    def checkv_files           = checkv_files           ? "--checkv_files ${checkv_files}"               : ''
    def quast_files            = quast_files            ? "--quast_files ${quast_files}"                 : ''
    def blast_files            = blast_files            ? "--blast_files ${blast_files}"                 : ''
    def annotation_files       = anno_files             ? "--annotation_files ${anno_files}"             : ''
    def bed_files              = bed_files              ? "--bed_files ${bed_files}"                     : ''
    def clusters_files         = clusters_tsv           ? "--clusters_files ${clusters_tsv}"             : ''
    def mapping_constrains     = mapping_constrains     ? "--mapping_constrains ${mapping_constrains}"   : ''
    def screen_files           = screen_files           ? "--screen_files ${screen_files}"               : ''
    def multiqc_dir            = multiqc_dir            ? "--multiqc_dir ${multiqc_dir}"                 : ''
    def comment_headers        = comment_headers        ? "--comment_dir ${comment_headers}"             : ''
    def custom_table_headers   = custom_table_headers   ? "--table_headers ${custom_table_headers}"      : ''

    """
    custom_multiqc_tables.py \\
        $args \\
        $clusters_summary_files \\
        $sample_metadata \\
        $checkv_files \\
        $quast_files \\
        $blast_files \\
        $bed_files \\
        $clusters_files \\
        $mapping_constrains \\
        $annotation_files \\
        $screen_files\\
        $comment_headers \\
        $custom_table_headers \\
        $multiqc_dir


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version | sed 's/Version: //g')
        yaml: \$(pip show pyyaml | grep Version | sed 's/Version: //g')
        plotly: \$(pip show plotly | grep Version | sed 's/Version: //g')
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
        plotly: \$(pip show plotly | grep Version | sed 's/Version: //g')
    END_VERSIONS
    """
}
