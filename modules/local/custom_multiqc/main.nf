process CUSTOM_MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d3/d36a0ca3cdc365286311ddd11caa877a0e9f6c39862fb700cb09da26d89da406/data'
        : 'community.wave.seqera.io/library/pip_multiqc_pandas:bffcc4521f62d9ac'}"

    input:
    path multiqc_files          , stageAs: "multiqc_files/?/*"
    path multiqc_config
    path clusters_summary_files
    path sample_metadata
    path checkv_files           , stageAs: "?/checkv/*"
    path quast_files            , stageAs: "?/quast/*"
    path blast_files            , stageAs: "?/blast/*"
    path mapping_constraints
    path anno_files             , stageAs: "?/annotation/*"
    path annotation_metadata
    path clusters_tsv           , stageAs: "?/clusters/*"
    path screen_files           , stageAs: "?/screen/*"
    path custom_table_headers

    output:

    path "contigs_overview.tsv"                , emit: contigs_table, optional: true
    path "contigs_overview-with-iterations.tsv", emit: contigs_it   , optional: true
    path "mapping_overview.tsv"                , emit: mapping_table, optional: true
    path "samples_overview.tsv"                , emit: samples_table, optional: true
    path "*multiqc_report.html"                , emit: report       , optional: true
    path "*_data"                              , emit: data         , optional: true
    path "*_plots"                             , emit: plots        , optional: true
    tuple val("${task.process}"), val('multiqc'), eval("multiqc --version | sed -e 's/multiqc, version //g'"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def multiqc_files_command = multiqc_files ? "--multiqc_files multiqc_files" : ''
    def multiqc_config_command = multiqc_config ? "--multiqc_config ${multiqc_config}" : ''
    def clusters_summary_files_command = clusters_summary_files ? "--clusters_summary ${clusters_summary_files}" : ''
    def sample_metadata_command = sample_metadata ? "--sample_metadata ${sample_metadata}" : ''
    def checkv_files_command = checkv_files ? "--checkv_files ${checkv_files}" : ''
    def quast_files_command = quast_files ? "--quast_files ${quast_files}" : ''
    def blast_files_command = blast_files ? "--blast_files ${blast_files}" : ''
    def annotation_files = anno_files ? "--annotation_files ${anno_files}" : ''
    def annotation_metadata_command = annotation_metadata ? "--annotation_metadata ${annotation_metadata}" : ''
    def clusters_files = clusters_tsv ? "--clusters_files ${clusters_tsv}" : ''
    def mapping_constraints_command = mapping_constraints ? "--mapping_constraints ${mapping_constraints}" : ''
    def screen_files_command = screen_files ? "--screen_files ${screen_files}" : ''
    def custom_table_headers_command = custom_table_headers ? "--table_headers ${custom_table_headers}" : ''

    """
    custom_multiqc.py \\
        ${args} \\
        ${multiqc_files_command} \\
        ${multiqc_config_command} \\
        ${clusters_summary_files_command} \\
        ${sample_metadata_command} \\
        ${checkv_files_command} \\
        ${quast_files_command} \\
        ${blast_files_command} \\
        ${annotation_files} \\
        ${annotation_metadata_command} \\
        ${clusters_files} \\
        ${mapping_constraints_command} \\
        ${screen_files_command} \\
        ${custom_table_headers_command}
    """

    stub:
    """
    touch contigs_overview.tsv
    touch contigs_overview-with-iterations.tsv
    touch mapping_overview.tsv
    touch samples_overview.tsv
    touch multiqc_report.html
    mkdir -p data
    mkdir -p plots
    """
}
