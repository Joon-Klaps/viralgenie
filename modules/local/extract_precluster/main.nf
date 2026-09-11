process EXTRACT_PRECLUSTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78':
        'biocontainers/biopython:1.78' }"

    input:
    tuple val(meta), path(kaiju_classifications)
    tuple val(meta2), path(kraken_classifications), path(kraken_report)
    tuple val(meta3), path(sequence)
    path(kaiju_db)

    output:
    tuple val(meta), path("*.{fa,fasta}"), path("*.json") , emit: sequences, optional: true
    tuple val(meta), path("*.resolved.txt")       , emit: resolved
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions
    tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kaiju = kaiju_classifications ? "--kaiju-classifications  <(sort -k2,2 ${kaiju_classifications})" : ''
    def kaiju_db_argument = kaiju_db ? "--database ${kaiju_db}" : ''
    def kraken = kraken_classifications ? "--kraken-classifications <(sort -k2,2  ${kraken_classifications})" : ''
    def kraken_report_argument = kraken_report ? "--kraken-report ${kraken_report}" : ''

    """
    extract_preclust.py \\
        $args \\
        ${kaiju} \\
        ${kraken} \\
        ${kraken_report_argument} \\
        ${kaiju_db_argument} \\
        --sequences ${sequence} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_taxid0000.fasta
    touch ${prefix}.resolved.txt
    cat <<-END_JSON > ${prefix}.json
    {
        "ntaxa": 1
    }
    END_JSON
    """
}
