process PRINSEQPLUSPLUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prinseq-plus-plus:1.2.3--hc90279e_1':
        'quay.io/biocontainers/prinseq-plus-plus:1.2.3--hc90279e_1' }"

    input:
    tuple val(meta), path(reads), path(fasta)

    output:
    tuple val(meta), path("*_good_out*")                  , emit: good_reads
    tuple val(meta), path("*_single_out*"), optional: true, emit: single_reads
    tuple val(meta), path("*_bad_out*")   , optional: true, emit: bad_reads
    tuple val(meta), path("*.log")                        , emit: log
    tuple val("${task.process}"), val('prinseqplusplus'), eval("prinseq++ --version | cut -f 2 -d ' '"), topic: versions, emit: versions_prinseqplusplus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastqs = reads ? meta.single_end ? "-fastq ${reads}" : "-fastq ${reads[0]} -fastq2 ${reads[1]}" : ''
    def fasta_arg = fasta ? "-fastq ${fasta} -FASTA" : ''

    """
    prinseq++ \\
        -threads $task.cpus \\
        ${fastqs} \\
        ${fasta_arg} \\
        -out_name ${prefix} \\
        -VERBOSE 1 \\
        $args \\
        | tee ${prefix}.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_good_out.fastq.gz
    touch ${prefix}.log
    """
}
