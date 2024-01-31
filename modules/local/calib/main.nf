process CALIB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://docker.io/jklaps/viralgenie:calib_pigz--d49838bb03867788' :
        'docker.io/jklaps/viralgenie:calib_pigz--ac030c8d900434e6' }"

    input:
    tuple val(meta), path(reads)
    val(tag_length)

    output:
    tuple val(meta), path("*.cluster"), emit: cluster
    tuple val(meta), path("*.")       , emit: reads
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: '--gzip-input'
    def args2     = task.ext.args2 ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def threads   = task.cpus <= 8  ? task.cpus : 8
    second_read   = meta.single_end ? '' : "-r ${reads[1]}"
    reads_joined  = meta.single_end ? reads : "${reads[0]} ${reads[1]}"
    output_prefix = meta.single_end ? "${prefix}_dedup." : "${prefix}_dedup1. ${prefix}_dedup2."

    """
    calib \\
        $args \\
        -f ${reads[0]} \\
        $second_read \\
        --threads $threads \\
        -l $tag_length \\
        -o ${prefix}.

    calib_cons \\
        $args2
        -c ${prefix}.cluster \\
        -q $reads_joined \\
        -o $output_prefix

    pigz -p $task.cpus *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calib: v.0.3.4
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args ?: '--gzip-input'
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def threads   = task.cpus <= 8  ? task.cpus : 8
    second_read   = meta.single_end ? '' : "-r ${reads[1]}"
    reads_joined  = meta.single_end ? reads : "${reads[0]} ${reads[1]}"
    output_prefix = meta.single_end ? "${prefix}_dedup." : "${prefix}_dedup1. ${prefix}_dedup2."
    """
    touch ${prefix}.cluster
    touch ${prefix}_dedup1.fastq.gz
    touch ${prefix}_dedup2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calib: v.0.3.4
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
