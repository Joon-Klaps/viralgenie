process SSPACE_BASIC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sspace_basic:2.1.1--hdfd78af_0':
        'biocontainers/sspace_basic:2.1.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(contigs)
    tuple val(distance), val(deviation), val(complement)

    output:
    tuple val(meta), path("*.{fa,fasta}") , emit: fasta
    tuple val(meta), path("*.library.txt"), emit: library
    tuple val(meta), path("*.dot")        , optional:true, emit: dot
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads = reads.join('\t')
    def version = "2.1.1" // version not available through CLI of tool
    """
    echo -e "${prefix}\t${reads}\t${distance}\t${deviation}\t${complement}" > ${prefix}.library.txt
    sspace_basic \\
        -l ${prefix}.library.txt \\
        -s ${contigs} \\
        $args \\
        -T $task.cpus \\
        -b ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sspace_base: ${version}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads = reads.join('\t')
    def version = "2.1.1" // version not available through CLI of tool
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sspace_base: ${version}
    END_VERSIONS
    """
}
