process SSPACE_BASIC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sspace_basic:2.1.1--hdfd78af_1':
        'biocontainers/sspace_basic:2.1.1--hdfd78af_1' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(contigs)
    tuple val(distance), val(deviation), val(complement)
    val(name)

    output:
    tuple val(meta), path("${prefix}.final.renamed.scaffolds.fa")   , emit: scaffolds
    tuple val(meta), path("${prefix}.final.scaffolds.fasta")        , emit: fasta
    tuple val(meta), path("${prefix}.library.txt")                  , emit: library
    tuple val(meta), path("${prefix}.logfile.txt")                  , emit: log
    tuple val(meta), path("${prefix}.summaryfile.txt")              , emit: summary
    tuple val(meta), path("dot/*.dot")                              , optional:true, emit: dot
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    rename = name ? "sed 's/>/>${name}_/g' ${prefix}.final.scaffolds.fasta > ${prefix}.final.renamed.scaffolds.fasta" : ''
    def version = "2.1.1" // version not available through CLI of tool
    """
    gunzip -f ${reads[0]}
    gunzip -f ${reads[1]}

    echo "${prefix} ${reads[0].baseName} ${reads[1].baseName} ${distance} ${deviation} ${complement}" > ${prefix}.library.txt

    sspace_basic \\
        -l ${prefix}.library.txt \\
        -s ${contigs} \\
        $args \\
        -T $task.cpus \\
        -b ${prefix}

    ${rename}

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
    touch ${prefix}.final.scaffolds.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sspace_base: ${version}
    END_VERSIONS
    """
}
