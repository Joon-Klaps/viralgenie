process SSPACE_BASIC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sspace_basic:2.1.1--hdfd78af_1'
        : 'biocontainers/sspace_basic:2.1.1--hdfd78af_1'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(contigs)
    tuple val(distance), val(deviation), val(complement)
    val name

    output:
    tuple val(meta), path("*.final.renamed.scaffolds.fa") , emit: scaffolds
    tuple val(meta), path("*.final.scaffolds.fasta")      , emit: fasta
    tuple val(meta), path("*.library.txt")                , emit: library
    tuple val(meta), path("*.logfile.txt")                , emit: log
    tuple val(meta), path("*.summaryfile.txt")            , emit: summary
    tuple val(meta), path("dot/*.dot")                    , emit: dot, optional: true
    tuple val("${task.process}"), val('sspace_base'), val("2.1.1"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def name_command = name ?: 'sspace'
    def unzip_contig = "${contigs.getExtension()}" == "gz" ? "gunzip -c ${contigs}" : "cat ${contigs}"
    // doesn't allow insertion with <() or accepts gunzipped input

    """
    gunzip -f ${reads[0]}
    gunzip -f ${reads[1]}
    ${unzip_contig} > "tmp.fasta"

    echo "${prefix} ${reads[0].baseName} ${reads[1].baseName} ${distance} ${deviation} ${complement}" > ${prefix}.library.txt

    sspace_basic \\
        -l ${prefix}.library.txt \\
        -s tmp.fasta \\
        ${args} \\
        -T ${task.cpus} \\
        -b ${prefix}

    sed 's/>/>${name_command}_/g' ${prefix}.final.scaffolds.fasta > ${prefix}.final.renamed.scaffolds.fa
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_joined = reads.join('\t')
    """
    echo "${args} ${prefix} ${reads_joined} ${distance} ${deviation} ${complement}" > ${prefix}.library.txt
    touch ${prefix}.final.renamed.scaffolds.fa
    touch ${prefix}.final.scaffolds.fasta
    touch ${prefix}.library.txt
    touch ${prefix}.logfile.txt
    touch ${prefix}.summaryfile.txt
    """
}
