process KRONA_CLEANUP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04'}"

    input:
    tuple val(meta), path(krona, stageAs: 'uncleaned.krona.txt')

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('sed'), eval("sed --version | sed '1!d;s/.* //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Copy the file to a new name
    cp ${krona} ${prefix}.txt
    # Remove ugly 'x__' prefixes for each of the taxonomic levels
    LEVELS=(d k p c o f g s)
    for L in "\${LEVELS[@]}"; do
        sed -i "s/\${L}__//g" ${prefix}.txt
    done
    # Remove underscores that are standing in place of spaces
    sed -i "s/_/ /g" ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
