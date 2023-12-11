process KAIJU_MERGEOUTPUTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.8.2--h5b5514e_1':
        'biocontainers/kaiju:1.8.2--h5b5514e_1' }"

    input:
    tuple val(meta), path(kaiju), path(kraken_classified), path(kraken_unclassified)
    path (db)

    output:
    tuple val(meta), path("*.txt"), emit: merged
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p tmp && \\
        cat ${kraken_classified} ${kraken_unclassified} > tmp/kraken.txt

    dbnodes=`find -L ${db} -name "*nodes.dmp"`

    kaiju-mergeOutputs \\
        -i <(sort -k2,2 ${kaiju}) \\
        -j <(sort -k2,2 tmp/kraken.txt) \\
        -o ${prefix}.tsv \\
        -t \$dbnodes \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
