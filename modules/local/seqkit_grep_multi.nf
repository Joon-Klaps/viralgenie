process SEQKIT_GREP_MULTI {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.5.1--h9ee0642_0':
        'biocontainers/seqkit:2.5.1--h9ee0642_0' }"

    // Modified from original nf-core file
    input:
    tuple val(meta), path(txt_pattern_files), path (sequence)

    output:
    tuple val(meta), path("*.{fa,fq}")  , emit: filter
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def suffix = task.ext.suffix ?: "${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"

    """
    for file in ${txt_pattern_files.join(' ')};
    do
        prefix=\$( basename \${file} | sed 's/\\.[^.]*\$//' )
        echo "Processing file: \${prefix} .."
        seqkit \\
            grep \\
            $args \\
            --threads $task.cpus \\
            -f \${file} \\
            ${sequence} \\
            -o \$prefix.${suffix} \\
            && echo "Done" || { echo "Failed"; exit 1; }
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // fasta or fastq. Exact pattern match .fasta or .fa suffix with optional .gz (gzip) suffix
    def suffix = task.ext.suffix ?: "${sequence}" ==~ /(.*f[astn]*a(.gz)?$)/ ? "fa" : "fq"

    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
