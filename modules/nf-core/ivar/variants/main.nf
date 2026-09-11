process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ivar:1.4.4--h077b44d_0' :
        'quay.io/biocontainers/ivar:1.4.4--h077b44d_0' }"

    input:
    tuple val(meta), path(bam)
    path  fasta
    path  fai
    path  gff
    val   save_mpileup

    output:
    tuple val(meta), path("*.tsv")    , emit: tsv
    tuple val(meta), path("*.mpileup"), optional:true, emit: mpileup
    tuple val("${task.process}"), val('ivar'), eval("ivar version | sed -n 's|iVar version \\(.*\\)|\\1|p'"), emit: versions_ivar, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/^samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def features = gff ? "-g $gff" : ""

    // `samtools mpileup` intermittently fails on network-backed storage, so retry a few
    // times before giving up instead of failing the whole run on a transient read error.
    def max_retries = 3
    """
    set -e  # Exit on error
    retry_count=0
    while [ \$retry_count -lt $max_retries ]; do
        samtools_log=\$(samtools mpileup --reference $fasta $args2 $bam --output ${prefix}.mpileup 2>&1)
        error_message=\$(echo "\$samtools_log" | grep 'E::' | head -n1 || true)
        if [ -n "\$error_message" ]; then
            echo "Samtools error: \$error_message"
            echo "Retrying (\$((\$retry_count + 1))/$max_retries)..."
            retry_count=\$((\$retry_count + 1))
        else
            break
        fi

        sleep 10
    done

    # Check if the maximum number of retries is reached
    if [ \$retry_count -eq $max_retries ]; then
            echo ""
            echo "Unable to solve: \$error_message"
            echo "Maximum number of retries reached. Exiting."
        exit 1
    else
        echo "Samtools mpileup successful"
        cat ${prefix}.mpileup | ivar variants \\
            -p $prefix \\
            -r $fasta \\
            $args \\
            $features
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_mpileup = save_mpileup ? "touch ${prefix}.mpileup" : ''
    """
    touch ${prefix}.tsv
    $touch_mpileup
    """
}
