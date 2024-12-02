process CUSTOM_MPILEUP {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8a/8a45adeea30a99d5e2b3775755942b42f32b5c0f2c7ebc91feb3ec5b3467ade4/data':
        'community.wave.seqera.io/library/argparse_pysam_pysamstats_pandas:986069f6be2d09b3' }"

    input:
    tuple val(meta), path(bam), path(ref)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    custom_mpileup.py \\
        $args \\
        --alignment ${bam} \\
        --reference ${ref} \\
        --output ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version: | sed 's/Version: //g')
        pysam: \$(pip show pysam | grep Version: | sed 's/Version: //g')
        pysamstats: \$(pip show pysamstats | grep Version: | sed 's/Version: //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(pip show pandas | grep Version: | sed 's/Version: //g')
        pysam: \$(pip show pysam | grep Version: | sed 's/Version: //g')
        pysamstats: \$(pip show pysamstats | grep Version: | sed 's/Version: //g')
    """
}
