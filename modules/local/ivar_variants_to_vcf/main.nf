//taken from https://github.com/nf-core/viralrecon/blob/master/modules/local/ivar_variants_to_vcf.nf
process IVAR_VARIANTS_TO_VCF {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e5b77713d1c484a5d154f87cc3b491ad175787cba119febd6ee0f36b61192a/data'
        : 'community.wave.seqera.io/library/python_pip_biopython_pandas_scipy:8c61b990627f38e4'}"

    input:
    tuple val(meta), path(tsv), path(fasta)
    path header

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // This script is bundled with the pipeline, in viralmetagenome/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar_variants_to_vcf.py \\
        ${tsv} \\
        ${prefix}.vcf \\
        --fasta ${fasta} \\
        ${args} \\
        > ${prefix}.variant_counts.log

    cat ${header} ${prefix}.variant_counts.log > ${prefix}.variant_counts_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.vcf
    touch ${prefix}.log
    """
}
