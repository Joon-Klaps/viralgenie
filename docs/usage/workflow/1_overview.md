# Workflow overview

nf-core/viralmetagenome takes in a set of reads and performs 5 major analyses, each of them are explained in more detail in the following sections:

1. [Preprocessing](2_preprocessing.md)
2. [Metagenomic diversity](3_metagenomic_diversity.md)
3. [Assembly & Polishing](4_assembly_polishing.md)
4. [Variant analysis & iterative refinement](5_variant_and_refinement.md)
5. [Consensus evaluation](6_consensus_qc.md)

By default all analyses are run.

:::tip{title="Skipping steps"}
All steps can be skipped and the pipeline can be run with only the desired steps. This can be done with the `--skip_preprocessing`, `--skip_read_classification`, `--skip_assembly`, `--skip_polishing`, `--skip_variant_calling`, `--skip_iterative_refinement`, `--skip_consensus_qc` flags.
:::

## Subway map

![viralmetagenome-workflow](../../images/metromap_style_pipeline_workflow_viralmetagenome.png)
