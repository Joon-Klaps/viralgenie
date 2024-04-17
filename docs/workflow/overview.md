## Workflow

Viralgenie takes in a set of reads and performs 5 major analyses, each of them are explained in more detail in the following sections:

1. [Preprocessing](preprocessing.md)
2. [Metagenomic diversity](metagenomic_diversity.md)
3. [Assembly & Polishing](assembly_polishing.md)
4. [Variant analysis & iterative refinement](variant_and_refinement.md)
5. [Consensus evaluation](consensus_qc.md)

By default all analyses are run.

!!! Tip "Skipping steps"
    All steps can be skipped and the pipeline can be run with only the desired steps. This can be done with the `--skip_preprocessing`, `--skip_metagenomic_diversity`, `--skip_assembly`, `--skip_polishing`, `--skip_variant_analysis`, `--skip_iterative_refinement`, `--skip_consensus_qc` flags.


## Subway map

![viralgenie-workflow](../images/metromap_style_pipeline_workflow_viralgenie.png)

--8<-- "README.md:31:67"
