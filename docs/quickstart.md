---
hide:
  - navigation
---
# Quick Start

Viralgenie needs two things:

 1. [Nextflow](https://www.nextflow.io/)
 2. [Docker](https://www.docker.com/resources/what-container/), [Singularity](https://docs.sylabs.io/guides/latest/user-guide/introduction.html), or [Conda](https://docs.conda.io/en/latest/)

```bash
# Install nextflow with conda
conda install nextflow
```

Run the pipeline with a small test dataset using Docker containers:
```bash
nextflow run Joon-Klaps/viralgenie \
    -profile test,docker
```

> For a more complete guide on how to set up Nextflow, Docker, Singularity, and Conda, see the [installation guide](installation.md).

## Running viralgenie with your own samples

```bash
nextflow run Joon-Klaps/viralgenie \
    -profile docker \
    --input my_samplesheet.csv
```

An input file contains the paths to the fastq files and the sample names. The input file can be in TSV, CSV, YAML, or JSON format.

=== "TSV"

    ```tsv title="input-samplesheet.tsv"
    sample	fastq_1	fastq_2
    sample1	AEG588A1_S1_L002_R1_001.fastq.gz	AEG588A1_S1_L002_R2_001.fastq.gz
    sample2	AEG588A5_S5_L003_R1_001.fastq.gz
    sample3	AEG588A3_S3_L002_R1_001.fastq.gz	AEG588A3_S3_L002_R2_001.fastq.gz
    ```

=== "CSV"

    ```csv title="input-samplesheet.csv"
    sample,fastq_1,fastq_2
    sample1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
    sample2,AEG588A5_S5_L003_R1_001.fastq.gz,
    sample3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
    ```

=== "YAML"

    ```yaml title="input-samplesheet.yaml"
    - sample: sample1
      fastq1: AEG588A1_S1_L002_R1_001.fastq.gz
      fastq2: AEG588A1_S1_L002_R2_001.fastq.gz
    - sample: sample2
      fastq1: AEG588A5_S5_L003_R1_001.fastq.gz
    - sample: sample3
      fastq1: AEG588A3_S3_L002_R1_001.fastq.gz
      fastq2: AEG588A3_S3_L002_R2_001.fastq.gz
    ```

=== "JSON"

    ```json title="samplesheet.json"
    [
        {
            "sample": "sample1",
            "fastq1": "AEG588A1_S1_L002_R1_001.fastq.gz",
            "fastq2": "AEG588A1_S1_L002_R2_001.fastq.gz"
        },
        {
            "sample": "sample2",
            "fastq1": "AEG588A5_S5_L003_R1_001.fastq.gz"
        },
        {
            "sample": "sample3",
            "fastq1": "AEG588A3_S3_L002_R1_001.fastq.gz",
            "fastq2": "AEG588A3_S3_L002_R2_001.fastq.gz"
        }
    ]
    ```

