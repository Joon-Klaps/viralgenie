# Usage

!!! warning "Please read this documentation on the [Viralgenie website](https://github.io/Joon-Klaps/viralgenie/usage)"

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run Joon-Klaps/viralgenie --input ./samplesheet.csv --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

!!! warning
Do **not** use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

=== "params.json"

    The above pipeline run specified with a params file in yaml format:

    ```bash
    nextflow run Joon-Klaps/viralgenie -profile docker -params-file params.yaml
    ```

    !!! info "`params.yaml` will contain:"
        ```json
        {
            input: "./samplesheet.csv",
            outdir: "./results/",
            host_k2_db: "./databases/kraken2/host",
            mapping_constrains: "./mapping_constrains.tsv",
            cluster_method: "mmseqs-linclust"
            ...
        }
        ```

=== "command line"

    ```bash
    nextflow run Joon-Klaps/viralgenie -profile docker \
        --input ./samplesheet.csv \
        --outdir ./results/ \
        --host_k2_db ./databases/kraken2/host \
        --mapping_constrains ./mapping_constrains.tsv \
        --cluster_method 'mmseqs-linclust' \
        ...
    ```

You can also generate such `YAML`/`JSON` files via [`nf-core launch`](https://nf-co.re/tools#launch-a-pipeline)  if `nf-core` is [installed](https://nf-co.re/tools#installation).
```console
nf-core launch Joon-Klaps/viralgenie
```
!!! Tip
    Use [`nf-core launch`](https://nf-co.re/tools#launch-a-pipeline) if it is the first time running the pipeline to explore all its features and options in an accesible way.

## Samplesheets

### Samplesheets: input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 2-3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet (i.e. if the `fastq_2` column is empty, the sample is assumed to be single-end).

An example samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 5 samples.

```csv title="input-samplesheet.csv"
sample,fastq_1,fastq_2
sample1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
sample2,AEG588A5_S5_L003_R1_001.fastq.gz,
sample3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
```

| Column    | Description                                                                                                                                       |
| --------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name, needs to be unique                                                                                                            |
| `fastq_1` | Full path (_not_ relative paths) to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `fastq_2` | Full path (_not_ relative paths) to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |

### Samplesheets: Mapping constrains

Viralgenie allows the user to specify the mapping constrains for the pipeline, meaning that the specified or all samples will be mapped against the specified references. This is done by providing a samplesheet to parameter `--mapping_constrains`.

An example mapping constrain samplesheet file consisting of references that need to have all samples mapped towards them and a selection may look something like the one below. This is for 5 references, 2 of them being a multi-fasta file.

```tsv title="constrains-samplesheet.tsv"
id	species	segment	samples	sequence	definition
Lassa-L-dataset	LASV	L		LASV_L.multi.fasta	Collection of LASV sequences used for hybrid capture bait design, all publicly availble sequences of the L segment clustered at 99.5% similarity
Lassa-S-dataset	LASV	S	sample1;sample3	LASV_S.multi.fasta	Collection of LASV sequences used for hybrid capture bait design, all publicly availble sequences of the S segment clustered at 99.5% similarity
NC038709.1	HAZV	L	sample1;sample2	L-NC_038709.1.fasta	Hazara virus isolate JC280 segment L, complete sequence.
NC038710.1	HAZV	M		M-NC_038710.1.fasta	Hazara virus isolate JC280 segment M, complete sequence.
NC038711.1	HAZV	S		S-NC_038711.1.fasta	Hazara virus isolate JC280 segment S, complete sequence.
```

| Column       | Description                                                                                         |
| ------------ | --------------------------------------------------------------------------------------------------- |
| `id`         | Reference identifier, needs to be unique'                                                           |
| `species`    | [Optional] Species name of the reference                                                            |
| `segment`    | [Optional] Segment name of the reference                                                            |
| `samples`    | [Optional] List of samples that need to be mapped towards the reference. If empty, map all samples. |
| `sequence`   | Full path (_not_ relative paths) to the reference sequence file.                                    |
| `definition` | [Optional] Definition of the reference sequence file.                                               |

!!! Tip
    The `samples` column is optional, if empty, all samples will be mapped towards the reference.

!!! Tip
    Multi-fasta files can be provided and all reads will mapped to all genomes but stats will not be reported separately in the final report.

### Samplesheets: metadata [optional]

A metadata file can also be provided to the pipeline. This file should contain information about the samples in the samplesheet, such as the condition, timepoint, or any other relevant information. This file should be a comma/tab-separated file with a header row and at least one column with the same name as the `sample` column in the samplesheet.

The file will only be used in the final report and will not affect the pipeline run.

```tsv title="metadata.tsv"
sample	sample_accession	secondary_sample_accession	study_accession	run_alias	library_layout
sample1	SAMN14154201	SRS6189918	PRJNA607948	vero76_Illumina.fastq	PAIRED
sample2	SAMN14154205	SRS6189924	PRJNA607948	veroSTAT-1KO_Illumina.fastq	PAIRED
```

Provide the metadata file with the argument `--metadata`.

## Updating the pipeline

```bash
nextflow pull Joon-Klaps/viralgenie
```

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull Joon-Klaps/viralgenie
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [Joon-Klaps/viralgenie releases page](https://github.com/Joon-Klaps/viralgenie/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

!!! Tip
    If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.


## Core Nextflow arguments

!!! note
    These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).


### The `-profile` parameter

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

!!! info
    We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.


The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

-   `test`
    -   A profile with a complete configuration for automated testing
    -   Includes links to test data so needs no other parameters
-   `docker`
    -   A generic configuration profile to be used with [Docker](https://docker.com/)
-   `singularity`
    -   A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
-   `podman`
    -   A generic configuration profile to be used with [Podman](https://podman.io/)
-   `shifter`
    -   A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
-   `charliecloud`
    -   A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
-   `apptainer`
    -   A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
-   `conda`
    -   A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customize the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default Viralgenie uses containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, Viralgenie provides some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running Viralgenie regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
