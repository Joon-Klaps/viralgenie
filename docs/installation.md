---
hide:
  - navigation
---
# Installation

Viralgenie uses Nextflow, and a [package/container management system](https://www.nextflow.io/docs/latest/container.html#containers) ([Docker](https://www.docker.com/resources/what-container/), [Singularity](https://docs.sylabs.io/guides/latest/user-guide/introduction.html) or [Conda](https://docs.conda.io/en/latest/)) so both need to be installed on the system where you launch your analysis.

!!! Tip "New to bioinformatics?"
    If the word "terminal" brings to mind an airport boarding area, you can become a little lost. [This blog post](https://www.nextflow.io/blog/2021/setup-nextflow-on-windows.html) (up until Configuring an Xserver ...) will help people with little bioinformatic experience set up Nextflow and Docker on a Windows computer.

## Software managers: Docker, Singularity, and Conda

Viralgenie can be run using either [Docker](https://www.docker.com/resources/what-container/), [Singularity](https://docs.sylabs.io/guides/latest/user-guide/introduction.html) or [Conda](https://docs.conda.io/en/latest/). The choice of container system is up to the user, but it is important to note that Docker and Singularity are the most reproducible. Nextflow supports more containers in addition to Docker and Singularity, such as Podman, Shifter, and Charliecloud. You can read the full list of supported containers and how to set them up [here](https://www.nextflow.io/docs/latest/container.html#containers).

When using these containers, Nextflow will use the manager for each process that is executed. In other words, Nextflow will be using `docker run` or `singularity exec` without the need for you to do anything else.

=== "Docker"

    Docker is a containerization system that allows you to package your code, tools, and data into a single image that can be run on most operating systems. It is the most widely used containerization system in bioinformatics.

    To install Docker, follow the instructions on the [Docker website](https://docs.docker.com/get-docker/).

    !!! warning
        Docker requires root access to run. If you do not have root access like, i.e. a user on an HPC or on a cloud - use Singularity instead.

=== "Singularity"

    Singularity is a containerization system that allows you to package your code, tools, and data into a single image that can be run on most operating systems. It is the most widely used containerization system in bioinformatics.

    To install Singularity, follow the instructions on the [Singularity website](https://sylabs.io/guides/latest/user-guide/quick_start.html).

    !!! warning
        Singularity is a great alternative to Docker but can be challenging to set up on an Apple silicon chip or any other ARM device. If you are using an ARM device, consider using Docker instead.

=== "Conda | Mamba"

    Conda is a package manager that allows you to install software packages and dependencies in isolated environments. It is a good choice if you are facing issues while installing Docker or Singularity.

    - To install Conda, follow the instructions on the [Conda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
    - To install Mamba, a faster alternative to Conda, follow the instructions on the [Mamba miniforge website](https://github.com/conda-forge/miniforge/blob/main/README.md#install).

    !!! warning Container systems are better than Conda
        Conda environments are great! However, Conda tools can easily become broken or incompatible due to dependency issues. For this reason, Conda is not as reproducible as Docker or Singularity containers. If you encounter issues with Conda, please try running the pipeline with Docker or Singularity first to see if the issue persists. In other words, if you have a container system, use it over Conda!

## Nextflow

Nextflow runs on most POSIX systems (Linux, macOS, etc) and requires Java 11 or later. It can be installed in several ways, including using the [Nextflow installer](https://www.nextflow.io/docs/latest/getstarted.html#installation) or [Bioconda](https://bioconda.github.io/).

=== "Nextflow installer"
    !!! Tip
        Unsure how to install Nextflow with these commands? Check out the [Nextflow installation documentation](https://www.nextflow.io/docs/latest/getstarted.html#installation) for more information.

    ```console
    # Make sure that Java v11+ is installed:
    java -version

    # Install Nextflow
    curl -fsSL get.nextflow.io | bash

    # Try a simple demo
    ./nextflow run hello
    ```

    !!! Tip
        Add Nextflow binary to your user's `PATH`:
        ```console
        mv nextflow ~/bin/
        ```
        Or to install it system-wide:
        ```console
        sudo mv nextflow /usr/local/bin/
        ```

=== "Bioconda"

    First, set up Bioconda according to the [Bioconda documentation](https://bioconda.github.io/#usage), notably setting up channels:
    ```console
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```

    A best practice with Conda is to create an environment and install the tools in it.
    Therefore you will prevent version conflicts and keep everything clean.
    To do so use the following command:

    ```bash
    conda create --name nextflow-env nextflow
    conda activate nextflow-env
    ```

    To deactivate the Conda environment, run the following command:

    ```bash
    conda deactivate
    ```

    > If you're already in the Conda environment you want to use, you can just install Nextflow directly:
    >
    > ```bash
    > conda install nextflow
    > ```

## Viralgenie

If you have both Nextflow and a software manager installed, you are all set! You can test the pipeline using the following command:

```bash
nextflow run Joon-Klaps/viralgenie \
    -profile test,<docker/singularity/.../institute>
```
!!! note
    With the argument `-profile <docker/singularity/.../institute>`, you can specify the container system you want to use. The `test` profile is used to run the pipeline with a small dataset to verify if everything is working correctly.

!!! Tip "Running Nextflow on a High performance computing (HPC) system?"
    You might not be the first person to run a Nextflow pipeline on your infrastructure! Check out the [nf-core configuration website](https://nf-co.re/configs) as it might already contain a specific configuration for your infrastructure.

!!! danger "Apple silicon (ARM)"
    If you are using an Apple silicon (ARM) machine, you may encounter issues. Most tools are not yet compatible with ARM architecture, therefore Conda will most likely fail. In this case, use Docker in combination with the profile `arm`.
    ```bash
    nextflow run Joon-Klaps/viralgenie \
        -profile test,docker,arm
    ```
    If you still encounter issues, you can set up a Nextflow Tower account and run the pipeline with [wave containers](https://www.nextflow.io/docs/latest/wave.html). In this config file, supply the following:
    ```groovy
    wave {
        enabled = true
        wave.strategy = ['dockerfile']
    }
    tower {
        accessToken = '<your access token>'
    }
    ```

