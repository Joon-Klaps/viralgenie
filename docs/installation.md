## Installation

Viralgenie uses Nextflow, and a [package/container management system](https://www.nextflow.io/docs/latest/container.html#containers) ([Docker](https://www.docker.com/resources/what-container/), [singularity](https://docs.sylabs.io/guides/latest/user-guide/introduction.html) or [conda](https://docs.conda.io/en/latest/)) so both need to be installed on the system where you launch your analysis.

### Software managers: Docker, singularity, and conda

Viralgenie can be run using either [Docker](https://www.docker.com/resources/what-container/), [singularity](https://docs.sylabs.io/guides/latest/user-guide/introduction.html) or [conda](https://docs.conda.io/en/latest/). The choice of container system is up to the user, but it is important to note that Docker and Singularity are the most reproducible. Nextflow supports more containers in addition to Docker and Singularity, such as Podman, Shifter, and Charliecloud. You can read the full list of supported containers and how to set them up [here](https://www.nextflow.io/docs/latest/container.html#containers).

=== "Docker"

    Docker is a containerisation system that allows you to package your code, tools and data into a single image that can be run on most operating systems. It is the most widely used containerisation system in bioinformatics.

    To install Docker, follow the instructions on the [Docker website](https://docs.docker.com/get-docker/).

    !!! warning
        Docker requires root access to run. If you do not have root access like, i.e. a user on a HPC or on a cloud - use Singularity instead.

=== "Singularity"

    Singularity is a containerisation system that allows you to package your code, tools and data into a single image that can be run on most operating systems. It is the most widely used containerisation system in bioinformatics.

    To install Singularity, follow the instructions on the [Singularity website](https://sylabs.io/guides/latest/user-guide/quick_start.html).


=== "Conda | Mamba"

    Conda is a package manager that allows you to install software packages and dependencies in isolated environments. It is a good choice if you are facing issues while installing Docker or Singularity.

    - To install Conda, follow the instructions on the [Conda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
    - To install Mamba, a faster alternative to Conda, follow the instructions on the [Mamba miniforge website](https://github.com/conda-forge/miniforge/blob/main/README.md#install).

    !!! warning
        Conda environments are not as reproducible as Docker or Singularity containers. If you encounter issues with Conda, please try running the pipeline with Docker or Singularity first to see if the issue persists.



### Nextflow

Nextflow runs on most POSIX systems (Linux, macOS, etc) and requires java 11 or later. It can be installed in several ways, including using the [Nextflow installer](https://www.nextflow.io/docs/latest/getstarted.html#installation) or [Bioconda](https://bioconda.github.io/).
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

    A best practice with conda is to create an environment and install the tools in it.
    Therefore you will prevent version conflicts and keep everything clean.
    To do so use the following command:

    ```bash
    conda create --name nextflow-env nextflow
    conda activate nextflow-env
    ```

    To deactivate the conda environment, run the following command:

    ```bash
    conda deactivate
    ```

    > If you're already in the conda environment you want to use, you can just install Nextflow directly:
    >
    > ```bash
    > conda install nextflow
    > ```

### Viralgenie

If you have both nextflow and a software manager installed, you are all set! You can test the pipeline using the following command:

```bash
nextflow run Joon-Klaps/viralgenie \
    -profile test,<docker/singularity/.../institute> \
```
!!! note
    With the argument `-profile <docker/singularity/.../institute>`, you can specify the container system you want to use. The `test` profile is used to run the pipeline with a small dataset to check if everything is working correctly.
