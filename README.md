[![Test pipeline CI](https://github.com/VIB-PSB/ksrates/actions/workflows/test_pipeline.yml/badge.svg)](https://github.com/VIB-PSB/ksrates/actions/workflows/test_pipeline.yml)
[![Push DockerHub CI](https://github.com/VIB-PSB/ksrates/actions/workflows/push_container.yml/badge.svg)](https://github.com/VIB-PSB/ksrates/actions/workflows/push_container.yml)
[![Documentation Status](https://readthedocs.org/projects/ksrates/badge/?version=latest)](https://ksrates.readthedocs.io/en/latest/?badge=latest)


VIB-UGent Center for Plant Systems Biology&mdash;[Evolutionary Systems Biology Lab](http://www.psb.ugent.be/esb)

# *ksrates*
*ksrates* is a tool to position whole-genome duplications\* (WGDs) relative to speciation events using substitution-rate-adjusted mixed paralog&ndash;ortholog distributions of synonymous substitutions per synonymous site (*K*<sub>S</sub>).

*\*or, more generally, whole-genome multiplications (WGMs), but we will simply use the more common WGD to refer to any multiplication*

## Overview

To position ancient WGD events with respect to speciation events in a phylogeny, the *K*<sub>S</sub> values of WGD paralog pairs in a species of interest are often compared with the *K*<sub>S</sub> values of ortholog pairs between this species and other species. For example, it is common practice to superimpose ortholog and paralog *K*<sub>S</sub> distributions in a mixed plot. However, if the lineages involved exhibit different substitution rates, such direct naive comparison of paralog and ortholog *K*<sub>S</sub> estimates can be misleading and result in phylogenetic misinterpretation of WGD signatures. 

*ksrates* is user-friendly command-line tool and [Nextflow](https://github.com/nextflow-io/nextflow) pipeline to compare paralog and ortholog *K*<sub>S</sub> distributions derived from genomic or transcriptomic sequences. *ksrates* estimates differences in synonymous substitution rates among the lineages involved and generates an adjusted mixed plot of paralog and ortholog *K*<sub>S</sub> distributions that allows to assess the relative phylogenetic positioning of presumed WGD and speciation events.

## Useful links

[Documentation](https://ksrates.readthedocs.io/)\
[Tutorial](https://ksrates.readthedocs.io/en/latest/usage.html)\
[FAQ](https://ksrates.readthedocs.io/en/latest/faqs.html)\
[Publication](https://doi.org/10.1093/bioinformatics/btab602)

## Quick start

*ksrates* can be executed using either a [Nextflow](https://github.com/nextflow-io/nextflow) pipeline (recommended) or a manual command-line interface. The latter is available via Docker and Apptainer containers, and as a Python package to be integrated into existing genomics toolsets and workflows. 

In the following sections we briefly describe how to install, configure and run the Nextflow pipeline and the basic usage of the command-line interface for the Docker or Apptainer containers. For detailed usage information, full tutorial and additional installation options, please refer the [full documentation](https://ksrates.readthedocs.io/).

### Example and test datasets

To illustrate how to use *ksrates*, a simple use case is provided for the analysis of WGD signatures in monocot plants, with oil palm (*Elaeis guineensis*) as the focal species:

- [`example`](example): a full dataset which contains the complete sequence data for the involved species and requires hours of computations depending on the available computing resources. We advise to run this dataset on a computer cluster with the *ksrates* Nextflow pipeline, which supports several types of HPC schedulers.

- [`test`](test): a test dataset that contains a small subset of the sequence data and takes only a few minutes to be run (also locally on a laptop). This is intended for a quick check of the tool installation or to get familiar with the tool, as the output is not scientifically meaningful.

See the *Usage* sections below and the [Tutorial](https://ksrates.readthedocs.io/en/latest/usage.html) for more detail.

### Nextflow pipeline

Instructions to execute the Nextflow pipeline.

#### Installation

1. Install [Nextflow](https://github.com/nextflow-io/nextflow), official instructions are [here](https://www.nextflow.io/docs/latest/getstarted.html), but briefly:
	1. If you do not have [Java](https://www.oracle.com/java/) installed, install [Java 11 or later](https://www.oracle.com/java/technologies/javase-downloads.html); on Linux you can use: 
	
           sudo apt-get install default-jdk

	2. Install Nextflow using either:
	
           wget -qO- https://get.nextflow.io | bash

	   or:
	
	       curl -fsSL https://get.nextflow.io | bash

	   It creates the `nextflow` executable file in the current directory. You may want to move it to a folder accessible from your `$PATH`, for example:
	    
           mv nextflow /usr/local/bin
    
2. Install either [Apptainer](https://apptainer.org/docs/admin/latest/installation.html) (recommended, but check [availability](https://ksrates.readthedocs.io/en/latest/installation.html#availability-and-dependencies)) or [Docker](https://docs.docker.com/get-docker/).

3. The *ksrates* Apptainer or Docker container will be automatically downloaded when you execute the *ksrates* Nextflow pipeline for the first time, and it will be stored and reused for any further executions (see [Nextflow pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html)). *ksrates* and its dependencies are already installed within the container.

#### Usage

We briefly illustrate here how to run the *ksrates* Nextflow pipeline on the `test` dataset.

1. Clone the ksrates repository from GitHub in order to access the test dataset:

           git clone https://github.com/VIB-PSB/ksrates
           cd ksrates/test

2. Inspect and possibly adapt the configuration files. The `test` directory contains:
      
      * A pre-filled *ksrates* configuration file (`config_files/config_elaeis.txt`).

      * A pre-filled *ksrates* expert configuration file (`config_files/config_expert.txt`).

      * A Nextflow configuration file template (`nextflow.config`). While the test dataset is very lightweight, real-case scenarios require configuration of this file: pipeline executor (i.e. local computer or computer cluster), allocated resources (e.g. number of CPUs) and *ksrates* Apptainer or Docker container.

        See the [full documentation](https://ksrates.readthedocs.io/) and the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more details. A general Nextflow configuration template is provided in the [docs](docs) directory on the GitHub repository.


3. Launch the *ksrates* Nextflow pipeline.

   > **Note:** If this is the first time you launch the pipeline, Nextflow will first download the *ksrates* Nextflow pipeline and the *ksrates* Apptainer or Docker container.
       
	      nextflow run VIB-PSB/ksrates --test -profile apptainer --config config_files/config_elaeis.txt --expert config_files/config_expert.txt

   The *ksrates* configuration file is specified through the `--config` parameter, while the *ksrates* expert configuration file is specified through the `--expert` parameter. A Nextflow configuration file named `nextflow.config` and located in the launching folder is automatically detected. The `--test` parameter is mandatory for the test dataset.

   > **Note:** To generate a new *ksrates* configuration file template for a new analysis, use the `--config` option to specify its file name or file path. If the specified file does not exist, the command will generate a configuration file template and exit. Edit the template (see the [full documentation](https://ksrates.readthedocs.io/) for more details) and then rerun the same command above to relaunch the pipeline.

### Command-line interface

Instructions to execute the manual command-line interface.

#### Installation
   
1. Install either [Apptainer](https://apptainer.org/docs/admin/latest/installation.html) (recommended, but check [availability](https://ksrates.readthedocs.io/en/latest/installation.html#availability-and-dependencies)) or [Docker](https://docs.docker.com/get-docker/), then pull the container from VIB-PSB Docker Hub. The container contains *ksrates* and all other software dependencies, so that nothing else needs to be installed.

	For the Apptainer container:
	    
        apptainer pull docker://vibpsb/ksrates
	    
	Or for the Docker container:
	    
        docker pull vibpsb/ksrates:latest

2. *ksrates* and its dependencies are already installed within the container.

#### Usage

Syntax to execute an individual command with command-line interface:

* Basic *ksrates* syntax:

      ksrates [OPTIONS] COMMAND [ARGS]

* Execute a *ksrates* command using the Apptainer container:

      apptainer exec docker://vibpsb/ksrates ksrates [OPTIONS] COMMAND [ARGS]

* Execute a *ksrates* command using the Docker container:

      docker run --rm -v $PWD:/temp -w /temp vibpsb/ksrates ksrates [OPTIONS] COMMAND [ARGS]
	

Example *ksrates* commands with basic syntax:

* Show help menu for usage and all available `COMMAND`s and `OPTIONS`:

      ksrates -h

* Generate a template configuration file for the focal species:

      ksrates generate-config config_files/config_elaeis.txt

* Show help menu for usage and `ARGS` of a specific `COMMAND`:

      ksrates orthologs-ks -h

* Run the ortholog *K*<sub>S</sub> analysis between two species specifying the number of CPUs:

      ksrates orthologs-ks config_files/config_elaeis.txt --expert config_files/config_expert.txt elaeis oryza --n-threads 4

Please see the [full documentation](https://ksrates.readthedocs.io/) for more details and for the complete pipeline.


## Support

If you come across a bug or have any question or suggestion, feel free to reach out on our GitHub [issue](https://github.com/VIB-PSB/ksrates/issues) page.


## Citation

If you publish results generated using *ksrates*, please cite:

Sensalari C., Maere S. and Lohaus R. (2021) *ksrates*: positioning whole-genome duplications relative to speciation events in *K*<sub>S</sub> distributions. *Bioinformatics*, btab602, [doi: https://doi.org/10.1093/bioinformatics/btab602](https://doi.org/10.1093/bioinformatics/btab602)
