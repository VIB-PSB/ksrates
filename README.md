VIB-UGent Center for Plant Systems Biology---[Evolutionary Systems Biology Lab](http://www.psb.ugent.be/esb)

# ksrates
*ksrates* is a tool to position whole-genome duplications\* (WGDs) relative to speciation events using substitution-rate-adjusted mixed paralog--ortholog distributions of synonymous substitutions per synonymous site (*K*<sub>S</sub>).

<font size="2">\* or, more generally, whole-genome multiplications (WGMs), but we will simply use the more common WGD to refer to any multiplication</font>

## Quick overview

To position ancient WGD events with respect to speciation events in a phylogeny, the *K*<sub>S</sub> values of WGD paralog pairs in a species of interest are often compared with the *K*<sub>S</sub> values of ortholog pairs between this species and other species. For example, it is common practice to superimpose ortholog and paralog *K*<sub>S</sub> distributions in a mixed plot. However, if the lineages involved exhibit different substitution rates, such direct naive comparison of paralog and ortholog *K*<sub>S</sub> estimates can be misleading and result in phylogenetic misinterpretation of WGD signatures. 

*ksrates* is user-friendly command-line tool and [Nextflow](https://github.com/nextflow-io/nextflow) pipeline to compare paralog and ortholog *K*<sub>S</sub> distributions derived from genomic or transcriptomic sequences. *ksrates* estimates differences in synonymous substitution rates among the lineages involved and generates an adjusted mixed plot of paralog and ortholog *K*<sub>S</sub> distributions that allows to assess the relative phylogenetic positioning of presumed WGD and speciation events.

For more details, see:
***TODO:*** paper link
and the documentation below.

<!--
## Concept

`ksrates` is a package for substitution rate-correction in mixed ortholog and paralog *K*<sub>S</sub> distributions. For more details see (TODO paper link) and documentation (TODO attach link).

Mixed *K*<sub>S</sub> distributions are one of the approaches applied to detect whole-genome multiplications (WGMs) and to locate them in a phylogeny. A mixed plot is composed of ortholog *K*<sub>S</sub> distributions - representing divergence events - overlapped onto paralog *K*<sub>S</sub> distributions - representing the duplication history of a species genome. The relative positions of the ortholog peaks and the WGMs peaks are informative about the order of the depicted evolutionary events, allowing to place the occurrence of a WGMs in a specific branch of the evolutionary history of the species.

The reliability of a mixed plot can be jeopardized in case of (remarkable) substitution rate differences between the involved species. In fact, since the *K*<sub>S</sub> value of a homolog pair depends on the substitution rate of the species, different distributions end up to be built on different *K*<sub>S</sub> scales. A direct overlap of distributions is therefore likely to lead to unreliable interpretations.

The *K*<sub>S</sub> rate-correction package offers a correction procedure that brings all the distributions to a common *K*<sub>S</sub> scale by compensating for the substitution rate differences relatively to one "main" species.
The resulting corrected plot provides a more reliable interpretation of the order of the evolutionary events and WGM placement in a phylogeny.
-->


## Documentation

[General documentation]() ***TODO:*** set link  
[Tutorial]() ***TODO:*** set link  
[FAQ]() ***TODO:*** set link  


## Quick start

*ksrates* can be executed using either a [Nextflow](https://github.com/nextflow-io/nextflow) pipeline (recommended) or a manual command-line interface. The latter is available via Docker and Singularity containers, and as a Python package to integrate into existing genomics toolsets and workflows. 

In the following sections we briefly describe how to install, configure and run the Nextflow pipeline and the basic usage of the command-line interface for the Docker or Singularity containers. For detailed usage information, a full tutorial and additional installation options, please see the [full documentation]()(***TODO:*** set link).

### Example datasets

To illustrate how to use *ksrates*, two example datasets are provided for a simple example use case analyzing WGD signatures in monocot plants with oil palm (*Elaeis guineensis*) as the focal species.

- [`example`](example): a full dataset which contains the complete sequence data for the focal species and two other species and may require hours of computations depending on the available computing resources. We advice to run this dataset on a compute cluster and using the *ksrates* Nextflow pipeline should make it fairly easy to configure this for a variety of HPC schedulers.

- [`test`](test): a small test dataset that contains only a small subset of the sequence data for each of the species from the `example` dataset and takes only a few minutes (***TODO:*** is this true on average?) to be run. This is intended for a quick check of the tool only and can be run locally, e.g. on a laptop. The results are not very meaningful.

See the Usage sections below, and the [Tutorial]()(***TODO:*** set link) for more detail. 

### Nextflow pipeline

#### Installation

1. Install [Nextflow](https://github.com/nextflow-io/nextflow), official instructions are [here](https://www.nextflow.io/docs/latest/getstarted.html), but briefly:
	1. If you do not have [Java](https://www.oracle.com/java/) installed, install [Java 8 or later](https://www.oracle.com/java/technologies/javase-downloads.html).
	2. Install Nextflow using either:
	
            wget -qO- https://get.nextflow.io | bash
	    or:
	
	        curl -fsSL https://get.nextflow.io | bash
	    It creates the `nextflow` executable file in the current directory. You may want to move it to a folder accessible from your `$PATH`, for example:
	    
            mv nextflow /usr/local/bin
    
2. Install either [Singularity](https://sylabs.io/guides/3.7/admin-guide/installation.html) (recommended, but see [here]() ***TODO:*** set link for more detail) or [Docker](https://docs.docker.com/get-docker/). This is needed to run the *ksrates* Singularity or Docker container which contain all other required software dependencies, so nothing else needs to be installed.

3. Install *ksrates*: When using Nextflow, *ksrates* and the *ksrates* Singularity or Docker container will be automatically downloaded simply when you execute the launch of the *ksrates* pipeline for the first time, and they will be stored and reused for any further executions (see [Nextflow pipeline sharing](https://www.nextflow.io/docs/latest/sharing.html)). Therefore, in this case it is not necessary to manually install *ksrates*, simply continue with the Usage section below.

<!--
*ksrates* can be installed locally or shipped through a container. The following instructions assume the use of a container (TODO: for more methods, see link to Read The Doc). 

The machine where the rate-correction will be performed (either local computer or remote cluster) needs to have [Singularity](https://sylabs.io/guides/3.7/admin-guide/installation.html#installing-singularity) or [Docker](https://docs.docker.com/get-docker/) installed.

The container will be downloaded from ``vibpsb/ksrates`` Docker Hub repository when launching the correction pipeline. Successive runs will use the local copy.

Since the rate-correction is implemented as a Nextflow pipeline, the only other dependency that must be installed is Nextflow (for more details see official [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#get-started) installation page):

    sudo apt-get install default-jdk
    wget -qO- https://get.nextflow.io | bash
    sudo mv nextflow /usr/bin
-->

#### Usage

We briefly illustrate here how to run the *ksrates* Nextflow pipeline with the `example` or `test` datasets.

1. Get the example datasets.
    1. Clone the repository to get the example datasets: ***TODO:*** or another data repo?

            git clone https://github.com/VIB-PSB/ksrates
	2. You may want to copy the dataset folder you want to use to another location, for example your home folder, and then change to that folder:
	    
            mv ksrates/example ~
            cd ~/example
        or: 
	    
            mv ksrates/test ~
            cd ~/test
	

2. Launch the *ksrates* Nextflow pipeline.
    (If this is the first time you launch the pipeline, Nextflow will first download *ksrates* and the *ksrates* Singularity or Docker container.)

	1. Running locally on a laptop/desktop:
	
		When using Singularity (recommended): ***TODO:*** check command
	
			nextflow run VIB-PSB/ksrates --config config_elaeis.txt -with-singularity docker://vibpsb/ksrates:latest
		Or when using Docker: ***TODO:*** check command
	
			nextflow run VIB-PSB/ksrates --config config_elaeis.txt -with-docker vibpsb/ksrates:latest
		  
		The required `--config` parameter specifies the (path to the) pipeline configuration file for the *ksrates* analyses to be run. If the specified file does not exist (at the given path) a new template configuration file will be generated and the pipeline exits. Edit and fill in the generated configuration file (see the [full documentation]()(***TODO:*** set link) for more detail) and then rerun the same command above to relaunch the pipeline. 
		Both example dataset folders already contain a pre-filled *ksrates* pipeline configuration file for the oil palm example use case, `config_elaeis.txt`, therefore the above Nextflow command should directly launch the pipeline.

	1. Running on a compute cluster:
	
			nextflow run VIB-PSB/ksrates --config config_elaeis.txt -c custom_nextflow.config
	
		The `--config` parameter is the same as above.
		
	    The `-c` parameter specifies a Nextflow configuration file. This file contains settings to configure the compute cluster to be used and the pipelines resources on it such as number of CPUs and amount of memory. It also now configures whether to use the *ksrates* Singularity or Docker container. 
	    Both example dataset folders already contain a template Nextflow configuration file called `custom_nextflow.config` that can be adapted to your resources. Other general template Nextflow configuration files can be found in the [doc](doc/source) folder in the repository.

    	If the Nextflow configuration file is simply named `nextflow.config`, the configuration file will be automatically recognized and used without having to specify it using the `-c` parameter.

	    Please see the [full documentation]()(***TODO:*** set link) and the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more detail on Nextflow configuration, e.g. for different HPC schedulers.


### Command-line interface

#### Installation
   
1. Install either [Singularity](https://sylabs.io/guides/3.7/admin-guide/installation.html) (recommended, but see [here]() ***TODO:*** set link for more detail) or [Docker](https://docs.docker.com/get-docker/). This is needed to run the *ksrates* Singularity or Docker container which contain *ksrates* and all other required software dependencies, so nothing else needs to be installed.
The *ksrates* Singularity or Docker container will be automatically downloaded simply when you execute a *ksrates* command on the publicly accessible container for the first time, and they will be stored and reused for any further command executions.

<!--
2. Download a *ksrates* container from the VIB-PSB Docker Hub:

	For the Singularity container (recommended): ***TODO:*** check command
	    
        singularity pull docker://vibpsb/ksrates
	    
	Or for the Docker container: ***TODO:*** full command
	    
        ... vibpsb/ksrates:latest
-->

#### Usage

We briefly illustrate here how to run *ksrates* using the Singularity or Docker container.

<!--
1. Get the example datasets.
(See Usage Nextflow pipeline above)

2. Execute individual *ksrates* commands.
-->

*ksrates* comes with a command-line interface. Its basic syntax is:  ***TODO:*** double-check this

	ksrates [OPTIONS] COMMAND [ARGS]...

To execute a *ksrates* command using the Singularity container the syntax is:

	singularity exec docker://vibpsb/ksrates ksrates [OPTIONS] COMMAND [ARGS]...

Or to execute a *ksrates* command using the Docker container the syntax is:

	docker run --rm -v $PWD:/temp -w /temp vibpsb/ksrates ksrates [OPTIONS] COMMAND [ARGS]...
	

Some example *ksrates* commands are:

Show usage and all available `COMMAND`s and `OPTIONS`:

	ksrates -h

Generate a template configuration file for a focal species (name):

	ksrates generate-config config_elaeis.txt

Show usage and `ARGS` for a specific `COMMAND`:

	ksrates orthologs-ks -h

Run the ortholog *K*<sub>S</sub> analysis between two species using four threads/CPU cores:

	ksrates orthologs-ks config_palm.txt elaeis oryza --n-threads 4

Please see the [full documentation]()(***TODO:*** set link) for more details and the complete set of commands.


## Support

If you come across a bug or have any question or suggestion, please [get in touch](). ***TODO:*** set link


## Citation

If you publish results generated using *ksrates*, please cite:

>__*ksrates*: positioning whole-genome duplications relative to speciation events using rate-adjusted mixed paralog--ortholog *K*<sub>S</sub> distributions__  
>Cecilia Sensalari, Steven Maere, and Rolf Lohaus  
>bioRxiv, February 28, 2021  

***TODO:*** add link


<!--
## Acknowledgements

***TODO:*** ...
-->
