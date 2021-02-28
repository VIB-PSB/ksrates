Usage
*****

This section illustrates how to run the correction pipeline on the use case dataset proposed in the :ref:`explained_example`, where the rate-correction is relative to oil palm (*Elaeis guineensis*). The use case dataset is stored in the GitHub repository under the ``example`` directory. The pipeline steps can be run either through Nextflow (suggested) or manually. In any case it is advised the use of a computing cluster. 

Clone the GitHub repository to get the use case dataset::

    git clone https://github.com/VIB-PSB/ksrates


Run example case as Nextflow pipeline (suggested)
=================================================

The rate-correction pipeline can be automatically run through Nextflow with few preparation steps.

1.  Access in a terminal the directory that will host the correction results (assumed here to be ``example``). ::

        cd ksrates/example
    
    Unzip the sequence data files::

        gunzip elaeis.fasta.gz oryza.fasta.gz asparagus.fasta.gz elaeis.gff3.gz

2.  Prepare the configuration files.

    The directory already contains a pre-filled *ksrates* configuration file (``config_elaeis.txt``) and a Nextflow configuration file template (``custom_nextflow.config``). If running on a cluster or within a container, please fill in the template as described in the :ref:`nextflow_config_section` section.

    .. note ::
        When running *ksrates* on a new dataset, the configuration files still have to be generated. For the Nextflow configuration file, please either consult Nextflow documentation or check the templates available in *ksrates* GitHub repository under ``doc/source`` (Singularity-oriented, Docker-oriented or container-independent templates).

        To generate a new *ksrates* configuration file, launch the pipeline (step 3) specifying a non-existing filename after the ``--config`` option. Not finding the file, the code produces a template to be filled in as described in :ref:`pipeline_config_section` section. After that, proceed with step 3.

3.  Launch *ksrates* through the following command line::

        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt [-c ./custom_nextflow.config]

    The ``--config`` option takes the *ksrates* configuration file, while ``-c`` takes the optional Nextflow configuration file. If the Nextflow-reserved ``nextflow.config`` name is used, the file is automatically recognized without explicitly calling it in the command line.
    
    The first time the command is launched it downloads the *ksrates* Nextflow pipeline from the ``VIB-PSB/ksrates`` GitHub repository; from then on it uses the local copy stored in the ``.nextflow`` directory.  If running in a container, the image is pulled from Docker Hub and stored locally for successive usage.  


Run example case as manual pipeline
===================================

The pipeline can otherwise be run by manually launching all the individual commands which it is composed of. This also allows to re-execute single desired steps.

The syntax to run a command depends on how the package is installed:

*   Local installation:: 

        ksrates <command> <options> <args>

*   Singularity container::

        singularity exec docker://vibpsb/ksrates ksrates <command> <options> <args>

    Singularity downloads the container image from Docker Hub in ``$HOME/.singularity/cache`` and from then on makes use of the local copy.

*   Docker container::

        docker run --rm -v $PWD:/temp -w /temp vibpsb/ksrates ksrates <command> <options> <args>

    The ``-v`` option mounts the current working directory in the container, while ``-w`` lets the command be run within this directory. The ``--rm`` option is given to remove the container after the command is executed to save disk space (note that the container *image* will not be removed).

    Docker pulls the container image from Docker Hub and from then on makes use of the local copy.

In order to submit the command as a job on a cluster, wrap the command in the appropriate syntax for the executor system (e.g. on a SGE: ``qsub -b y``) and provide enough memory and computational power.

An overview of the commands is available by accessing the package help menu (``ksrates -h``)::

    generate-config       Generates configuration file.
    init                  Initializes rate-adjustment.
    orthologs-analysis    Computes missing ortholog peaks.
    orthologs-correction  Performs substitution rate-adjustment.
    orthologs-ks          Performs ortholog Ks estimate.
    paralogs-analyses     Reconstructs WGM traces in paralog distribution.
    paralogs-ks           Performs paralog Ks estimate.
    plot-orthologs        Plots ortholog distributions used to adjust.
    plot-paralogs         Plots rate-adjusted mixed distribution.
    plot-tree             Plots phylogram with branches equal to rates.

The execution order of the single commands to compose the functional pipeline is the following. Here it is assumed a local installation and it is shown the relative syntax.

1.  Open in a terminal the directory that will host the correction results (assumed here to be ``example``). ::

        cd ksrates/example

2.  The ``example`` directory already contains a pre-filled configuration file (``config_elaeis.txt``).

    .. note ::
        To generate a new configuration file for your own analyses, run the following command and fill in the template as described in :ref:`pipeline_config_section` section::

            ksrates generate-config config_elaeis.txt

3.  Run the initialization script to obtain the ortholog trios for the correction (``ortholog_triplets_elaeis.tsv``) and to extract the species pairs to be directed to the ortholog ``wgd`` analysis (``ortholog_pairs_elaeis.txt``)::

        ksrates init config_elaeis.txt

    This step also generates ``wgd_runs_elaeis.txt`` in the launching directory, which lists all the commands to be run in steps 4 and 5. 

4.  Launch the paralog ``wgd`` analysis to estimate the paranome (and optionally anchor pair) *K*:sub:`S` values. The command generates the output file ``elaeis.ks.tsv`` (and ``elaeis.ks_anchors.tsv``). ::

        ksrates paralogs-ks config_elaeis.txt [--n-threads 4]
   
5.  Launch the ortholog ``wgd`` analysis to estimate the ortholog *K*:sub:`S` values *for each required species pair* listed in ``ortholog_pairs_elaeis.txt``. The command generates the output files in the ``ortholog_distributions`` directory: for example, for species pair palm-rice the output file is named ``elaeis_oryza.ks.tsv``, with names in case-insensitive alphabetic order. ::
 
        ksrates orthologs-ks config_elaeis.txt elaeis oryza [--n-threads 4]
        ksrates orthologs-ks config_elaeis.txt elaeis asparagus [--n-threads 4]
        ksrates orthologs-ks config_elaeis.txt oryza asparagus [--n-threads 4]

6.  Compute the peak of each ortholog *K*:sub:`S` distribution and store it in a local database::
    
        ksrates orthologs-analysis config_elaeis.txt
    
7.  Plot all the ortholog distributions used for the rate-adjustment::
    
        ksrates plot-orthologs config_elaeis.txt
    
8.  Perform the rate-adjustment. *Pre-requisites: all paralog and ortholog pipelines (step 4 and 5) and peak estimates (step 6) must have been already completed.* ::
    
        ksrates orthologs-correction config_elaeis.txt
    
9.  Plot the adjusted mixed distribution::
    
        ksrates plot-paralogs config_elaeis.txt
    
10. Plot the input tree with branch lengths equal to *K*:sub:`S` distances::
    
        ksrates plot-tree config_elaeis.txt

11. Plot the adjusted mixed distribution together with the detected WGM traces::
    
        ksrates paralogs-analyses config_elaeis.txt
    
    The method used for detecting WGM traces depends on the analysis settings in the configuration file: if ``colinearity`` is turned on, then the anchor *K*:sub:`S` clustering is performed, otherwise an exponential-lognormal mixture model is performed. Additional methods can be executed upon specification in the expert configuration file (see :ref:`expert_config_section`).
    
