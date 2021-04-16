Usage
*****

This section illustrates how to run *ksrates* on the use case dataset proposed in the :ref:`explained_example`, where the rate-adjustment is relative to the focal species oil palm (*Elaeis guineensis*). The use case dataset is stored in the GitHub repository under the ``example`` directory. The pipeline steps can be run either through Nextflow (recommended) or manually. In either case it is advised to use a computing cluster. 

Clone the GitHub repository to get the use case dataset::

    git clone https://github.com/VIB-PSB/ksrates


Run example case as a Nextflow pipeline (recommended)
=====================================================

The *ksrates* pipeline can be automatically run through Nextflow with a few preparation steps.

1.  Access in a terminal the directory that will host the rate-adjustment results (assumed here to be ``example``) and unzip the sequence data files in there:: ::

        cd ksrates/example
        gunzip elaeis.fasta.gz oryza.fasta.gz asparagus.fasta.gz elaeis.gff3.gz

2.  Prepare the configuration files.

    The directory already contains a pre-filled *ksrates configuration file* (``config_elaeis.txt``) and a *Nextflow configuration file* template (``custom_nextflow.config``). If running on a cluster or within a container, please fill in the template as described in the :ref:`nextflow_config_section` section.

    .. note ::
        When running *ksrates* on a new dataset, the configuration files still have to be generated.
        
        For the *Nextflow configuration file*, please either consult the Nextflow documentation or check the templates available in the *ksrates* GitHub repository under ``doc/source`` (Singularity-targeted, Docker-targeted or container-independent templates).

        To generate a new *ksrates configuration file*, launch the pipeline (step 3 below) specifying a non-existing filename after the ``--config`` option. Not finding the file, the code produces a template to be filled in as described in :ref:`pipeline_config_section` section. After that, repeat step 3 again.

3.  Launch *ksrates* through the following command line::

        nextflow run VIB-PSB/ksrates --config ./config_elaeis.txt -c ./custom_nextflow.config

    The ``--config`` option takes the *ksrates configuration file*, while ``-c`` takes the optional *Nextflow configuration file*. If the Nextflow-reserved ``nextflow.config`` name is used, this latter file is automatically recognized without explicitly calling it in the command line.
    
    The first time the command is launched it downloads the *ksrates* Nextflow pipeline from the ``VIB-PSB/ksrates`` GitHub repository; from then on it uses the local copy stored in the ``.nextflow`` directory. If running a container, the image is pulled from Docker Hub and stored locally for successive usage.  


.. _`manual_pipeline`:

Run example case as a manual pipeline
=====================================

The pipeline can otherwise be run by manually launching all the individual commands which it is composed of. This also allows to re-execute single desired steps.

The syntax to run a command depends on how the package is installed:

*   Local installation:: 

        ksrates [OPTIONS] COMMAND [ARGS]

*   Singularity container:

        Open an interactive container where to launch commands with the syntax indicated in the local installation above::

            singularity shell docker://vibpsb/ksrates

        Or launch a single command through the container::

            singularity exec docker://vibpsb/ksrates ksrates [OPTIONS] COMMAND [ARGS]

    Singularity downloads the container image from Docker Hub in ``$HOME/.singularity/cache`` and from then on makes use of the local copy.

*   Docker container:

        Open an interactive container where to launch commands with the syntax indicated in the local installation above::

            docker run -it --rm -v $PWD:/temp -w /temp vibpsb/ksrates

        Or launch a single command through the container::

            docker run --rm -v $PWD:/temp -w /temp vibpsb/ksrates ksrates [OPTIONS] COMMAND [ARGS]

    The ``--rm`` option is given to remove the container after the command is executed to save disk space (note that the container *image* will not be removed). The ``-v`` option mounts the current working directory in the container, while ``-w`` lets the command be run within this directory. 

    Docker pulls the container image from Docker Hub and from then on makes use of the local copy.

In order to submit the command as a job on a compute cluster, wrap the command in the appropriate syntax for the cluster executor system/HPC scheduler (e.g. ``qsub`` for a Sun Grid Engine (SGE) or compatible cluster or a PBS/Torque family scheduler). It is strongly recommended to run the *K*:sub:`S` paralog and orthologs estimation steps (see commands below) on a compute cluster.

An overview of the commands is available by accessing the package help menu (``ksrates -h``)::

    generate-config       Generates configuration file.
    init                  Initializes rate-adjustment.
    orthologs-adjustment  Performs ortholog substitution rate-adjustment.
    orthologs-analysis    Computes ortholog divergence times Ks estimates.
    orthologs-ks          Performs ortholog Ks estimation.
    paralogs-analyses     Detects WGD signatures in paralog Ks distribution.
    paralogs-ks           Performs paralog Ks estimation.
    plot-orthologs        Generates ortholog Ks distributions plot.
    plot-paralogs         Generates rate-adjusted mixed Ks plot.
    plot-tree             Generates phylogram with Ks-unit branch lengths.

The order of execution of the single commands to run the whole workflow is the following. We assume here a local installation without the use of a *ksrates* container.

1.  Access in a terminal the directory that will host the rate-adjustment results (assumed here to be ``example``) and unzip the sequence data files in there:: ::

        cd ksrates/example
        gunzip elaeis.fasta.gz oryza.fasta.gz asparagus.fasta.gz elaeis.gff3.gz

2.  The ``example`` directory already contains a pre-filled configuration file (``config_elaeis.txt``).

    .. note ::
        To generate a new configuration file for your own analyses, run the following command and fill in the template as described in :ref:`pipeline_config_section` section::

            ksrates generate-config config_filename.txt

3.  Run the initialization script to obtain the ortholog trios for the rate-adjustment (``rate_adjustment/elaeis/ortholog_trios_elaeis.tsv``) and to extract the species pairs to be run through the *wgd* ortholog *K*:sub:`S` analysis (``rate_adjustment/elaeis/ortholog_pairs_elaeis.txt``)::

        ksrates init config_elaeis.txt

    This step also generates ``wgd_runs_elaeis.txt`` in the launching directory, which lists all the commands to be run in steps 4 and 5. 

4.  Launch the *wgd* paralog *K*:sub:`S` analysis to estimate the whole-paranome *K*:sub:`S` values (``paralogs_distributions/wgd_elaies/elaeis.ks.tsv``) and optionally the anchor pair *K*:sub:`S` values (``paralogs_distributions/wgd_elaies/elaeis.ks_anchors.tsv``)::

        ksrates paralogs-ks config_elaeis.txt --n-threads 4

    Using multiple threads to parallelize the analysis will reduce the compute time. The ``--n-threads`` option configures the number of threads to use (set this according to your available resources, i.e. CPUs/cores; we recommend a value around 10 and thus the use of a compute cluster).

5.  Launch the *wgd* ortholog *K*:sub:`S` analysis to estimate the ortholog *K*:sub:`S` values *for each required species pair*. These are listed in ``rate_adjustment/elaeis/ortholog_pairs_elaeis.txt``::

        ksrates orthologs-ks config_elaeis.txt elaeis asparagus --n-threads 4
        ksrates orthologs-ks config_elaeis.txt elaeis oryza --n-threads 4
        ksrates orthologs-ks config_elaeis.txt oryza asparagus --n-threads 4

    Using multiple threads to parallelize the analysis will reduce the compute time. The ``--n-threads`` option configures the number of threads to use (set this according to your available resources, i.e. CPUs/cores; we recommend a value around 10 and thus the use of a compute cluster).

    The output files are generated in the directory ``ortholog_distributions``, e.g. the first command generates the file ``ortholog_distributions/wgd_asparagus_elaeis/asparagus_elaeis.ks.tsv``. The two species names are in case-insensitive alphabetical order.

6.  Estimate the mode and associated standard deviation for each ortholog *K*:sub:`S` distribution::
    
        ksrates orthologs-analysis config_elaeis.txt

    The results are stored in a local database, namely a TSV file called by default ``ortholog_peak_db.tsv`` and generated by default in the launching directory (see :ref:`pipeline_config_section`).

7.  Plot the ortholog *K*:sub:`S` distributions for each focal species--other species pair (and each of their trios)::
    
        ksrates plot-orthologs config_elaeis.txt

    The command generates a PDF file for each species pair with the three ortholog *K*:sub:`S` distributions obtained from each of the species trios the species pair is involved in. Note that if multiple trios/outgroups exist, the file is a multi-page PDF showing one trio per page. The two species names are in case-insensitive alphabetical order. In this example case there is only the *E. guineensis*--*O. sativa* species pair, thus the correspondent PDF file generated is ``rate_adjustment/elaeis/orthologs_elaeis_oryza.pdf``.
     
8.  Perform the rate-adjustment. **Pre-requisite**: all *wgd* paralog and ortholog *K*:sub:`S` analyses (steps 4 and 5) and ortholog *K*:sub:`S` distribution mode estimates (step 6) must be completed. ::
    
        ksrates orthologs-adjustment config_elaeis.txt

    The branch-specific *K*:sub:`S` contributions and the rate-adjusted ortholog *K*:sub:`S` mode estimates are collected in ``rate_adjustment/elaeis/adjustment_table_elaeis.tsv``.

9.  Plot the adjusted mixed paralog--ortholog *K*:sub:`S` distribution plot (``rate_adjustment/elaeis/mixed_elaeis_adjusted.pdf``)::

        ksrates plot-paralogs config_elaeis.txt
    
10. Plot the phylogram based on the input phylogenetic tree with branch lengths equal to the *K*:sub:`S` distances estimated from the ortholog *K*:sub:`S` distirbutions (``rate_adjustment/elaeis/tree_elaeis_distances.pdf``)::
    
        ksrates plot-tree config_elaeis.txt

11. Plot the adjusted mixed paralog--ortholog *K*:sub:`S` distribution with inferred WGD components::
    
        ksrates paralogs-analyses config_elaeis.txt
    
    The method(s) used for detecting WGD signatures depends on the paralog analysis settings in the *ksrates* configuration file(s): if ``colinearity`` is turned on, the anchor *K*:sub:`S` clustering is performed (``rate_adjustment/elaeis/mixed_elaeis_anchor_clusters.pdf``), otherwise an exponential-lognormal mixture model is performed (``rate_adjustment/elaeis/mixed_species_elmm.pdf``). Additional methods can be executed upon specification in the *ksrates* expert configuration file (``rate_adjustment/elaeis/mixed_species_lmm_paranome.pdf`` and ``rate_adjustment/elaeis/mixed_species_lmm_colinearity.pdf``) (see :ref:`expert_config_section`).
