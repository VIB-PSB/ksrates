Usage
*****

This section illustrates how to run *ksrates* on the use case dataset proposed in the :ref:`explained_example`, where the rate-adjustment is relative to the focal species oil palm (*Elaeis guineensis*). The use case dataset is stored in the GitHub repository under the ``example`` directory. The pipeline steps can be run either through Nextflow (recommended) or manually. In either case it is advised to use a computing cluster. 

.. note::
    WSL2 users can enter the Windows file system from the terminal through e.g. ``cd mnt/c/Users/your_username``.


.. _`nextflow_pipeline`:

Run example case as a Nextflow pipeline (recommended)
=====================================================

The *ksrates* pipeline can be automatically run through Nextflow with a few preparation steps.

1.  Clone the GitHub repository to get the ``example`` dataset, access the subdirectory and unzip the sequence data files in there::

        git clone https://github.com/VIB-PSB/ksrates
        cd ksrates/example
        gunzip sequences/*

2.  Prepare the configuration files.

    The directory already contains a pre-filled *ksrates configuration file* for focal species ``elaeis`` (``config_files/config_elaeis.txt``), a pre-filled *ksrates expert configuration file* (``config_files/config_expert.txt``) and a *Nextflow configuration file* template (``nextflow.config``) to be filled in as described in the :ref:`nextflow_config_section` section. For more details, refer to the  :ref:`config_sections` section.

    .. note ::
        To generate a new *ksrates configuration file* for your own analyses, launch the pipeline (step 3 below) specifying the desired non-existing filename after the ``--config`` option. By not finding the file, the code produces a template to be filled in as described in :ref:`pipeline_config_section` section. After that, repeat step 3 again.

3.  Launch *ksrates* through the following command line::

        nextflow run VIB-PSB/ksrates -profile apptainer --config config_files/config_elaeis.txt --expert config_files/config_expert.txt

    .. note::
       As from `ksrates` ``v2.0.0``, the Nextflow pipeline has been ported to DSL2 syntax and requires at least Nextflow version ``22.03.0-edge``. Refer to the :ref:`installation page <install_nextflow>` to learn how to get the latest Nextflow version. You can also launch a specific (e.g. previous) Nextflow version through the ``NXF_VER`` environmental variable in the command line::

            NXF_VER=24.10.5 nextflow run VIB-PSB/ksrates <args>

    The first time the command is executed, Nextflow downloads a local copy of the *ksrates* Nextflow pipeline from the ``VIB-PSB/ksrates`` GitHub repository and stores it in the ``$HOME/.nextflow`` directory.
    Parameter ``-profile`` specifies which container will be pulled from Docker Hub (either Apptainer or Docker).

    .. note::
        Since the Apptainer image is by default stored in the *launching folder* under ``work/singularity``, it is recommended to specify a "centralized" destination path through ``apptainer.cacheDir`` in the :ref:`Nextflow configuration file <nextflow_config_section>`.


    The *ksrates configuration file* is specified through the ``--config`` parameter, while the *ksrates expert configuration file* is specified through the ``--expert`` parameter.
    The *Nextflow configuration file* is automatically detected when named with the Nextflow-reserved ``nextflow.config`` filename and when located in the launching directory; alternatively, the user can provide a custom file by specifying its name or path using the ``-C`` option (see `Nextflow documentation <https://www.nextflow.io/docs/latest/cli.html#hard-configuration-override>`__).


.. _`manual_pipeline`:

Run example case as a manual pipeline
=====================================

The pipeline can otherwise be run by manually launching all the individual commands which it is composed of. This also allows to re-execute single desired steps.

The syntax to run a command depends on how the package is installed:

*   Local installation:: 

        ksrates [OPTIONS] COMMAND [ARGS]

*   Apptainer container:

        Open an interactive container where to launch commands with the syntax indicated in the local installation above::

            apptainer shell docker://vibpsb/ksrates

        Or launch a single command through the container::

            apptainer exec docker://vibpsb/ksrates ksrates [OPTIONS] COMMAND [ARGS]

        .. note::
            WSL2 users need option ``-B`` to mount the Windows file system in the container (e.g. ``-B /mnt/c/Users/your_username``).

    Apptainer downloads the container image from Docker Hub in ``$HOME/.apptainer/cache`` and from then on makes use of the local copy.

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
    orthologs-ks-cleanup  Delete all ortholog BLAST tables.
    paralogs-analyses     Detects WGD signatures in paralog Ks distribution.
    paralogs-ks           Performs paralog Ks estimation.
    paralogs-ks-multi     Performs paralog Ks estimation for all species.
    plot-orthologs        Generates ortholog Ks distributions plot.
    plot-paralogs         Generates rate-adjusted mixed Ks plot.
    plot-tree             Generates phylogram with Ks-unit branch lengths.

The order of execution of the single commands to run the whole workflow is the following. We assume here a local installation without the use of a *ksrates* container.

1.  Clone the GitHub repository to get the ``example`` dataset, access the subdirectory and unzip the sequence data files in there::

        git clone https://github.com/VIB-PSB/ksrates
        cd ksrates/example
        gunzip sequences/*

2.  The directory already contains a pre-filled configuration file for focal species ``elaeis`` (``config_files/config_elaeis.txt``) and a pre-filled expert configuration file (``config_files/config_expert.txt``).

    .. note ::
        To generate a new configuration file for your own analyses, run the following command to produce a template to be filled in as described in :ref:`pipeline_config_section` section::

            ksrates generate-config path/to/config_filename.txt

3.  Run the initialization script to obtain the ortholog trios for the rate-adjustment (``rate_adjustment/elaeis/ortholog_trios_elaeis.tsv``) and to extract the species pairs to be run through the *wgd* ortholog *K*:sub:`S` analysis (``rate_adjustment/elaeis/ortholog_pairs_elaeis.txt``)::

        ksrates init config_files/config_elaeis.txt --expert config_files/config_expert.txt

    This step also generates ``wgd_runs_elaeis.txt`` in the launching directory, which drafts all the commands to be run in steps 4 and 5. 

4.  Launch the *wgd* paralog *K*:sub:`S` analysis to estimate the paralog *K*:sub:`S` values for the focal species::

        ksrates paralogs-ks config_files/config_elaeis.txt --expert config_files/config_expert.txt --n-threads 4

    The output files are generated in the ``paralog_distributions/wgd_elaies`` directory, i.e. ``/elaeis.ks.tsv`` for whole-paranome, ``elaeis.ks_anchors.tsv`` for anchor pairs and ``elaeis.ks_recret_top2000.tsv`` for reciprocally retained gene families.

    Using multiple threads to parallelize the analysis will reduce the compute time. The ``--n-threads`` option configures the number of threads to use (set this according to your available resources, i.e. CPUs/cores; e.g. 10 or more cores running on a compute cluster).

5.  Launch the *wgd* ortholog *K*:sub:`S` analysis to estimate the ortholog *K*:sub:`S` values *for each required species pair*. These are listed in ``rate_adjustment/elaeis/ortholog_pairs_elaeis.txt``::

        ksrates orthologs-ks config_files/config_elaeis.txt --expert config_files/config_expert.txt elaeis asparagus --n-threads 4
        ksrates orthologs-ks config_files/config_elaeis.txt --expert config_files/config_expert.txt elaeis oryza --n-threads 4
        ksrates orthologs-ks config_files/config_elaeis.txt --expert config_files/config_expert.txt oryza asparagus --n-threads 4

    The output files are generated in the ``ortholog_distributions`` directory, e.g. the first command generates file ``wgd_asparagus_elaeis/asparagus_elaeis.ks.tsv``. The two species names are in case-insensitive alphabetical order.

    Using multiple threads to parallelize the analysis will reduce the compute time. The ``--n-threads`` option configures the number of threads to use (set this according to your available resources, i.e. CPUs/cores; e.g. 10 or more cores running on a compute cluster).

6.  Estimate the mode and associated standard deviation for each ortholog *K*:sub:`S` distribution::
    
        ksrates orthologs-analysis config_files/config_elaeis.txt --expert config_files/config_expert.txt

    The results are stored in a local database, namely a TSV file called by default ``ortholog_peak_db.tsv`` and generated by default in the launching directory (see :ref:`pipeline_config_section`).

7.  Plot the ortholog *K*:sub:`S` distributions for each focal species--other species pair (and each of their trios)::
    
        ksrates plot-orthologs config_files/config_elaeis.txt --expert config_files/config_expert.txt

    The command generates a PDF file for each species pair with the three ortholog *K*:sub:`S` distributions obtained from each of the species trios the species pair is involved in. Note that if multiple trios/outgroups exist, the file is a multi-page PDF showing one trio per page. The two species names are in case-insensitive alphabetical order. In this example case there is only the *E. guineensis*--*O. sativa* species pair, thus the correspondent PDF file generated is ``rate_adjustment/elaeis/orthologs_elaeis_oryza.pdf``.
     
8.  Perform the rate-adjustment. **Pre-requisite**: all *wgd* paralog and ortholog *K*:sub:`S` analyses (steps 4 and 5) and ortholog *K*:sub:`S` distribution mode estimates (step 6) must be completed. ::
    
        ksrates orthologs-adjustment config_files/config_elaeis.txt --expert config_files/config_expert.txt

    The branch-specific *K*:sub:`S` contributions and the rate-adjusted ortholog *K*:sub:`S` mode estimates are collected in ``rate_adjustment/elaeis/adjustment_table_elaeis.tsv``.

9.  Plot the adjusted mixed paralog--ortholog *K*:sub:`S` distribution plot (``rate_adjustment/elaeis/mixed_elaeis_adjusted.pdf``)::

        ksrates plot-paralogs config_files/config_elaeis.txt --expert config_files/config_expert.txt
    
10. Plot the phylogram based on the input phylogenetic tree with branch lengths equal to the *K*:sub:`S` distances estimated from the ortholog *K*:sub:`S` distirbutions (``rate_adjustment/elaeis/tree_elaeis_distances.pdf``)::
    
        ksrates plot-tree config_files/config_elaeis.txt --expert config_files/config_expert.txt

11. Plot the adjusted mixed paralog--ortholog *K*:sub:`S` distribution with inferred WGD components::
    
        ksrates paralogs-analyses config_files/config_elaeis.txt --expert config_files/config_expert.txt
    
    The methods used for detecting WGD signatures depend on the paralog analysis settings in the *ksrates* configuration files. For more details please refer to section :ref:`paralogs_analyses`.


Finally, the two following commands are not strictly part of the workflow:

12. Remove all BLAST ``.tsv`` files generated by ``orthologs-ks`` in order to free disk space::

        ksrates orthologs-ks-cleanup path/to/ortholog_distributions --expert config_files/config_expert.txt

    The command only acts within the provided path to the ``ortholog_distributions`` directory, for example removing ``wgd_asparagus_elaeis/asparagus_elaeis.blast.tsv`` and all the other analogous files.

13. Run the paralog *K*:sub:`S` analysis for *all* species provided in the Newick tree, and not only for the focal species::

        ksrates paralogs-ks-multi config_files/config_elaeis.txt --expert config_files/config_expert.txt --n-threads 4
    
    For example it will generate ``paralog_distributions/wgd_asparagus`` and ``paralog_distributions/wgd_oryza`` with all related paralog output files.

Practical considerations
========================

When dealing with large input phylogenies it is useful to know that *ksrates* can be used iteratively, by starting with a small dataset and subsequently adding additional species to finetune the phylogenetic positioning of any hypothesized WGDs.
For such iterative analyses the pipeline can reuse data from previous runs, and will only perform additional calculations on the extended dataset where needed.

When *ksrates* is run, the ortholog *K*:sub:`S` values for each species pair in the input phylogenetic tree and the associated ortholog *K*:sub:`S` modes are stored in a local database.
When the *ksrates* pipeline is subsequently rerun with additional species included in the input phylogeny, *ksrates* will skip the ortholog *K*:sub:`S` calculations for any species pair for which an ortholog *K*:sub:`S` mode has already been stored. The database consists of two tabular files (``ortholog_peak_db.tsv`` and ``ortholog_ks_list_db.tsv``, see :ref:`other_output` for more details) generated/accessed by default in the working directory. A custom path location can be otherwise specified in the :ref:`pipeline_config_section`.

In case a user doesn't want to reuse an existing ortholog *K*:sub:`S` mode of a particular species pair and wants instead to re-estimate it from the same input data but using e.g. a different number of bootstrap iterations or KDE bandwidth, the line concerning the mode has to be manually deleted from the ``ortholog_peak_db.tsv`` database file. The successive *ksrates* pipeline will re-estimate the mode according to the new parameters by starting from the previously computed ortholog *K*:sub:`S` estimates for the species pair concerned, thereby skipping the onerous ortholog *K*:sub:`S` estimation step.
