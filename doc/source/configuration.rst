.. _`config_sections`:

Configuration files
*******************

*ksrates* requires at the minimum the *ksrates* configuration file that contains all the settings necessary for a given analysis. Optionally, an additional *ksrates* configuration file with expert and debug settings can be specified. If using the *ksrates* Nextflow pipeline, an optional Nextflow configuration file is usually necessary to configure, for example, compute cluster and container settings.


.. _`pipeline_config_section`:

*ksrates* configuration file
============================

Description
-----------

The analysis configuration file is composed of a first section defining the species used, their phylogenetic relationships and their sequence data files, a second section concerning general analysis settings and a third section concerning mainly figure parameters. Below is an example for an analysis with oil palm (*Elaeis guineensis*) as the focal species. ::

    [SPECIES]
    focal_species = elaeis
    newick_tree = ((elaeis, oryza), asparagus);
    latin_names = elaeis:Elaeis guineensis, oryza:Oryza sativa, asparagus:Asparagus officinalis

    fasta_filenames = elaeis:elaeis.fasta, oryza:oryza.fasta, asparagus:asparagus.fasta
    gff_filenames = elaeis:elaeis.gff3

    peak_database_path = ortholog_peak_db.tsv
    ks_list_database_path = ortholog_ks_list_db.tsv


    [ANALYSIS SETTING]
    paranome = yes
    colinearity = yes

    gff_feature = mrna
    gff_attribute = id

    max_number_outspecies = 4
    consensus_peak_for_multiple_outgroups = mean among outgroups


    [PARAMETERS]
    x_axis_max_limit_paralogs_plot = 5
    bin_width_para = 0.1
    y_axis_limit_paralogs_plot = None

    num_iterations = 200
    divergence_colors =  Red, MediumBlue, Goldenrod, Crimson, ForestGreen, Gray, SaddleBrown, Black

    x_axis_max_limit_orthologs_plots = 5
    bin_width_ortho = 0.1

    max_ks_para = 5
    max_ks_ortho = 10

The [SPECIES] section includes:

.. TODO: are spaces tolerated in species name?

* **focal_species**: name of the focal species. A *K*:sub:`S` paralog distribution is generated for this species and its *K*:sub:`S`-scale is used as the rate-adjustment reference. It is advised to use a short name (for example, the genus or family name) or common abbreviation here. For example, "elaeis" or "eguineensis" instead of "Elaeis guineensis".
* **newick_tree**: phylogenetic relationships among the involved species as Newick format (more info on the ETE toolkit `website <http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees>`__). It has to contain the focal species as named in parameter `focal_species`. It is advised to use short names or abbreviations for the other species as well.
* **latin_names**: list of associations between each of the species names in parameter `newick_tree` and their scientific names, which will be used in legends and plot titles. The association is made with a colon (':').
* **fasta_filenames**: list of associations between each of the species names in parameter `newick_tree` and the corresponding paths to their FASTA files. The association is made with a colon (':').
* **gff_filenames**: association between the focal species as named in parameter `focal_species` and the path to the GFF3 file for the focal species (only required for collinearity analysis). The association is made with a colon (':').
* **peak_database_path**: path to the database of ortholog *K*:sub:`S` distribution peaks. If the file is not present yet, it will be automatically generated.
* **ks_list_database_path**: path to the database of ortholog *K*:sub:`S` lists. If the file is not present, it will be automatically generated.

The [ANALYSIS SETTING] section includes:

* **paranome**: whether to build/plot the whole-paranome *K*:sub:`S` distribution of the focal species \[yes/no\]
* **colinearity**: whether to build/plot the anchor pair *K*:sub:`S` distribution of the focal species \[yes/no\]
* **gff_feature**: parsing keyword from the third column of the GFF file (e.g. gene, mrna...). Case insensitive.
* **gff_attribute**: parsing keyword from the ninth column of the GFF file (e.g. id, name...). Case insensitive. 
* **max_number_outspecies**: maximum number of trios/outspecies allowed to adjust a divergent pair; if None, all possible outspecies obtained from the phylogenetic tree will be used to form trios and adjust the pair. For more details see below. [Default: 4]
* **consensus_peak_for_multiple_outgroups**: when a divergent pair is adjusted by two or more outgroups, it is possible to get a consensus value considering either the mean among all the rate-adjustments for that pair ("mean among outgroups") or only one rate-adjustment, mostly likely coming from a close and slowly evolving species ("best outgroup"). [Default: "mean among outgroups"]

The [PARAMETERS] section includes:

* For the mixed *K*:sub:`S` distributions plot

    * **x_axis_max_limit_paralogs_plot**: highest value of the x axis in the mixed distribution plot [Default: 5]
    * **bin_width_para**: bin width in paralog *K*:sub:`S` distribution histogram. By default there are ten bins per unit [Default: 0.1]
    * **y_axis_limit_paralogs_plot**: customized highest value of the y axis in the mixed plot [Default: None]
    
* For ortholog divergence *K*:sub:`S`

    * **num_iterations**: number of bootstrap iterations for mode/median estimate [default: 200]
    * **divergence_colors**: list of colors assigned to the divergence nodes: all divergence lines coming from the same divergence node share the same color [Default: 8 colors]
    
* For the ortholog *K*:sub:`S` distribution plots

    * **x_axis_max_limit_orthologs_plots**: highest value of the x axis in the ortholog distribution plots [Default: 5]
    * **bin_width_ortho**: bin width in ortholog *K*:sub:`S` distribution histogram. By default there are ten bins per unit [Default: 0.1]
    
* *K*:sub:`S` value thresholds

    * **max_ks_para**: maximum value accepted for paralog *K*:sub:`S` from data table [Default: 5]
    * **max_ks_ortho**: maximum value accepted for ortholog *K*:sub:`S` from data table [Default: 10]


Guidelines to set the maximum number of outgroups per rate-adjustment
---------------------------------------------------------------------

``max_num_outspecies`` is a parameter used to limit the amount of outgroup species used to adjust a species pair; without that, all possible outgroups would be taken. Having multiple rate-adjustments on the same divergence can provide stronger support for the rate-adjusted plot and is therefore advised to adjust with at least 3 or 4 outgroups to have more reliable results.

However, the more the outgroups, the more the number of ortholog distributions that will have to be computed by the `wgd` ortholog pipeline, which is a quite computationally demanding step. Setting a maximum amount of outgroups lowers the number of rate-adjustments and can therefore save time and resources. It is a good option in case the tree has a complex structure that would collect an unnecessary large number of outgroups or in case the user wants to have a quicker, although somewhat less reliable, result. Note that another option to lower the number of ortholog distributions is to start with a simpler tree structure.

In case ``mean among outgroup`` is set for the consensus rate-adjustment value, it is advised to use at least 3 or better 4 outgroups to adjust a species pair in order to buffer the weight of misleading outliers when computing the mean.


Guidelines to set the consensus method for multiple rate-adjustments
--------------------------------------------------------------------

A consensus value for the rate-adjustment is needed when multiple rate-adjustments are performed for a species pair. The pipeline computes two consensus strategies, but then generates the divergence lines in the mixed plot according to the method specified in the configuration file under ``consensus_peak_for_multiple_outgroups``.

* ``mean among outgroups``: with this option, the final rate-adjustment of a species pair is the mean of the rate-adjustments obtained from all the used outgroups. It is the default method because it avoids to rely on a single voice that could be biased (e.g. bad quality data).
* ``best outgroup``: with this option, only the rate-adjustment obtained from the best outgroup is considered for the final rate-adjustment of a species pair. The best outgroup is the one with the smallest OC segment, which is computed during the relative rate detection. The OC segment is a combined measure of how close is the outgroup and how slow is its rate; the smaller the OC segment, the better can the outgroup detect the relative rates. The OC is stored in ``correction_table_species.tsv``. If one outgroup shows a remarkably slower OC than the others, then it can be worth it to re-run the pipeline (or just the plotting of the mixed distribution) by setting in the configuration file the ``best outgroup`` method. However, it's first better to check the quality of the rate-adjustment result coming from it, especially if the outgroup species has transcriptome data: its ortholog distributions in ``orthologs_species1_species2.pdf`` should have clear peaks in order to give a reliable rate-adjustment.


.. _`nextflow_config_section`:

Nextflow configuration file
===========================

It is a configuration file used to set the communication with the cluster system, the use of a container and to define parameters or variables for the Nextflow pipeline. For a more complete description please refer to `Nextflow documentation <https://www.nextflow.io/docs/latest/config.html#configuration>`_. The user can download a configuration file template from the GitHub repository documentation and adapt it according to their resources and requirements. Below is explained the basic file structure::

    singularity {
        enabled = true
        cacheDir = ''
    }
    docker.enabled = true

    executor.name = ''

    process {
        container = ''

        withName: 'processName' {
            clusterOptions = ''
            beforeScript = ''
        }
    }

    env.SOME_ENV_VARIABLE = ''

* The **singularity** and **docker** scopes deal with container-related specifications:

    * **enable** enables or disables the use of a container
    * **cacheDir** defines the directory where to download and store the Singularity image file from Docker Hub

* The **executor** scope defines the cluster system type (e.g. SGE) which the jobs are submitted to
* The **process** scope defines the container image and the pipeline configuration on the cluster:

    * **container** defines the *ksrates* container image (from Docker Hub or from a local copy if already downloaded).

        * to pull a Singularity container from Docker Hub: ``docker://vibpsb/ksrates:latest``
        * to pull a Docker container from Docker Hub: ``vibpsb/ksrates:latest``

    * **withName** defines settings for individual processes in the Nextflow pipeline; ``wgdParalogs`` and ``wgdOrthologs`` are the most computationally demanding and it is advised to assign them a higher computational power than the other processes.
    * **clusterOption** defines cluster options (allocated memory, number of threads...)
    * **beforeScript** can be used to load required dependencies in the cluster; it is necessary only if the container is not available, provided that the cluster has all dependencies installed

* The **env** scope defines variables exported in the workflow environment


.. _`expert_config_section`:

Expert configuration file
=========================

It is an optional configuration file containing expert parameters for fine-tuning the analysis or for development purposes. The file can be generated using the following template and it is automatically detected when launching the command line (it must be called `config_expert.txt`). ::

    [EXPERT PARAMETERS]
    
    logging_level = info
    peak_stats = mode
    kde_bandwidth_modifier = 0.4
    plot_correction_arrows = no
    max_mixture_model_iterations = 300
    num_mixture_model_initializations = 10
    extra_paralogs_analyses_methods = no
    max_mixture_model_components = 5
    max_ks_for_mixture_model = 5
    max_gene_family_size = 200

* **logging_level**: the logging message level to be shown in the screen (critical, error, warning, info, debug, notset) [Default: info]
* **peak_stats**: the statistics measure that is used to get a representative peak *K*:sub:`S` value of an ortholog distribution or of an anchor *K*:sub:`S` cluster (options: mode or median) [Default: mode]
* **kde_bandwidth_modifier**: modifier to adjust the fitting of the KDE curve on the underlying paranome or anchor *K*:sub:`S` distribution. The kde Scott's factor computed by SciPy tends to produce an overly smooth KDE curve, especially with steep WGD peaks, and therefore it is reduced by multiplying it by a modifier. Decreasing the modifier leads to tighter fits, increasing it leads to smoother fits and setting it at 1 gives the default kde factor. Note that a too small factor is likely to take into account data noise [Default: 0.4]
* **plot_correction_arrows**: flag to turn on or off the presence of rate-adjustment arrows, which start from the original ortholog peak position and end on the rate-adjusted position
* **max_mixture_model_iterations**: maximum number of EM iterations during mixture modeling [Default: 300] 
* **num_mixture_model_initializations**: number of times the EM algorithm is initialized (either for the random initialization in exp-log mixture model or for k-means in lognormal mixture model)
* **max_mixture_model_components**: maximum number of components considered during the execution of mixture models
* **max_ks_for_mixture_model**: upper limit for the Ks range considered during the execution of mixture models 
* **max_gene_family_size**: maximum number of members in a paralog gene family to be taken into account during Ks estimate (larger families will probably increase the computation time, but they may also provide a significant contribute for the Ks distribution) [Default: 200]