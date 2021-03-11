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
    gff_filename = elaeis:elaeis.gff3

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

    num_bootstrap_iterations = 200
    divergence_colors =  Red, MediumBlue, Goldenrod, Crimson, ForestGreen, Gray, SaddleBrown, Black

    x_axis_max_limit_orthologs_plots = 5
    bin_width_ortho = 0.1

    max_ks_para = 5
    max_ks_ortho = 10

The [SPECIES] section includes:

.. TODO: are spaces tolerated in species name?

* **focal_species**: name of the focal species. A *K*:sub:`S` paralog distribution is generated for this species and its *K*:sub:`S`-scale is used as the rate-adjustment reference. It is advised to use a short name (for example, the genus or family name) or common abbreviation here. For example, "elaeis" or "eguineensis" instead of "Elaeis guineensis".
* **newick_tree**: phylogenetic relationships among the involved species as Newick format using leaf nodes only (more info e.g. on this `website <https://evolution.genetics.washington.edu/phylip/newicktree.html>`__). It has to contain the focal species as named in parameter `focal_species`. It is advised to use short names or abbreviations for the other species as well.
* **latin_names**: list of associations between each of the species names in parameter `newick_tree` and their scientific names, which will be used in legends and plot titles. The association is made with a colon (':').
* **fasta_filenames**: list of associations between each of the species names in parameter `newick_tree` and the corresponding paths to their FASTA files. The association is made with a colon (':').
* **gff_filename**: association between the focal species as named in parameter `focal_species` and the path to the GFF3 file for the focal species (only required for collinearity analysis). The association is made with a colon (':').
* **peak_database_path**: path to the database of ortholog *K*:sub:`S` distribution peaks. If the file is not present yet, it will be automatically generated.
* **ks_list_database_path**: path to the database of ortholog *K*:sub:`S` lists. If the file is not present, it will be automatically generated.

The [ANALYSIS SETTING] section includes:

* **paranome**: whether to build/plot the whole-paranome *K*:sub:`S` distribution of the focal species (options: "yes" and "no"). [Default: "yes"]
* **colinearity**: whether to build/plot the anchor pair *K*:sub:`S` distribution of the focal species (options: "yes" and "no"). [Default: "no"]
* **gff_feature**: parsing keyword from the third column of the GFF file (e.g. gene, mrna...). Case insensitive.
* **gff_attribute**: parsing keyword from the ninth column of the GFF file (e.g. id, name...). Case insensitive. 
* **max_number_outspecies**: maximum number of trios/outspecies allowed to adjust a divergent pair; if None, all possible outspecies obtained from the phylogenetic tree will be used to form trios and adjust the pair. For more details see below. [Default: 4]
* **consensus_peak_for_multiple_outgroups**: when a divergent pair is adjusted by two or more outgroups, it is possible to get a consensus value considering either the mean among all the rate-adjustments for that pair ("mean among outgroups") or only one rate-adjustment, mostly likely coming from a close and slowly evolving species ("best outgroup"). [Default: "mean among outgroups"]

The [PARAMETERS] section includes:

* For the mixed *K*:sub:`S` distributions plot

    * **x_axis_max_limit_paralogs_plot**: highest value of the x axis in the mixed distribution plot. [Default: 5]
    * **bin_width_para**: bin width in paralog *K*:sub:`S` distribution histogram. By default there are ten bins per unit. [Default: 0.1]
    * **y_axis_limit_paralogs_plot**: customized highest value of the y axis in the mixed plot. [Default: None]
    
* For ortholog divergence *K*:sub:`S`

    * **num_bootstrap_iterations**: number of bootstrap iterations for mode/median estimate. [default: 200]
    * **divergence_colors**: list of colors assigned to the divergence nodes: all divergence lines coming from the same divergence node share the same color. [Default: 8 colors]
    
* For the ortholog *K*:sub:`S` distribution plots

    * **x_axis_max_limit_orthologs_plots**: highest value of the x axis in the ortholog distribution plots. [Default: 5]
    * **bin_width_ortho**: bin width in ortholog *K*:sub:`S` distribution histogram. By default there are ten bins per unit. [Default: 0.1]
    
* *K*:sub:`S` value thresholds

    * **max_ks_para**: maximum value accepted for paralog *K*:sub:`S` from data table. [Default: 5]
    * **max_ks_ortho**: maximum value accepted for ortholog *K*:sub:`S` from data table. [Default: 10]


Guidelines to set the maximum number of outgroups per rate-adjustment
---------------------------------------------------------------------

``max_num_outspecies`` is a parameter used to limit the amount of outgroup species used to adjust a species pair; without that, all possible outgroups would be taken. Having multiple rate-adjustments on the same divergence can provide stronger support for the rate-adjusted plot and is therefore advised to adjust with at least 3 or 4 outgroups to have more reliable results.

However, the more the outgroups, the more the number of ortholog distributions that will have to be computed by the `wgd` ortholog pipeline, which is a quite computationally demanding step. Setting a maximum amount of outgroups lowers the number of rate-adjustments and can therefore save time and resources. It is a good option in case the tree has a complex structure that would collect an unnecessary large number of outgroups or in case the user wants to have a quicker, although somewhat less reliable, result. Note that another option to lower the number of ortholog distributions is to start with a simpler tree structure.

In case ``mean among outgroup`` is set for the consensus rate-adjustment value, it is advised to use at least 3 or better 4 outgroups to adjust a species pair in order to buffer the weight of misleading outliers when computing the mean.


Guidelines to set the consensus method for multiple rate-adjustments
--------------------------------------------------------------------

A consensus value for the rate-adjustment is needed when multiple rate-adjustments are performed for a species pair. The pipeline computes two consensus strategies, but then generates the divergence lines in the mixed plot according to the method specified in the configuration file under ``consensus_peak_for_multiple_outgroups``.

* ``mean among outgroups``: with this option, the final rate-adjustment of a species pair is the mean of the rate-adjustments obtained from all the used outgroups. It is the default method because it avoids to rely on a single voice that could be biased (e.g. bad quality data).
* ``best outgroup``: with this option, only the rate-adjustment obtained from the best outgroup is considered for the final rate-adjustment of a species pair. The best outgroup is the one with the smallest OC segment, which is also computed through principles of the relative rate test. The OC segment is a combined measure of how close is the outgroup and how low is its rate; the smaller the OC segment, the better can the outgroup detect the branch-specific Ks contributions of the two ingroups. The OC is stored in ``adjustment_table_species.tsv``. If one outgroup shows a remarkably smaller OC than the others, then it can be worth it to re-run the pipeline (or just the plotting of the mixed distribution) by setting in the configuration file the ``best outgroup`` method. However, it's first better to check the quality of the rate-adjustment result coming from it, especially if the outgroup species has transcriptome data: its ortholog distributions in ``orthologs_species1_species2.pdf`` should have clear peaks in order to give a reliable rate-adjustment.


.. _`nextflow_config_section`:

Nextflow configuration file
===========================

The Nextflow configuration file is used to configure various settings for the *ksrates* Nextflow pipeline, including the settings to use and configure resources on a compute cluster and to use the *ksrates* Singularity or Docker container. We provide a few general template Nextflow configuration files for the *ksrates* Nextflow pipeline in the `doc <https://github.com/VIB-PSB/ksrates/blob/master/doc/source>`_ directory in the GitHub repository. These can be adapted to a user's specific resources and requirements. Below, we briefly explain some of the basic key settings. For a more complete description please refer to the `Nextflow documentation <https://www.nextflow.io/docs/latest/config.html#configuration>`_. ::

    singularity {
        enabled = true
        cacheDir = ''
    }
    docker {
        enabled = false
    }

    executor.name = ''
								
    process {
        container = ''

        withName: 'processName' {
            clusterOptions = ''
            beforeScript = ''
        }
    }

    env {
    	SOME_ENV_VARIABLE = ''
    	ANOTHER_ENV_VARIABLE = ''
    }

* The **singularity** and **docker** scopes configure container type usage and execution:

    * **enable** enables or disables the use of the respective container
    * **cacheDir** the directory where remote the Singularity image from Docker Hub is stored. When using a computing cluster it must be a shared folder accessible to all computing nodes.

* The **executor.name** setting defines the system type or HPC scheduler to be used (e.g. ``sge``, ``local``)
* The **process** scope defines the configuration for the processes of the *ksrates* pipeline:

    * **container** defines the Singularity or Docker *ksrates* container image to be used (from Docker Hub or from a local copy if already downloaded):

        * to pull a Singularity container from Docker Hub: ``docker://vibpsb/ksrates:latest``
        * to pull a Docker container from Docker Hub: ``vibpsb/ksrates:latest``
        
      Has to match enabled container type (see above).

    * **withName** defines settings for individual processes in the *ksrates* Nextflow pipeline, for example:
    
    	* **clusterOption** any native configuration option accepted by your cluster submit command. You can use it to request non-standard resources or use settings that are specific to your cluster and not supported out of the box by Nextflow.
    	* **beforeScript** allows you to execute a custom (Bash) snippet before the main process script is run. This may be useful to initialise the underlying cluster environment or for other custom initialisation, for example it can be used to load required dependencies if one of the container is not used, provided that the cluster has those dependencies installed.
    	
      For a complete list of available settings, see the `Nextflow documentation <https://www.nextflow.io/docs/latest/process.html#process-directives>`_. The processes ``wgdParalogs`` and ``wgdOrthologs`` are the most computationally demanding and it is advised to assign them a higher computational power than the other processes.

* The **env** scope allows the definition one or more variable that will be exported in the environment where the workflow tasks will be executed.


.. _`expert_config_section`:

Expert configuration file
=========================

This is an optional configuration file that contains several \"expert\" parameters for fine-tuning the analysis or for development/debug purposes. The file has to be named `config_expert.txt` and is then automatically detected when launching *ksrates*. The following can be used as a template::

    [EXPERT PARAMETERS]
    
    logging_level = info
    distribution_peak_estimate = mode
    kde_bandwidth_modifier = 0.4
    plot_adjustment_arrows = no
    max_mixture_model_iterations = 300
    num_mixture_model_initializations = 10
    extra_paralogs_analyses_methods = no
    max_mixture_model_components = 5
    max_mixture_model_ks = 5
    max_gene_family_size = 200

* **logging_level**: the lowest logging/verbosity level of messages printed to the console/logs (increasing severity: *notset*, *debug*, *info*, *warning*, *error*, *critical*). Messages less severe than *level* will be ignored; *notset* causes all messages to be processed. [Default: "info"]
* **distribution_peak_estimate**: the statistical method used to obtain a single ortholog *K*:sub:`S` estimate for the divergence time of a species pair from its ortholog distribution or to obtain a single paralog *K*:sub:`S` estimate from an anchor *K*:sub:`S` cluster or from lognormal components in mixture models (options: mode or median). [Default: "mode"]
* **kde_bandwidth_modifier**: modifier to adjust the fitting of the KDE curve on the underlying whole-paranome or anchor *K*:sub:`S` distribution. The KDE Scott's factor internally computed by SciPy tends to produce an overly smooth KDE curve, especially with steep WGD peaks, and therefore it is reduced by multiplying it by a modifier. Decreasing the modifier leads to tighter fits, increasing it leads to smoother fits, and setting it to 1 gives the default KDE factor. Note that a too small factor is likely to take into account data noise. [Default: 0.4]
* **plot_adjustment_arrows**: flag to toggle the plotting of rate-adjustment arrows below the adjusted mixed paralog--ortholog *K*:sub:`S` plot. These arrows start from the original unadjusted ortholog divergence *K*:sub:`S` estimate and end on the rate-adjusted estimate (options: "yes" and "no"). [Default: "yes"]
* **max_mixture_model_iterations**: maximum number of EM iterations for mixture modeling. [Default: 300]
* **num_mixture_model_initializations**: number of times the EM algorithm is initialized (either for the random initialization in the exponential-lognormal mixture model or for k-means in the lognormal mixture model). [Default: 10]
* **max_mixture_model_components**: maximum number of components considered during execution of the mixture models. [Default: 5]
* **max_mixture_model_ks**: upper limit for the Ks range in which the exponential-lognormal and lognormal-only mixture models are performed. [Default: 5]
* **max_gene_family_size**: maximum number of members that any paralog gene family can have to be included in *K*:sub:`S` estimation. Large gene families increase the run time and are often composed of unrelated sequences grouped together by shared protein domains or repetitive sequences. But this is not always the case, so one may want to check manually the gene families in file ``paralog_distributions/wgd_<focal species>/<focal species>.mcl.tsv`` and increase (or even decrease) this number. [Default: 200]
