Input and output description
****************************

Input files
===========

* FASTA files of coding regions (CDS) or of coding transcripts for each species in the dataset.
* GFF3 file for the focal species, allows collinearity (syntenic) analysis and anchor pair paralog *K*:sub:`S` distribution.
* *ksrates* configuration file(s) (for more details see :ref:`config_sections`).


.. _`output_files`:

Output files and directory organization
=======================================

* ``rate_adjustment/species`` directory collects the output files of the substitution rate-adjustment relative to the focal species.

    Figures:

        * Rate-adjusted mixed paralog--ortholog *K*:sub:`S` distribution plot in PDF format (``mixed_species_adjusted.pdf``).
        * Input tree with branch length set to *K*:sub:`S` distances estimated from ortholog *K*:sub:`S` distributions (``tree_species_distances.pdf``).
        * Multi-panel figure(s) of the ortholog *K*:sub:`S` distributions used to adjust a divergent species pair (``orthologs_species1_species2.pdf``).
        * Rate-adjusted mixed anchor pair--ortholog *K*:sub:`S` distribution clustered for inference of putative WGDs, with only significant clusters retained (``mixed_species_anchor_clusters.pdf``).
        * Rate-adjusted mixed paralog--ortholog *K*:sub:`S` distribution with superimposed exponential-lognormal mixture model inference of putative WGDs (``mixed_species_elmm.pdf``).
        * Rate-adjusted mixed paralog-- or anchor pair--ortholog *K*:sub:`S` distribution with superimposed lognormal-only mixture model for inference of putative WGDs (``mixed_species_lmm_colinearity.pdf`` and ``mixed_species_lmm_paranome.pdf``).
        * Unadjusted naive mixed paralog--ortholog *K*:sub:`S` distribution plot in PDF format (``mixed_species_unadjusted.pdf``).
        * Original input phylogenetic tree in PDF format with fixed branch length (``tree_species.pdf``)

    Files:

        * Rate-adjustment results in tab-separated format: raw results for each trio (``adjustment_table_species_all.tsv``) and final results for each divergent pair after finding a consensus value in case of multiple outgroups (``adjustment_table_species.tsv``).
        * Original input phylogenetic tree in ASCII format and list of sister species and outgroup species per node (``tree_species.txt``).
        * List of trios used for substitution rate-adjustment (``ortholog_trios_species.tsv``).
        * List of species pairs to be submitted to *wgd* ortholog runs (``ortholog_pairs_species.txt``).


* ``rate_adjustment/species/paralogs_analyses`` directory collects secondary output files produced during the inference of putative WGD signals through mixture modeling (see also section :ref:`paralogs_analyses`).

    From anchor *K*:sub:`S` clustering:

        * Anchor pair *K*:sub:`S` distribution with highlighted clusters of segment pair medians (``anchor_clusters_species_medians.pdf``).
        * Rate-adjusted mixed anchor pair--ortholog *K*:sub:`S` distributions clustered for inference of putative WGDs, with all inferred clusters (``mixed_species_anchor_clusters_unfiltered.pdf``).

    From exponential-lognormal mixture modeling:
    
        * Plots showing the kernel density estimation (KDE) and spline obtained from the log-transformed whole-paranome *K*:sub:`S` distribution (``elmm_species_kde_spline.pdf``).
        * Plots showing the peaks detected in the spline (``elmm_species_peaks.pdf``).
        * Multi-panel figure showing fitted mixture models obtained with data-driven and hybrid initializations (``elmm_species_models_data_driven.pdf``).
        * Multi-panel figure showing the best-fitted mixture model obtained for each number of components with random initialization (``elmm_species_models_random.pdf``).
        * TSV and TXT files collecting component parameters (``elmm_species_parameters.tsv`` and ``elmm_species_parameters.txt``) (more details on the file format in section :ref:`elmm`).

    From lognormal-only mixture modeling:

        * Multi-panel figure showing the best-fitted mixture model on whole-paranome and anchor pair *K*:sub:`S` distributions obtained for each number of components (``lmm_species_all_models_paranome.pdf`` and ``lmm_species_all_models_colinearity.pdf``).
        * TSV and TXT files collecting component parameters (``lmm_species_parameters_colinearity.tsv``, ``lmm_species_parameters_colinearity.txt``, ``lmm_species_parameters_paranome.tsv`` and ``lmm_species_parameters_paranome.txt``) (more details on the file format in section :ref:`lmm`).


* ``rate_adjustment/species/log_XXXXXXXX`` directory: when launching *ksrates* as a Nextflow pipeline, each execution generates a log directory named with a unique 8-character ID stated at the beginning of a Nextflow run. Details about how the processes of the workflow are proceeding and about encountered warnings or errors are stored in log files collected in this directory:

    * ``setup_adjustment.log`` shows the progress in checking input files and setting up species trios and pairs for rate-adjustment. 
    * ``wgd_paralogs.log`` shows the progress in estimating paralog *K*:sub:`S` values.
    * ``set_orthologs.log`` states whether ortholog *K*:sub:`S` data are already available or are missing for each species pair.
    * ``estimate_peak.log`` shows the progress in updating the ortholog *K*:sub:`S` databases from already existing ortholog *K*:sub:`S` data.
    * ``wgd_orthologs_species1_species2.log`` shows the progress in estimating ortholog *K*:sub:`S` values for a species pair.
    * ``plot_ortholog_distributions.log`` shows the progress in plotting the ortholog *K*:sub:`S` distributions.
    * ``rate_adjustment.log`` shows the progress in performing the actual rate-adjustment step.
    * ``paralogs_analyses.log`` shows the progress in analyzing the paralog distribution to detect potential WGD signatures through anchor *K*:sub:`S` clustering, exponential-lognormal mixture modeling and/or lognormal-only mixture modeling. 


* ``paralog_distributions/wgd_species`` directory contains the files generated during the paralog *K*:sub:`S` estimate for the focal species:

    * ``species.blast.tsv`` lists the paralog BLAST homology hits.
    * ``species.mcl.tsv`` lists the paralog gene families, one family per line from the largest to the smallest.
    * ``species.ks.tsv`` and  ``species.ks_anchors.tsv`` are tabular format files listing paralog or anchor pair hits (column 1) together with their *K*:sub:`S` estimate (column 9). Other pieces of information include alignment coverage, identity and length (columns 2 to 5) and gene family, tree node and weight (column 7, 10 and last column). For more details, see `wgd documentation <https://wgd.readthedocs.io/en/latest/methods.html?highlight=some%20information>`__.

    .. figure:: _images/ks_tsv.png
        :align: center
        :width: 800

    * ``species_i-adhore`` directory contains i-ADHoRe output files used during anchor *K*:sub:`S` clustering (see section :ref:`anchor_ks_clustering`)


* ``ortholog_distributions/wgd_species1_species2`` directory contains the files generated during the one-to-one ortholog *K*:sub:`S` estimate of a species pair:

    * ``species1_species2.blast.tsv`` lists the ortholog BLAST homology hits. When the execution of the *wgd* ortholog run is over it is possible to delete this file to save disk space.
    * ``species1_species2.orthologs.tsv`` lists the one-to-one ortholog reciprocal best hits between the two species, one hit per line.
    * ``species1_species2.ks.tsv`` lists the one-to-one ortholog reciprocal best hits (column 1) together with their *K*:sub:`S` estimate (column 9). The tabular file format is identical to the paralog ``.ks.tsv`` file described above. However, gene family, tree node and weight columns are of less interest for orthologs since each family is composed of only two members.


* Generated directly in the launching directory:

    * Databases storing the ortholog *K*:sub:`S` lists (``ks_list_database_path.txt``) and the estimated divergence time *K*:sub:`S` estimate (``peak_database_path.txt``) of the ortholog *K*:sub:`S` distributions. Their location can be customised in the configuration file.
    * List of commands to launch the ortholog *wgd* runs in the manual pipeline (``wgd_runs_species.txt``). Note that this file is not generated if launching the Nextflow pipeline.
    * The ``work`` directory is automatically generated by Nextflow to handle process organization and communication between processes (for more details see Nextflow documentation, e.g. the Get started `page <https://www.nextflow.io/docs/latest/getstarted.html#your-first-script>`__).


Note on *wgd* output files
==========================

If a *ksrates* Nextflow pipeline run is prematurely interrupted for some reasons (e.g. cancelled by the user or crashed) while one or more *wgd* runs were still ongoing, the latter will leave temporary directories and incomplete files within ``paralog_distributions`` and/or ``ortholog_distributions``. Such leftovers must be manually removed before relaunching the Nextflow pipeline to avoid that the next run continues the task from incomplete data. For safety, if the pipeline encounters some leftovers it will immediately stop and return an error message in the Nextflow log files (``wgd_paralogs.log`` and/or ``wgd_orthologs_species1_species2.log``).
