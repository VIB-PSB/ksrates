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

.. note::
    In the following listings of directory and file names, ``species`` is used as a placeholder for the actual (informal) name of the focal species (e.g. ``elaeis``) as specified in the *ksrates* configuration file.


Main output
-----------

* ``rate_adjustment/species``: this directory collects the output files of the substitution rate-adjustment relative to the focal species.

    Figures:

        * Rate-adjusted mixed paralog--ortholog *K*:sub:`S` distribution plot in PDF format (``mixed_species_adjusted.pdf``).
        * Input phylogenetic tree in PDF format with branch length set to *K*:sub:`S` distances estimated from ortholog *K*:sub:`S` distributions (``tree_species_distances.pdf``).
        * Rate-adjusted mixed anchor pair--ortholog *K*:sub:`S` distribution clustered for inference of putative WGDs, with only significant clusters retained (``mixed_species_anchor_clusters.pdf``).
        * Rate-adjusted mixed paralog--ortholog *K*:sub:`S` distribution with superimposed exponential-lognormal mixture model inference of putative WGDs (``mixed_species_elmm.pdf``).
        * Rate-adjusted mixed paralog-- and anchor pair--ortholog *K*:sub:`S` distributions with superimposed lognormal-only mixture model for inference of putative WGDs (``mixed_species_lmm_paranome.pdf`` and ``mixed_species_lmm_colinearity.pdf``).
        * Multi-panel figure(s) of the ortholog *K*:sub:`S` distributions used to adjust a divergent species pair (``orthologs_species1_species2.pdf``).
        * Unadjusted naive mixed paralog--ortholog *K*:sub:`S` distribution plot in PDF format (``mixed_species_unadjusted.pdf``).
        * Original input phylogenetic tree in PDF format with fixed branch lengths (``tree_species.pdf``).

    Files:

        *   Raw rate-adjustment results for each trio (``adjustment_table_species_all.tsv``). Tabular format.

            Each row shows the result for a species pair (column 2 ``Focal_Species`` and 3 ``Sister_Species``) diverging at a certain node (column 1 ``Node``) and adjusted with the outgroup in column 3 ``Out_Species``. The rate-adjusted mode with associated standard deviation are given in column 4 ``Adjusted_Mode`` and 5 ``Adjusted_Mode_SD``; for comparison the unadjusted original mode with associated standard deviation is provided in column 6 ``Original_Mode`` and 7 ``Original_Mode_SD``. The branch-specific *K*:sub:`S` contributions for the divergent species pair are listed in column 8 ``Ks_Focal`` and 9 ``Ks_Sister``; the *K*:sub:`S`distance of the outgroup to the divergence event of the species pair is listed in column 10 ``Ks_Out``.

            .. figure:: _images/adj_table_all.png
                :align: center
                :width: 800

                Raw rate-adjustment results on a divergent pair using four outgroups.

        *   Final rate-adjustment results for each divergent species pair after finding a consensus value in case of multiple outgroups (``adjustment_table_species.tsv``). Tabular format.
        
            Each row shows the result for a species pair (column 2 ``Focal_Species`` and 3 ``Sister_Species``) diverging at a certain node (column 1 ``Node``). Columns 4--7 report the consensus obtained by taking the *mean* of multiple outgroups (if available): rate-adjusted mode with standard deviation in column 4 ``Adjusted_Mode_Mean`` and 5 ``Adjusted_Mode_SD_Mean``, branch-specific *K*:sub:`S` contributions for the divergent species pair in column 6 ``Ks_Focal_Mean`` and 7 ``Ks_Sister_Mean``. Columns 8-11 report the consensus obtained when considering only the *best outgroup*: rate-adjusted mode with standard deviation in column 8 ``Adjusted_Mode_Best`` and 9 ``Adjusted_Mode_SD_Best``, *K*:sub:`S` contributions for the divergent species pair in column 10 ``Ks_Focal_Best`` and 11 ``Ks_Sister_Best``. For comparison the unadjusted original mode with associated standard deviation is provided in column 12 ``Original_Mode`` and 13 ``Original_Mode_SD``.

            .. figure:: _images/adj_table_consensus.png
                :align: center
                :width: 800

                Consensus result for the divergent pair obtained from the four raw rate-adjustments.
        
        * Original input phylogenetic tree in ASCII format and list of sister species and outgroup species per node (``tree_species.txt``).
        * List of trios used for substitution rate-adjustment (``ortholog_trios_species.tsv``).
        * List of species pairs for which ortholog *K*:sub:`S` distributions are estimated using *wgd* (``ortholog_pairs_species.txt``).


* ``rate_adjustment/species/paralogs_analyses``: this directory collects secondary output files produced during the inference of putative WGD signals through mixture modeling (see :ref:`paralogs_analyses`).

    From anchor *K*:sub:`S` clustering:

        * Anchor pair *K*:sub:`S` distribution with highlighted clusters of segment pair medians (``anchor_clusters_species_medians.pdf``).
        * Rate-adjusted mixed anchor pair--ortholog *K*:sub:`S` distributions clustered for inference of putative WGDs, with all inferred clusters (``mixed_species_anchor_clusters_unfiltered.pdf``).

    From exponential-lognormal mixture modeling:
    
        * Plots showing the kernel density estimation (KDE) and spline obtained from the log-transformed whole-paranome *K*:sub:`S` distribution (``elmm_species_kde_spline.pdf``).
        * Plots showing the peaks detected in the spline (``elmm_species_peaks.pdf``).
        * Multi-panel figure showing fitted mixture models obtained with data-driven and hybrid initializations (``elmm_species_models_data_driven.pdf``).
        * Multi-panel figure showing the best-fitted mixture model obtained for each number of components with random initialization (``elmm_species_models_random.pdf``).
        * TSV and TXT files collecting component parameters (``elmm_species_parameters.tsv`` and ``elmm_species_parameters.txt``) (see :ref:`elmm` for more details on the file format).

    From lognormal-only mixture modeling:

        * Multi-panel figure showing the best-fitted mixture model on whole-paranome and anchor pair *K*:sub:`S` distributions obtained for each number of components (``lmm_species_all_models_paranome.pdf`` and ``lmm_species_all_models_colinearity.pdf``).
        * TSV and TXT files collecting component parameters (``lmm_species_parameters_colinearity.tsv``, ``lmm_species_parameters_colinearity.txt``, ``lmm_species_parameters_paranome.tsv`` and ``lmm_species_parameters_paranome.txt``) (see :ref:`lmm` for more details on the file format).


Nextflow log files
------------------

* ``rate_adjustment/species/log_XXXXXXXX``: when launching *ksrates* as a Nextflow pipeline, each execution generates a log directory named with a unique 8-character ID stated at the beginning of a Nextflow run. Details about how the processes of the workflow are proceeding and about encountered warnings or errors are stored in log files collected in this directory:

    * ``setup_adjustment.log`` shows the progress in checking input files and setting up species trios and pairs for rate-adjustment. 
    * ``wgd_paralogs.log`` shows the progress in estimating paralog *K*:sub:`S` values.
    * ``set_orthologs.log`` states whether ortholog *K*:sub:`S` data are already available or are missing for each species pair.
    * ``estimate_peak.log`` shows the progress in updating the ortholog *K*:sub:`S` databases from already existing ortholog *K*:sub:`S` data.
    * ``wgd_orthologs_species1_species2.log`` shows the progress in estimating ortholog *K*:sub:`S` values for a species pair.
    * ``plot_ortholog_distributions.log`` shows the progress in plotting the ortholog *K*:sub:`S` distributions.
    * ``rate_adjustment.log`` shows the progress in performing the actual rate-adjustment step.
    * ``paralogs_analyses.log`` shows the progress in analyzing the paralog distribution to detect potential WGD signatures through anchor *K*:sub:`S` clustering, exponential-lognormal mixture modeling and/or lognormal-only mixture modeling. 


*K*:sub:`S` estimate output (*wgd*)
-----------------------------------

* ``paralog_distributions/wgd_species``: this directory contains the files generated during the *wgd* paralog *K*:sub:`S` estimation run for the focal species:

    * ``species.blast.tsv`` lists the paralog BLAST homology hits in tabular output format (``-outfmt 6``) 
    * ``species.mcl.tsv`` lists the paralog gene families, one family per line from the largest to the smallest family with the gene IDs of individual family members separated by tabs.
    *   ``species.ks.tsv`` and  ``species.ks_anchors.tsv`` are tabular format files listing the *K*:sub:`S` estimate (column 9 ``Ks``) for every paralog and anchor pair found, respectively. Other noteworthy data per pair includes the alignment coverage, identity and length (columns 2 to 5: ``AlignmentCoverage``, ``AlignmentIdentity``, ``AlignmentLength`` and ``AlignmentLengthStripped``), the gene family (column 7 ``Family``), the node in the gene family's tree (column 10 ``Node``), and the weight associated with the pair's *K*:sub:`S` estimate (column 15 ``WeightOutliersExcluded``). For more details, see the *wgd* `documentation <https://wgd.readthedocs.io/en/latest/methods.html?highlight=some%20information>`__.

        .. figure:: _images/ks_tsv.png
            :align: center
            :width: 800

            File section showing the structure of the ``.ks.tsv`` format.

    * ``species_i-adhore``: this directory contains the i-ADHoRe output files necessary for the anchor *K*:sub:`S` clustering (see :ref:`anchor_ks_clustering`).


* ``ortholog_distributions/wgd_species1_species2``: these directories contain the files generated during the *wgd* one-to-one ortholog *K*:sub:`S` estimation for each species pair:

    * ``species1_species2.blast.tsv`` lists the ortholog BLAST homology hits.
    
      .. note::
          When the *wgd* ortholog *K*:sub:`S` estimation analysis is finished it is possible to delete this file to save disk space.
        
    * ``species1_species2.orthologs.tsv`` lists the one-to-one ortholog (i.e. the reciprocal best BLAST hits) between the two species, one ortholog pair per line.
    * ``species1_species2.ks.tsv`` lists the *K*:sub:`S` estimate (column 9 ``Ks``) for every one-to-one ortholog pair found. The tabular file format is identical to the paralog ``.ks.tsv`` file described above. However, the gene family, tree node and weight columns can be ignored since each ortholog "family" is composed of only two members.


Other output
------------

* Generated directly in the directory from where *ksrates* is launched:

    * ``ortholog_peak_db.tsv`` is a tabular data file storing the *K*:sub:`S` mode estimate from the ortholog *K*:sub:`S` distribution of species pairs. The name and location can be customised in the *ksrates* configuration file.
    * ``ortholog_ks_list_db.tsv`` is a tabular data file storing the ortholog *K*:sub:`S` value lists of species pairs.  The name and location can be customised in the *ksrates* configuration file.
    * ``wgd_runs_species.txt`` contains a list of *ksrates* commands to launch the *wgd* paralog and ortholog analysis when using the manual pipeline (see :ref:`manual_pipeline`). Note that this file is not generated if using the *ksrates* Nextflow pipeline.
    * ``work``: when using the *ksrates* Nextflow pipeline this directory is automatically generated by Nextflow to handle process organization and communication between processes (for more details, see the Nextflow documentation, e.g. `here <https://www.nextflow.io/docs/latest/getstarted.html#your-first-script>`__).


Note on *wgd* output files
==========================

If a *ksrates* Nextflow pipeline run is prematurely interrupted for some reasons (e.g. cancelled by the user or crashed) while one or more *wgd* runs were still ongoing, the latter will leave temporary directories and incomplete files within ``paralog_distributions`` and/or ``ortholog_distributions`` (e.g. BLAST files). Such leftovers are by default automatically detected and removed at the end of the workflow as a safety measure to avoid that the next run continues the task from incomplete data.

It is possible to preserve the leftover files for investigating what caused the pipeline to crash (see ``preserve`` parameter in :ref:`nextflow_config_section`). In this case it will be later necessary to manually remove the leftovers before relaunching the pipeline, otherwise the workflow will immediately stop and return an error message in the Nextflow log files (``wgd_paralogs.log`` and/or ``wgd_orthologs_species1_species2.log``).
