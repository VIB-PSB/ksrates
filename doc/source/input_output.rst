Input and output description
****************************

Input files
===========

* FASTA files of coding regions (CDS) or of coding transcripts for each species in the dataset.
* GFF3 file for the focal species, allows collinearity (syntenic) analysis and anchor pair paralog *K*:sub:`S` distribution.
* *ksrates* configuration file(s) (for more details see :ref:`config_sections`).


Output folders and folder organization
======================================

*ksrates* generates the following output folders and files within the folder where it is executed from:

* ``rate_adjustment/species``: contains the main and secondary output files of the substitution rate-adjustment (figures and data files) for the focal species, see next section.
* ``paralog_distributions/wgd_species``: contains the files of the ``wgd`` paralog *K*:sub:`S` estimation for the focal species (``species.ks.tsv``)
* ``ortholog_distributions/wgd_species1_species2``: contains the files of the ``wgd`` one-to-one ortholog *K*:sub:`S` estimation of a species pair (``species1_species2.ks.tsv``)
* ``rate_adjustment/species/log_XXXXXXXX``: a folder collecting all log files produced by processes of the *ksrates* Nextflow pipeline. The log files collect standard and error output of the commands being executed. Each Nextflow run produces a new log folder named with a unique 8-character ID and stated at beginning of the Nextflow run.


Output files
============

Main output files:
------------------

.. TODO: rename "corrected" and "rates" in output file names

* Rate-adjusted mixed paralog--ortholog *K*:sub:`S` distribution plot in PDF format (``mixed_species_adjusted.pdf``)
* Rate-adjusted mixed anchor pair--ortholog *K*:sub:`S` distribution clustered to infer putative WGD components (``mixed_species_anchor_clusters.pdf``)
* Rate-adjusted mixed paralog--ortholog *K*:sub:`S` distribution with superimposed exponential-lognormal mixture model (``mixed_species_elmm.pdf``)
* Rate-adjusted mixed paralog-- or anchor pair--ortholog *K*:sub:`S` distribution with superimposed lognormal-only mixture model (``mixed_species_lmm_colinearity.pdf`` and ``mixed_species_lmm_paranome.pdf``)
* Rate-adjustment results in tab-separated format: raw results for each trio (``adjustment_table_species_all.tsv``) and final results for each divergent pair after finding a consensus value in case of multiple outgroups (``adjustment_table_species.tsv``)
* Input tree with branch length set to *K*:sub:`S` distances estimated from ortholog *K*:sub:`S` distributions (``tree_species_distances.pdf``)

Secondary output files:
-----------------------

* Original input phylogenetic tree in PDF format (``tree_species.pdf``)
* Original input phylogenetic tree in ASCII format and list of sister species and outgroup species per node (``tree_species.txt``)
* List of trios used for substitution rate-adjustment (``ortholog_trios_species.tsv``)
* List of species pairs to be submitted to ``wgd`` ortholog pipeline (``ortholog_pairs_species.txt``)
* Databases storing the *K*:sub:`S` list and the estimated divergence time Ks estimate of each ortholog distribution needed for the rate-adjustment
* Un-adjusted naive mixed paralog--ortholog *K*:sub:`S` distribution plot in PDF format (``mixed_species_unadjusted.pdf``)
* Multi-panel figure(s) of the ortholog *K*:sub:`S` distributions used to adjust a divergent species pair (``orthologs_species1_species2.pdf``)
* Rate-adjusted mixed anchor pair paralog--ortholog *K*:sub:`S` distribution clustered to infer putative WGD components, with *all* inferred clusters (``mixed_species_anchor_clusters_unfiltered.pdf``)

Nextflow log files:
-------------------

* When running *ksrates* as a Nextflow pipeline, each step of the pipeline produces a log file that collects output and error messages from the executed commands. Details about how the pipeline proceeds and about errors are stored there, instead of being printed on screen. See below.


Note on wgd output files
========================

If a *ksrates* Nextflow pipeline run is prematurely interrupted for some reasons (e.g. cancelled by the user or crashed) while the ``wgd`` pipelines were still running, the latter will leave temporary folders and incomplete files. Such leftovers must be manually removed before starting a new Nextflow run to avoid that the next run continues the task from incomplete data. For safety, if the pipeline encounters some leftovers it will immediately stop and return an error in the Nextflow log files.

Tip to save disk space: when the execution of the ``wgd`` ortholog runs is over it is possible to delete the ``.blast.tsv`` file in their ortholog distribution folder, since it is of no use anymore and can take up quite some space.
