import sys
import os
import logging
import ksrates.fc_manipulate_trees as fcTree
import ksrates.fc_check_input as fcCheck
import ksrates.fc_configfile as fcConf
from ksrates.fc_rrt_correction import _ADJUSTMENT_TABLE
from ksrates.utils import init_logging
import pandas


def plot_tree_rates(config_file, correction_table_file, nextflow_flag):
    # INPUT
    config = fcConf.Configuration(config_file)
    init_logging("Generating PDF of input tree with branch length equal to Ks distances", config.get_logging_level())
    logging.info("Loading parameters and input files")

    # GET PARAMETERS and INPUT FILES
    species = config.get_species()
    peak_db_path = config.get_ortho_db()
    newick_tree = config.get_newick_tree() # as Tree object by ete3
    latin_names = config.get_latin_names(newick_tree)

    # Get correction results TSV file
    default_path_correction_table_file = os.path.join("rate_adjustment", f"{species}", f"{_ADJUSTMENT_TABLE.format(species)}")
    correction_table_file = fcCheck.get_argument_path(correction_table_file, default_path_correction_table_file, "Correction table file")
    if correction_table_file == "": # it means that the correction_table is not present or available yet
        logging.warning(f"Rate-adjustment data not available yet: PDF figure of phylogenetic tree not generated.")
        logging.info(f"Exiting")
        sys.exit(1) # exit 1 because plot_tree is executed at the end of the Nextflow pipeline and correction_table should exits
    else:
        with open(correction_table_file, "r") as f:
            correction_table = pandas.read_csv(f, sep="\t")
            species_in_correction_table = fcTree.counts_expected_line_number_in_correction_table(species, newick_tree, latin_names)
            missing_required_rates = set(sorted(species_in_correction_table)) - set(sorted(correction_table["Sister_Species"]))

            if correction_table.shape[0] == 0:
                logging.warning(f"Rate-adjustment data not available yet: PDF figure of phylogenetic tree not generated.")
                logging.info(f"Exiting")
                sys.exit(1) # exit 1 because plot_tree is executed at the end of the Nextflow pipeline and correction_table should be completed
            elif len(missing_required_rates) != 0:
                # Having a complete correction table is strictly required for building the tree,
                # because all the branch-specific Ks contributions contained in there are needed.
                # Other extra-table branch contributions may be additionally required to fill in all the branch lengths,
                # but their absence is tolerated.
                logging.warning(f"The Ks distances between {latin_names[species]} and the following species are missing from the rate-adjustment table. Please compute them before building the tree:")
                for name in sorted(missing_required_rates):
                    logging.warning(f" - {name}")
                logging.warning(f"Exiting")
                sys.exit(1) # exit 1 because plot_tree is executed at the end of the Nextflow pipeline and correction_table should be completed

    # Getting the choice on how to deal with the presence of multiple corrections for the same divergent pair
    # due to the use of multiple trios/outgroup during correction
    # Available options:
    #  - 'mean among outgroups': taking the average of the corrected peaks
    #  - 'best outgroup': taking the corrected peak coming from the best outgroup, which is the one with smallest OC segment
    consensus_peak_for_multiple_outgroups = config.get_consensus_peak_for_multiple_outgroups()
    # Getting the statistical measure for how to determine the representative value of an ortholog distribution
    peak_stats = config.get_peak_stats() # default is mode (other option, median)

    # Loading the ortholog peak database for the tree picture: it can contain data to compute all the branch lengths.
    try:
        with open(peak_db_path, "r") as f:
            ortholog_db = pandas.read_csv(f, sep="\t", index_col=0)
            if ortholog_db.shape[0] == 0:
                logging.warning(f"Ortholog peak database is present by doesn't contain any data, will be ignored when plotting the tree")
                ortholog_db = pandas.DataFrame()
    except Exception:
        logging.warning(f"Ortholog peak database [{peak_db_path}] not found or empty, will be ignored when plotting the tree")
        ortholog_db = pandas.DataFrame() # since there is no ortholog database, just assign an empty dataframe to this variable

    logging.info("")
    fcTree.plotting_tree(species, latin_names, newick_tree, correction_table, consensus_peak_for_multiple_outgroups, ortholog_db, peak_stats, nextflow_flag)
    logging.info("")
    logging.info(f"Saved PDF tree figure [{fcTree._TREE_BRANCH_DISTANCES.format(species)}]")
    logging.info("")
    logging.info("All done")