import sys
import os
import pandas
import logging
from pandas import DataFrame
from ast import literal_eval
import ksrates.fc_check_input as fcCheck
import ksrates.fc_kde_bootstrap as fcPeak
import ksrates.fc_configfile as fcConf
from ksrates.utils import init_logging


def compute_peaks(config_file, expert_config_file, ortholog_pairs_file):
    # INPUT
    config = fcConf.Configuration(config_file, expert_config_file)
    init_logging("Computing ortholog distribution peaks with related error", config.get_logging_level())
    logging.info("Loading parameters and input files")

    # Get parameters and input files
    species = config.get_species()
    tree = config.get_newick_tree()
    latin_names = config.get_latin_names()
    max_ks_ortho = config.get_max_ks_ortho()
    n_iter = config.get_num_iteration()
    bin_width_ortho = config.get_bin_width_ortho()
    x_lim_ortho = config.get_x_lim_ortho()

    default_path_ortholog_pairs_file = os.path.join("rate_adjustment", f"{species}", f"ortholog_pairs_{species}.tsv")
    ortholog_pairs_file = fcCheck.get_argument_path(ortholog_pairs_file, default_path_ortholog_pairs_file,  "Ortholog pairs file")
    if ortholog_pairs_file == "":
        logging.error(f"Ortholog pairs file not found at default position [{default_path_ortholog_pairs_file}].")
        logging.error("Exiting")
        sys.exit(1)
    with open(ortholog_pairs_file, "r") as f:
        header = f.readline() # ignore first line due to headers
        ortholog_pairs = f.readlines()

    db_path = config.get_ortho_db()
    ks_list_db_path = config.get_ks_db()

    try:
        with open(db_path, "r") as f:
            db = pandas.read_csv(f, sep="\t", index_col=0)
    except Exception:
        logging.warning(f"Ortholog peak database [{db_path}] not found or empty: a new one is now generated.")

        db = DataFrame()
        with open(db_path, "w+") as outfile:
            outfile.write('\tSpecies1\tSpecies2\tMode\tMode_SD\n')

    try:
        with open(ks_list_db_path, "r") as f:
            ks_list_db = pandas.read_csv(f, sep="\t", index_col=0)
            # When imported from csv format, all the Ks lists in the df are read as
            # plain text (strings) and must be converted back to value lists
            ks_list_db.loc[:, 'Ks_Values'] = ks_list_db.loc[:, 'Ks_Values'].apply(literal_eval)
    except Exception:
        logging.warning(f"Ortholog Ks list database [{ks_list_db_path}] not found or empty: a new one is now generated.")

        ks_list_db = DataFrame()
        with open(ks_list_db_path, "w+") as outfile:
            outfile.write('\tSpecies1\tSpecies2\tKs_Values\n')

    logging.info("")

    # -----------------------------------------------------------------------------

    if len(ortholog_pairs) == 0:
        logging.info(f"There are no ortholog species pairs listed in [{ortholog_pairs_file}].") 
        logging.info(f"Nothing to do.")
        sys.exit(0) # exit code 0 because no actual errors were thrown (the pair list
                    # is empty if such ortholog data were already obtained in previous runs)

    # Checking if the ortholog pair files contains correct format and valid species names
    format_error, species_name_error = False, False
    for pair_string in ortholog_pairs:
        # File parsing support both tab space (inserted by the script which generates the file, init.py),
        # but also accepts a single space in case the user manually write in the file.
        if len(pair_string.rstrip().split("\t")) == 2:
            pair = pair_string.rstrip().split("\t")
        elif len(pair_string.rstrip().split(" "))==2:
            pair = pair_string.rstrip().split(" ")
        else:
            format_error = True
            continue
        # good format, now let's check if the species names are valid
        if pair[0] not in tree.get_leaf_names() or pair[1] not in tree.get_leaf_names():
            species_name_error = True

    if format_error:
        logging.error(f"Format error in [{ortholog_pairs_file}]. Please check that each line contains two species names separated by a tabular space or by a single space.")
    if species_name_error:
        logging.error(f"One or more invalid (misspelled?) species names are listed in [{ortholog_pairs_file}]")
    if format_error or species_name_error:
        logging.error(f"Exiting")
        sys.exit(1)


    # CHECKING FOR WHICH SPECIES PAIRS THE ORTHOLOG PEAK IS MISSING OR THE KS LIST IS MISSING

    list_of_compute_peak_results = [] # will be the flag for unsuccessful peak computations
    for pair_string in ortholog_pairs:
        if len(pair_string.rstrip().split("\t")) == 2: # e.g. [ ['sp1', 'sp2'], ['sp3','sp4'] ]
            pair = pair_string.rstrip().split("\t")
        elif len(pair_string.rstrip().split(" "))==2:
            pair = pair_string.rstrip().split(" ")
        pair.sort(key=str.lower)
        sp1, sp2 = pair[0], pair[1]
        latin_name_list = sorted([latin_names[sp1], latin_names[sp2]], key=str.casefold)
        latinSp1, latinSp2 = latin_name_list[0], latin_name_list[1]

        # Checking if there are temporary folders in the wgd folder of the two species (in case, their peak in the DB or the Ks list that will be extracted could be incomplete) 
        if os.path.isdir(os.path.join("ortholog_distributions", f"wgd_{sp1}_{sp2}", f"{sp1}_{sp2}.ks_tmp")) or os.path.isdir(os.path.join("ortholog_distributions", f"wgd_{sp1}_{sp2}", f"{sp1}_{sp2}.blast_tmp")):
            logging.warning(f"One or more temporary folders have been found in [ortholog_distributions/wgd_{sp1}_{sp2}]: the ortholog data may be incomplete!")
            logging.warning("It is advised to delete the temporary folders and the file that was being created (.blast.tsv or .ks.tsv) and re-compute the data.")
            logging.warning("")

        if f"{latinSp1}_{latinSp2}" not in db.index: # flag for presence/absence of species pair in the ortholog peak database
            flag_not_in_peak_db = True 
        else:
            flag_not_in_peak_db = False

        if f"{latinSp1}_{latinSp2}" not in ks_list_db.index: # same but for ks list database
            flag_not_in_ks_db = True
        else:
            flag_not_in_ks_db = False

        if flag_not_in_peak_db is False and flag_not_in_ks_db is False:
            logging.info(f"{latinSp1} and {latinSp2}: data already present in ortholog database")
            
        if flag_not_in_peak_db or flag_not_in_ks_db:    # if the pair is missing in the ortholog peak or ks databases
            compute_peak_failed_flag = fcPeak.estimate_peak(sp1, sp2, latinSp1, latinSp2, max_ks_ortho, n_iter, x_lim_ortho, bin_width_ortho, ks_list_db_path, db_path, flag_not_in_peak_db=flag_not_in_peak_db, flag_not_in_ks_db=flag_not_in_ks_db)
            list_of_compute_peak_results.append(compute_peak_failed_flag)
        logging.info("")

    if True in list_of_compute_peak_results:
        logging.warning(f"Number of failed peak computations: {list_of_compute_peak_results.count(True)}")
        logging.warning("")
    logging.info("All done")
