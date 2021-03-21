import sys
import os
import pandas as pd
from pandas import DataFrame
from math import sqrt, floor
import logging
import ksrates.fc_rrt_correction as fcCorrect
import ksrates.fc_check_input as fcCheck
import ksrates.fc_configfile as fcConf
from ksrates.utils import init_logging


def correct(config_file, trios_file):
    # INPUT
    config = fcConf.Configuration(config_file)
    init_logging("Rate-adjustment of ortholog Ks distributions", config.get_logging_level())
    logging.info("Loading parameters and input files")

    # Get parameters from configfile
    species_of_interest = config.get_species()
    latin_names = config.get_latin_names()
    db_path = config.get_ortho_db()

    default_path_trios_file = os.path.join("rate_adjustment", f"{species_of_interest}", f"ortholog_trios_{species_of_interest}.tsv")
    trios_file = fcCheck.get_argument_path(trios_file, default_path_trios_file, "Trios TSV file")
    if trios_file == "":
        logging.error(f"Trios TSV file not found at default position [{default_path_trios_file}].")
        logging.error("Exiting")
        sys.exit(1)

    try:
        with open(db_path, "r") as f:
            db = pd.read_csv(f, sep="\t", index_col=0)
    except Exception:
        logging.error(f"Ortholog peak database [{db_path}] not found or empty\n\
            -> correction will be skipped.")
        sys.exit(0)

    # Getting the statistical measure for how to determine the representative value of an ortholog distribution
    peak_stats = config.get_peak_stats() # default is mode (other option, median)

    # Getting the choice on how to deal with the presence of multiple corrections for the same divergent pair
    # due to the use of multiple trios/outgroup during correction
    consensus_peak_for_multiple_outgroups = config.get_consensus_peak_for_multiple_outgroups()

    # -------------------------------------------------------------------

    all_trios_correction_array = []  # will contain the corrected peak for each trio
    all_pairs_array = []    # will contain the corrected peak for each divergent pair after getting a consensus from multiple outspecies; 
                            # (both "best OC" and "multiple outspecies" strategies)
    sisters_per_node = {}   # keys=nodes, values=list with sisters

    # FILLING IN THE DATAFRAME FOR ALL TRIOS
    # It lists the correction results for each outgroup that has been used
    logging.info("")
    logging.info(f"Performing rate-adjustment of each divergent pair by using one or more outgroups:")
    with open(trios_file, "r") as f:
        trios = pd.read_csv(f, sep="\t")

    for __, row in trios.iterrows():
        node = row['Node']
        species, sister, out = row['Species'], row['Sister_Species'], row['Out_Species']
        latinSpecies, latinSister, latinOut = latin_names[species], latin_names[sister], latin_names[out]
        logging.info(f" - Adjusting the peak for [{latinSpecies}] and [{latinSister}] with outspecies [{latinOut}]")

        species_sister = "_".join(sorted([latinSpecies, latinSister], key=str.casefold))  # e.g. A.filiculoides_S.cucullata
        species_out = "_".join(sorted([latinSpecies, latinOut], key=str.casefold))
        sister_out = "_".join(sorted([latinSister, latinOut], key=str.casefold))

        if species_sister in db.index and species_out in db.index and sister_out in db.index:
            rate_species, rate_species_sd, rate_sister, rate_sister_sd = fcCorrect.decompose_ortholog_ks(db, species_sister, species_out, sister_out, peak_stats)
            correct_peak, correct_sd = fcCorrect.compute_corrected_ks_species_sister(rate_species, rate_species_sd)
            # OC_segment is a measure of better/worse outgroup choices for the decomposition into branch-specific Ks contributions; see documentation.
            OC_segment = db.loc[species_out]['Ortholog_Mode'] - rate_species

            all_trios_correction_array.append([node, latinSpecies, latinSister, latinOut, round(OC_segment, 6), round(correct_peak, 6), round(correct_sd, 6), round(rate_species, 6), round(rate_sister, 6)])

            if node not in sisters_per_node:
                sisters_per_node[node] = []
            if latinSister not in sisters_per_node[node]:
                sisters_per_node[node].append(latinSister)

        else:  # missing ortholog data
            logging.warning(f"Couldn't process trio [{latinSpecies}, {latinSister}, {latinOut}]:")
            if species_sister not in db.index: logging.warning(f" - [{species_sister}] not in ortholog peak database.")
            if species_out not in db.index: logging.warning(f" - [{species_out}] not in ortholog peak database.")
            if sister_out not in db.index: logging.warning(f" - [{sister_out}] not in ortholog peak database.")

    # Generating file with correction data for each trio.
    all_trios_correction_df = DataFrame.from_records(all_trios_correction_array, columns=["Node", "Species", "Sister_Species",
                                                        "Out_Species", "Ks_OC", "Mode", "SD", "Ks_Species", "Ks_Sister"])
    with open(os.path.join("rate_adjustment", f"{species_of_interest}", f"{fcCorrect._ADJUSTMENT_TABLE_ALL.format(species_of_interest)}"),
            "w+") as outfile:
        outfile.write(all_trios_correction_df.to_csv(sep="\t", index=False))
    logging.info(f"Rate-adjustment results for each trio saved in TSV format [{fcCorrect._ADJUSTMENT_TABLE_ALL.format(species_of_interest)}]")
    logging.info("")

    # FILLING IN THE DATAFRAME FOR ALL ORTHOLOG PAIRS WITH FOCAL SPECIES
    # (ALL DIVERGENCE EVENTS OF FOCAL SPECIES)
    if len(sisters_per_node) != 0:
        logging.info(f"Finding a consensus value in case multiple outgroups have been used to adjust a divergent pair [strategy: {consensus_peak_for_multiple_outgroups}]")

    # Get the headers for the file with data to plot the tree
    df_tested_outgroups_per_sister = all_trios_correction_df.loc[all_trios_correction_df['Node'] == 1]
    number_of_lines = len(df_tested_outgroups_per_sister.index)

    for node in sisters_per_node:
        # to select all the lines with "Node" field equal to 0 (or to 1, 2...)
        node_df = all_trios_correction_df.loc[all_trios_correction_df['Node'] == node]

        for sister in sisters_per_node[node]:
            # FIRST STRATEGY: taking the MEAN among corrected peaks from all outgroups
            # TODO: take the median too?
            peak_list = node_df.loc[node_df['Sister_Species'] == sister, ['Mode']]
            sd_list = node_df.loc[node_df['Sister_Species'] == sister, ['SD']]
            sd_list = sd_list["SD"].values.tolist()

            peak_mean = float(peak_list.mean())
            # Computing st.dev of the mean peak by following error propagation rules (sum of squares, square root; division)
            sd_err_prop = 0
            for sd in sd_list:
                sd_err_prop += pow(sd, 2)
            sd_err_prop = sqrt(sd_err_prop) / len(sd_list)
            # Getting the mean rate_species and the mean rate_sister out of the corrections (when multiple trios/outgroups are used)
            rate_species_list = node_df.loc[node_df['Sister_Species'] == sister, ['Ks_Species']]
            rate_species_mean = float(rate_species_list.mean())
            rate_sister_list = node_df.loc[node_df['Sister_Species'] == sister, ['Ks_Sister']]
            rate_sister_mean = float(rate_sister_list.mean())        

            # SECOND STRATEGY: taking only the corrected peak from the "BEST" outgroup (lowest OC value)
            # TODO: use DataFrame.at?
            oc_list = node_df.loc[node_df['Sister_Species'] == sister, ['Ks_OC']]
            oc_best_value = float(oc_list.min())
            peak_best_oc = node_df.loc[node_df['Ks_OC'] == oc_best_value, ['Mode']]
            peak_best_oc_float = float(peak_best_oc.mean()) # trick to make it a float number
            sd_best_oc = node_df.loc[node_df['Ks_OC'] == oc_best_value, ['SD']]
            sd_best_oc_float = float(sd_best_oc.mean())  # trick to make it a float number
            # Getting the rate_species and Ks_Sister associated to the correction with the best outgroup 
            rate_species_best_out = node_df.loc[node_df['Ks_OC'] == oc_best_value, ['Ks_Species']]
            rate_species_best_out_float = float(rate_species_best_out.mean()) # trick to make it a float number
            rate_sister_best_out = node_df.loc[node_df['Ks_OC'] == oc_best_value, ['Ks_Sister']]
            rate_sister_best_out_float = float(rate_sister_best_out.mean()) # trick to make it a float number

            species_sister = "_".join(sorted([latin_names[species_of_interest], sister], key=str.casefold))

            orig_mode = db.loc[species_sister]['Ortholog_Mode']
            orig_mode_sd = db.loc[species_sister]['Ortholog_Mode_SD']
            orig_median = db.loc[species_sister]['Ortholog_Median']
            orig_median_sd = db.loc[species_sister]['Ortholog_Median_SD']

            all_pairs_array.append([node, latin_names[species_of_interest], sister, round(peak_mean, 6), round(sd_err_prop, 6), round(rate_species_mean, 6), round(rate_sister_mean, 6),
                                    round(peak_best_oc_float, 6), round(sd_best_oc_float, 6), round(rate_species_best_out_float, 6), round(rate_sister_best_out_float, 6),
                                    round(orig_mode, 6), round(orig_mode_sd, 6), round(orig_median, 6), round(orig_median_sd, 6)])


    # Generating file with correction data for each divergent pair,
    # namely after obtaining a consensus value for the results coming from using different outspecies on the same divergent pair.
    all_pairs_df = DataFrame.from_records(all_pairs_array, columns=["Node", "Species", "Sister_Species",
                                                        "Peak_MeanOut", "Peak_MeanOut_SD", "Rate_Species_MeanOut", "Rate_Sister_MeanOut",
                                                        "Peak_BestOut", "Peak_BestOut_SD", "Rate_Species_BestOut", "Rate_Sister_BestOut", 
                                                        "Original_Mode", "Original_Mode_SD", "Original_Median", "Original_Median_SD"])
    with open(os.path.join("rate_adjustment", f"{species_of_interest}", f"{fcCorrect._ADJUSTMENT_TABLE.format(species_of_interest)}"),
            "w+") as outfile:
        outfile.write(all_pairs_df.to_csv(sep="\t", index=False))
        logging.info(f"Rate-adjustment results as consensus values saved in TSV format [{fcCorrect._ADJUSTMENT_TABLE.format(species_of_interest)}]")

    logging.info("")
    logging.info("All done")
