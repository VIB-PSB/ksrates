import pandas
import sys
import os
import logging
from ast import literal_eval
from ksrates.utils import init_logging
import ksrates.fc_plotting as fcPlot
import ksrates.fc_extract_ks_list as fc_extract_ks_list
import ksrates.fc_check_input as fcCheck
import ksrates.fc_configfile as fcConf
import ksrates.fc_consolidate_paralog_ks as fc_consolidate_paralog_ks
from ksrates.fc_rrt_correction import _ADJUSTMENT_TABLE
from ksrates.fc_plotting import _MIXED_ADJUSTED_PLOT_FILENAME, _MIXED_UNADJUSTED_PLOT_FILENAME, _OTHER_MIXED_PLOTS_SUBDIR
from ksrates.fc_wgd import _OUTPUT_KS_FILE_PATTERN_PARA, _OUTPUT_KS_FILE_PATTERN_ANCHORS, _OUTPUT_KS_FILE_PATTERN_RR_OMCL
from wgd_ksrates.viz import filter_compute_weights

def plot_paralogs_distr(config_file, expert_config_file, correction_table_file, paralog_tsv_file, anchors_ks_tsv_file, rec_ret_tsv_file):
    # INPUT
    config = fcConf.Configuration(config_file, expert_config_file)
    init_logging("Generating mixed paralog and ortholog distributions", config.get_logging_level())
    logging.info("Loading parameters and input files")

    # GET PARAMETERS and INPUT FILES
    species = config.get_species()
    latin_names = config.get_latin_names()
    # Get analysis type
    paranome_analysis = config.get_paranome()
    colinearity_analysis = config.get_colinearity()
    reciprocal_retention_analysis = config.get_reciprocal_retention()
    # Get parameters related to reciprocal retention pipeline
    num_gfs = config.get_num_reciprocal_retention_gfs(reciprocal_retention_analysis) # Number of top gene families
    bottom = config.use_bottom_gfs_instead_of_top(reciprocal_retention_analysis) # Use actually the BOTTOM GFs instead of the top ones (number of bottom GFs remains defined by "top" variable)
    rank_type = config.get_reciprocal_retention_rank_type(reciprocal_retention_analysis) # Rank type (only "lambda" supported)  
    ks_list_paralog_db_path = config.get_paralog_ks_database() # TSV file listing paralog Ks values

    # By default the pipeline uses the TOP-ranked reciprocally retained GFs.
    # However, for comparison purposes, the user might want to use the BOTTOM GFs (e.g. the bottom 2000 ones)
    if bottom == False:
        top_or_bottom = "top"
    else:
        top_or_bottom = "bottom"

    if not colinearity_analysis and not paranome_analysis and not reciprocal_retention_analysis:
        logging.error('At least one of the "paranome" or "collinearity" or "reciprocal_retention" parameters in the configuration file needs to be set to "yes".')
        logging.error("Exiting.")
        sys.exit(1)

    # GET PARALOG KS DATA, EITHER FROM DATABASE OF FROM ORIGINAL KS TSV FILES

    # Assume there is no database or no species data in it
    paralog_db_available = False
    paranome_ks_list_db = None

    # If database was configured, load from it paralog Ks lists
    if ks_list_paralog_db_path:
        try:
            # Create database file if it doesn't exist yet
            if not os.path.isfile(ks_list_paralog_db_path):
                logging.info(f"Paralog Ks database not found at [{ks_list_paralog_db_path}]. Creating...")
                fc_consolidate_paralog_ks.initialize_paralog_db(ks_list_paralog_db_path)

            with open(ks_list_paralog_db_path, "r") as f:
                paranome_ks_list_db = pandas.read_csv(f, sep="\t", index_col=0)

            # If species not in database, add it using its TSV files, then reload
            species_in_db = species in paranome_ks_list_db.index
            if not species_in_db:
                logging.info(f"Species [{species}] not found in database. Populating it from TSV files...")
                ks_data = fc_consolidate_paralog_ks.extract_paralog_ks_from_tsv(
                    species, paranome_analysis, colinearity_analysis, reciprocal_retention_analysis,
                    num_gfs=num_gfs, rank_type=rank_type, bottom=bottom
                )
                fc_consolidate_paralog_ks.write_to_paralog_db(species, ks_data, ks_list_paralog_db_path)
                with open(ks_list_paralog_db_path, "r") as f:
                    paranome_ks_list_db = pandas.read_csv(f, sep="\t", index_col=0)

            # When imported from csv format, all the Ks lists in the df are read as
            # plain text (strings) and must be converted back to value lists
            logging.info("Loading paralog Ks list database...")
            for col in ['Ks_paranome', 'Ks_paranome_weights', 'Ks_anchors', 'Ks_anchors_weights',
                        'Ks_reciprocally_retained', 'Ks_reciprocally_retained_weights']:
                if col in paranome_ks_list_db.columns:
                    paranome_ks_list_db.loc[:, col] = paranome_ks_list_db.loc[:, col].apply(
                        lambda x: literal_eval(x) if x is not None and x != 'None' else None)

            paralog_db_available = True

        except Exception as e:
            logging.warning(f"Could not use paralog Ks database [{ks_list_paralog_db_path}]: {str(e)}. Will read TSV files instead.")
            paralog_db_available = False

    # If the species' Ks data is not available in the database, use the original paralog TSV files as fallback
    # Determine which TSV files are needed based on:
    #   1. Whether the analysis type is enabled AND
    #   2. Whether data for that type is missing from the database
    # Example: TSV file IS required if the analysis type is enabled BUT its data is missing from database
    # Example: paralog_tsv_file is NOT required if paranome_analysis=True AND database contains valid Ks_paranome data
    paralog_tsv_file_required = paranome_analysis and not (paralog_db_available and paranome_ks_list_db is not None and species in paranome_ks_list_db.index and paranome_ks_list_db.loc[species, 'Ks_paranome'] is not None)
    anchors_tsv_file_required = colinearity_analysis and not (paralog_db_available and paranome_ks_list_db is not None and species in paranome_ks_list_db.index and paranome_ks_list_db.loc[species, 'Ks_anchors'] is not None)
    recret_tsv_file_required = reciprocal_retention_analysis and not (paralog_db_available and paranome_ks_list_db is not None and species in paranome_ks_list_db.index and paranome_ks_list_db.loc[species, 'Ks_reciprocally_retained'] is not None)

    if paralog_tsv_file_required:
        default_path_paralog_tsv_file = os.path.join("paralog_distributions", f"wgd_{species}", _OUTPUT_KS_FILE_PATTERN_PARA.format(species))
        paralog_tsv_file = fcCheck.get_argument_path(paralog_tsv_file, default_path_paralog_tsv_file, "Paralog Ks TSV file")
        if paralog_tsv_file == "":
            logging.error(f"Paralog Ks TSV file not found at default position [{default_path_paralog_tsv_file}].")
    if anchors_tsv_file_required:
        default_path_anchors_tsv_file = os.path.join("paralog_distributions", f"wgd_{species}", _OUTPUT_KS_FILE_PATTERN_ANCHORS.format(species))
        anchors_ks_tsv_file = fcCheck.get_argument_path(anchors_ks_tsv_file, default_path_anchors_tsv_file, "Anchor pair Ks TSV file")
        if anchors_ks_tsv_file == "":
            logging.error(f"Anchor pair Ks TSV file not found at default position [{default_path_anchors_tsv_file}].")
    if recret_tsv_file_required:
        default_path_rec_ret_tsv_file = os.path.join("paralog_distributions", f"wgd_{species}", _OUTPUT_KS_FILE_PATTERN_RR_OMCL.format(species, top_or_bottom, num_gfs))
        rec_ret_tsv_file = fcCheck.get_argument_path(rec_ret_tsv_file, default_path_rec_ret_tsv_file, "Reciprocally retained paralog Ks TSV file")
        if rec_ret_tsv_file == "":
            logging.error(f"Reciprocally retained paralog Ks TSV file not found at default position [{default_path_rec_ret_tsv_file}].")

    # Exit if both the database source and the original Ks TSV source have failed
    if (paralog_tsv_file_required and paralog_tsv_file == "") or (anchors_tsv_file_required and anchors_ks_tsv_file == "") or (recret_tsv_file_required and rec_ret_tsv_file == ""):
        logging.error("Exiting")
        sys.exit(1)


    # Creating folders for output files
    output_folder = os.path.join("rate_adjustment", f"{species}")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        logging.info(f"Generating output folder [{output_folder}]")

    # Get correction results TSV file
    # If correction_table is (still) missing, it will be equal to empty string (""), but the script will not exit
    default_path_correction_table_file = os.path.join(output_folder, f"{_ADJUSTMENT_TABLE.format(species)}")
    correction_table_file = fcCheck.get_argument_path(correction_table_file, default_path_correction_table_file, "Rate-adjustment table file")
    if correction_table_file == "": # it means that the correction_table is not present or available yet
        logging.warning("Rate-adjustment data are not available yet, only paralog distribution will be plotted.")
        correction_table_available = False
    else:
        with open(correction_table_file, "r") as f:
            correction_table = pandas.read_csv(f, sep="\t")
            if correction_table.shape[0] == 0:
                logging.warning(f"Rate-adjustment table file is present but doesn't contain any data: rate-adjusted divergences will not be plotted.")
                correction_table_available = False
            else:
                correction_table_available = True

    # Getting the choice on how to deal with the presence of multiple corrections for the same divergent pair
    # due to the use of multiple trios/outgroup during correction
    # Available options:
    #  - 'mean among outgroups': taking the average of the corrected peaks
    #  - 'best outgroup': taking the corrected peak coming from the best outgroup, which is the one with smallest OC segment
    consensus_peak_for_multiple_outgroups = config.get_consensus_peak_for_multiple_outgroups()

    # Get other parameters
    max_ks_para = config.get_max_ks_para()
    min_ks_anchors = config.get_min_ks_anchors()
    bin_width_para = config.get_bin_width_para()
    bin_list = fcPlot.get_bins(max_ks_para, bin_width_para)
    x_max_lim = config.get_x_max_lim()
    y_lim = config.get_y_lim()  # by default it's "None"
    color_list = config.get_color_list()  # colors for vertical lines depicting speciation events
    kde_bandwidth_modifier = config.get_kde_bandwidth_modifier() # get the modifier to adjust KDE fitting
    plot_correction_arrows = config.plot_correction_arrows()
    peak_stats = config.get_peak_stats() # default is mode (other option, median)

    # -----------------------------------------------------------------------------

    # GENERATING MIXED DISTRIBUTION PLOT
    logging.info("")
    ax_uncorr_list, ax_corr_list = [], []
    ax_uncorr_include_para, ax_corr_include_para= [], []
    ax_uncorr_include_col, ax_corr_include_col= [], []
    ax_uncorr_include_rr, ax_corr_include_rr= [], []

    # Generate the mixed plot with a single data type (e.g. paranome) for any required data_type
    if paranome_analysis:
        logging.info(f"Plotting paranome Ks distribution for species [{species}]")
        fig_uncorr_para, ax_uncorr_para = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "un-corrected", correction_table_available, plot_correction_arrows,
                                    paranome_data=paranome_analysis, num_gfs=num_gfs, rank_type=rank_type)
        fig_corr_para, ax_corr_para = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "corrected", correction_table_available, plot_correction_arrows,
                                    paranome_data=paranome_analysis, num_gfs=num_gfs, rank_type=rank_type)
        ax_uncorr_include_para.append(ax_uncorr_para)
        ax_corr_include_para.append(ax_corr_para)
        ax_uncorr_list.append(ax_uncorr_para)
        ax_corr_list.append(ax_corr_para)

    if colinearity_analysis:
        logging.info(f"Plotting anchor pair Ks distribution for species [{species}]")
        fig_uncorr_col, ax_uncorr_col = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "un-corrected", correction_table_available, plot_correction_arrows,
                                    colinearity_data=colinearity_analysis, num_gfs=num_gfs, rank_type=rank_type)
        fig_corr_col, ax_corr_col = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "corrected", correction_table_available, plot_correction_arrows,
                                    colinearity_data=colinearity_analysis, num_gfs=num_gfs, rank_type=rank_type)
        ax_uncorr_include_col.append(ax_uncorr_col)
        ax_corr_include_col.append(ax_corr_col)
        ax_uncorr_list.append(ax_uncorr_col)
        ax_corr_list.append(ax_corr_col)

    if reciprocal_retention_analysis:
        logging.info(f"Plotting reciprocally retained paralog Ks distribution for species [{species}]")
        fig_uncorr_rr, ax_uncorr_rr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "un-corrected", correction_table_available, plot_correction_arrows,
                                    reciprocal_retention_data=reciprocal_retention_analysis, num_gfs=num_gfs, rank_type=rank_type)
        fig_corr_rr, ax_corr_rr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "corrected", correction_table_available, plot_correction_arrows,
                                    reciprocal_retention_data=reciprocal_retention_analysis, num_gfs=num_gfs, rank_type=rank_type)
        ax_uncorr_include_rr.append(ax_uncorr_rr)
        ax_corr_include_rr.append(ax_corr_rr)
        ax_uncorr_list.append(ax_uncorr_rr)
        ax_corr_list.append(ax_corr_rr)

    # Generate the mixed plot for any required pair of data types (e.g. paranome and anchors)
    if paranome_analysis and colinearity_analysis:
        logging.info(f"Plotting paranome and anchor pairs Ks distributions for species [{species}]")
        fig_uncorr_para_col, ax_uncorr_para_col = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "un-corrected", correction_table_available, plot_correction_arrows,
                                    paranome_data=paranome_analysis, colinearity_data=colinearity_analysis,
                                    num_gfs=num_gfs, rank_type=rank_type)
        fig_corr_para_col, ax_corr_para_col = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "corrected", correction_table_available, plot_correction_arrows,
                                    paranome_data=paranome_analysis, colinearity_data=colinearity_analysis,
                                    num_gfs=num_gfs, rank_type=rank_type)
        ax_uncorr_include_para.append(ax_uncorr_para_col)
        ax_corr_include_para.append(ax_corr_para_col)
        ax_uncorr_include_col.append(ax_uncorr_para_col)
        ax_corr_include_col.append(ax_corr_para_col)
        ax_uncorr_list.append(ax_uncorr_para_col)
        ax_corr_list.append(ax_corr_para_col)

    if paranome_analysis and reciprocal_retention_analysis:
        logging.info(f"Plotting paranome and reciprocally retained paralog Ks distributions for species [{species}]")
        fig_uncorr_para_rr, ax_uncorr_para_rr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "un-corrected", correction_table_available, plot_correction_arrows,
                                    paranome_data=paranome_analysis, reciprocal_retention_data=reciprocal_retention_analysis,
                                    num_gfs=num_gfs, rank_type=rank_type)
        fig_corr_para_rr, ax_corr_para_rr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "corrected", correction_table_available, plot_correction_arrows,
                                    paranome_data=paranome_analysis, reciprocal_retention_data=reciprocal_retention_analysis,
                                    num_gfs=num_gfs, rank_type=rank_type)
        ax_uncorr_include_para.append(ax_uncorr_para_rr)
        ax_corr_include_para.append(ax_corr_para_rr)
        ax_uncorr_include_rr.append(ax_uncorr_para_rr)
        ax_corr_include_rr.append(ax_corr_para_rr)
        ax_uncorr_list.append(ax_uncorr_para_rr)
        ax_corr_list.append(ax_corr_para_rr)

    if colinearity_analysis and reciprocal_retention_analysis:
        logging.info(f"Plotting anchor pairs and reciprocally retained paralog Ks distributions for species [{species}]")
        fig_uncorr_col_rr, ax_uncorr_col_rr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "un-corrected", correction_table_available, plot_correction_arrows,
                                    colinearity_data=colinearity_analysis, reciprocal_retention_data=reciprocal_retention_analysis,
                                    num_gfs=num_gfs, rank_type=rank_type)
        fig_corr_col_rr, ax_corr_col_rr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "corrected", correction_table_available, plot_correction_arrows,
                                    colinearity_data=colinearity_analysis, reciprocal_retention_data=reciprocal_retention_analysis,
                                    num_gfs=num_gfs, rank_type=rank_type)
        ax_uncorr_include_col.append(ax_uncorr_col_rr)
        ax_corr_include_col.append(ax_corr_col_rr)
        ax_uncorr_include_rr.append(ax_uncorr_col_rr)
        ax_corr_include_rr.append(ax_corr_col_rr)
        ax_uncorr_list.append(ax_uncorr_col_rr)
        ax_corr_list.append(ax_corr_col_rr)

    # Generate the mixed plot for the three data types altogether if required by configuration
    if paranome_analysis and colinearity_analysis and reciprocal_retention_analysis:
        logging.info(f"Plotting paranome, anchor pairs and reciprocally retained paralog Ks distributions for species [{species}]")
        fig_uncorr_para_col_rr, ax_uncorr_para_col_rr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "un-corrected", correction_table_available, plot_correction_arrows,
                                    paranome_data=paranome_analysis, colinearity_data=colinearity_analysis,
                                    reciprocal_retention_data=reciprocal_retention_analysis, num_gfs=num_gfs, rank_type=rank_type)
        fig_corr_para_col_rr, ax_corr_para_col_rr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, 
                                    "corrected", correction_table_available, plot_correction_arrows,
                                    paranome_data=paranome_analysis, colinearity_data=colinearity_analysis,
                                    reciprocal_retention_data=reciprocal_retention_analysis, num_gfs=num_gfs, rank_type=rank_type)
        ax_uncorr_include_para.append(ax_uncorr_para_col_rr)
        ax_corr_include_para.append(ax_corr_para_col_rr)
        ax_uncorr_include_col.append(ax_uncorr_para_col_rr)
        ax_corr_include_col.append(ax_corr_para_col_rr)
        ax_uncorr_include_rr.append(ax_uncorr_para_col_rr)
        ax_corr_include_rr.append(ax_corr_para_col_rr)
        ax_uncorr_list.append(ax_uncorr_para_col_rr)
        ax_corr_list.append(ax_corr_para_col_rr)

    # PLOTTING THE BACKGROUND PARALOG DISTRIBUTION(S)
    if paranome_analysis:
        paranome_list = None
        paranome_weights = None
        # Try database first if available
        if paralog_db_available and paranome_ks_list_db is not None and species in paranome_ks_list_db.index:
            db_paranome = paranome_ks_list_db.loc[species, 'Ks_paranome']
            if db_paranome is not None:
                # Reconstruct dataframe from database and recalculate weights
                paranome_df = pandas.DataFrame(db_paranome)
                paranome_df_recalc = filter_compute_weights(paranome_df, min_ks=0.005, max_ks=max_ks_para)
                paranome_list = paranome_df_recalc["Ks"].to_list()
                paranome_weights = paranome_df_recalc["WeightOutliersExcluded"].to_list()
                logging.info(f"Using paranome Ks list from database for [{species}] (weights recalculated for Ks <= {max_ks_para})")
        # Fall back to TSV file if database not available or empty
        if paranome_list is None:
            # Get paranome Ks values within the requested range and recalculate their associated weight
            paranome_list, paranome_weights = fc_extract_ks_list.ks_list_from_tsv(paralog_tsv_file, max_ks_para, "paralogs")
        
        for ax_uncorr in ax_uncorr_include_para:
            hist_paranome = fcPlot.plot_histogram("Whole-paranome", ax_uncorr, paranome_list, bin_list, bin_width_para,
                                max_ks_para, kde_bandwidth_modifier, weight_list=paranome_weights)
        for ax_corr in ax_corr_include_para:
            fcPlot.plot_histogram("Whole-paranome", ax_corr, paranome_list, bin_list, bin_width_para,
                            max_ks_para, kde_bandwidth_modifier, weight_list=paranome_weights)

    if colinearity_analysis:
        anchors_list = None
        anchors_weights = None
        # Try database first if available
        if paralog_db_available and paranome_ks_list_db is not None and species in paranome_ks_list_db.index:
            db_anchors = paranome_ks_list_db.loc[species, 'Ks_anchors']
            if db_anchors is not None:
                # Reconstruct dataframe from database and recalculate weights
                anchors_df = pandas.DataFrame(db_anchors)
                anchors_df_recalc = filter_compute_weights(anchors_df, min_ks=min_ks_anchors, max_ks=max_ks_para)
                anchors_list = anchors_df_recalc["Ks"].to_list()
                anchors_weights = anchors_df_recalc["WeightOutliersExcluded"].to_list()
                logging.info(f"Using anchor Ks list from database for [{species}] (weights recalculated for {min_ks_anchors} <= Ks <= {max_ks_para})")
        # Fall back to TSV file if database not available or empty
        if anchors_list is None:
            # Get anchor pair Ks values within the requested range (using min_ks_anchors) and recalculate their associated weight
            anchors_list, anchors_weights = fc_extract_ks_list.ks_list_from_tsv(anchors_ks_tsv_file, max_ks_para, "anchor pairs", min_ks=min_ks_anchors)

        if len(anchors_list) == 0:
            logging.warning(f"No anchor pairs found! Maybe check your (gene) IDs between "
                            f"anchor pairs file [{_OUTPUT_KS_FILE_PATTERN_ANCHORS.format(species)}] and"
                            f"whole-paranome file [{_OUTPUT_KS_FILE_PATTERN_PARA.format(species)}]")
        for ax_uncorr in ax_uncorr_include_col:
            hist_anchors = fcPlot.plot_histogram("Anchor pairs", ax_uncorr, anchors_list, bin_list, bin_width_para, max_ks_para,
                                kde_bandwidth_modifier, color=fcPlot.COLOR_ANCHOR_HISTOGRAM, weight_list=anchors_weights)
        for ax_corr in ax_corr_include_col:
            fcPlot.plot_histogram("Anchor pairs", ax_corr, anchors_list, bin_list, bin_width_para, max_ks_para,
                                kde_bandwidth_modifier, color=fcPlot.COLOR_ANCHOR_HISTOGRAM, weight_list=anchors_weights)

    if reciprocal_retention_analysis:
        rec_ret_list = None
        rec_ret_weights = None
        # Try database first if available
        if paralog_db_available and paranome_ks_list_db is not None and species in paranome_ks_list_db.index:
            db_recret = paranome_ks_list_db.loc[species, 'Ks_reciprocally_retained']
            if db_recret is not None:
                # Reconstruct dataframe from database and recalculate weights
                recret_df = pandas.DataFrame(db_recret)
                recret_df_recalc = filter_compute_weights(recret_df, min_ks=0.005, max_ks=max_ks_para)
                rec_ret_list = recret_df_recalc["Ks"].to_list()
                rec_ret_weights = recret_df_recalc["WeightOutliersExcluded"].to_list()
                logging.info(f"Using reciprocally retained Ks list from database for [{species}] (weights recalculated for 0.005 <= Ks <= {max_ks_para})")
        # Fall back to TSV file if database not available or empty
        if rec_ret_list is None:
            # Get recret Ks values within the requested range and recalculate their associated weight
            rec_ret_list, rec_ret_weights = fc_extract_ks_list.ks_list_from_tsv(rec_ret_tsv_file, max_ks_para, "reciprocally retained")
        
        for ax_uncorr in ax_uncorr_include_rr:
            hist_rec_ret = fcPlot.plot_histogram("Reciprocally retained paralogs", ax_uncorr, rec_ret_list, bin_list, bin_width_para,
                                max_ks_para, kde_bandwidth_modifier, color=fcPlot.COLOR_REC_RET_HISTOGRAM, weight_list=rec_ret_weights)
        for ax_corr in ax_corr_include_rr:
                fcPlot.plot_histogram("Reciprocally retained paralogs", ax_corr, rec_ret_list, bin_list, bin_width_para,
                                    max_ks_para, kde_bandwidth_modifier, color=fcPlot.COLOR_REC_RET_HISTOGRAM, weight_list=rec_ret_weights)

    # Setting plot height based on tallest histogram among the paralog distribution(s)
    if y_lim is None:
        # Set the height based on paranome distribution when this is for sure the tallest
        if paranome_analysis:
            for ax_uncorr, ax_corr in zip(ax_uncorr_include_para, ax_corr_include_para):
                fcPlot.set_mixed_plot_height(ax_uncorr, y_lim, hist_paranome)
                fcPlot.set_mixed_plot_height(ax_corr, y_lim, hist_paranome)

        # Set height based on anchor pair distributions if that's the only one in the plot
        if colinearity_analysis:
            fcPlot.set_mixed_plot_height(ax_uncorr_col, y_lim, hist_anchors)
            fcPlot.set_mixed_plot_height(ax_corr_col, y_lim, hist_anchors)
            
        # Set height based on reciprocally retained paralog distribution if that's the only one in the plot 
        if reciprocal_retention_analysis:
            fcPlot.set_mixed_plot_height(ax_uncorr_rr, y_lim, hist_rec_ret)
            fcPlot.set_mixed_plot_height(ax_corr_rr, y_lim, hist_rec_ret)

        # Set height to tallest distribution between anchor pair and reciprocally retained gene families
        if colinearity_analysis and reciprocal_retention_analysis:
            fcPlot.set_mixed_plot_height(ax_uncorr_col_rr, y_lim, hist_anchors, hist_rec_ret)
            fcPlot.set_mixed_plot_height(ax_corr_col_rr, y_lim, hist_anchors, hist_rec_ret)


    # PLOTTING THE ORTHOLOG DIVERGENCE LINES on the paralog distribution
    if correction_table_available:
        logging.info("Plotting ortholog divergence lines in the mixed plot")
        for ax_uncorr, ax_corr in zip(ax_uncorr_list, ax_corr_list):
            fcPlot.plot_divergences(correction_table, peak_stats, consensus_peak_for_multiple_outgroups, ax_uncorr, ax_corr, color_list, plot_correction_arrows)

    logging.info("")
    logging.info(f"Saving PDF figures of mixed plots")
    if paranome_analysis:
        fcPlot.save_mixed_plot(fig_corr_para, fig_uncorr_para, ax_corr_para, ax_uncorr_para, species, correction_table_available, paranome=paranome_analysis,
                            colinearity=False, reciprocal_retention=False,
                            filename_adjusted=_MIXED_ADJUSTED_PLOT_FILENAME.format(species, "paranome"),
                            filename_unadjusted=_MIXED_UNADJUSTED_PLOT_FILENAME.format(species, "paranome"))

    if colinearity_analysis:
        fcPlot.save_mixed_plot(fig_corr_col, fig_uncorr_col, ax_corr_col, ax_uncorr_col, species, correction_table_available, paranome=False,
                            colinearity=colinearity_analysis, reciprocal_retention=False,
                            filename_adjusted=_MIXED_ADJUSTED_PLOT_FILENAME.format(species, "anchors"),
                            filename_unadjusted=_MIXED_UNADJUSTED_PLOT_FILENAME.format(species, "anchors"))

    if reciprocal_retention_analysis:
        fcPlot.save_mixed_plot(fig_corr_rr, fig_uncorr_rr, ax_corr_rr, ax_uncorr_rr, species, correction_table_available, paranome=False,
                            colinearity=False, reciprocal_retention=reciprocal_retention_analysis,
                            filename_adjusted=_MIXED_ADJUSTED_PLOT_FILENAME.format(species, f"recret"),
                            filename_unadjusted=_MIXED_UNADJUSTED_PLOT_FILENAME.format(species, f"recret"))
   
    if paranome_analysis and colinearity_analysis:
        fcPlot.save_mixed_plot(fig_corr_para_col, fig_uncorr_para_col, ax_corr_para_col, ax_uncorr_para_col, species, correction_table_available,
                            paranome=paranome_analysis, colinearity=colinearity_analysis, reciprocal_retention=False,
                            filename_adjusted=os.path.join(_OTHER_MIXED_PLOTS_SUBDIR, _MIXED_ADJUSTED_PLOT_FILENAME.format(species, "paranome_anchors")),
                            filename_unadjusted=os.path.join(_OTHER_MIXED_PLOTS_SUBDIR, _MIXED_UNADJUSTED_PLOT_FILENAME.format(species, "paranome_anchors")))
    
    if paranome_analysis and reciprocal_retention_analysis:
        fcPlot.save_mixed_plot(fig_corr_para_rr, fig_uncorr_para_rr, ax_corr_para_rr, ax_uncorr_para_rr, species, correction_table_available, 
                            paranome=paranome_analysis, colinearity=False, reciprocal_retention=reciprocal_retention_analysis,
                            filename_adjusted=os.path.join(_OTHER_MIXED_PLOTS_SUBDIR, _MIXED_ADJUSTED_PLOT_FILENAME.format(species, f"paranome_recret")),
                            filename_unadjusted=os.path.join(_OTHER_MIXED_PLOTS_SUBDIR, _MIXED_UNADJUSTED_PLOT_FILENAME.format(species, f"paranome_recret")))

    if colinearity_analysis and reciprocal_retention_analysis:
        fcPlot.save_mixed_plot(fig_corr_col_rr, fig_uncorr_col_rr, ax_corr_col_rr, ax_uncorr_col_rr, species, correction_table_available, paranome=False,
                            colinearity=reciprocal_retention_analysis, reciprocal_retention=reciprocal_retention_analysis,
                            filename_adjusted=os.path.join(_OTHER_MIXED_PLOTS_SUBDIR, _MIXED_ADJUSTED_PLOT_FILENAME.format(species, f"anchors_recret")),
                            filename_unadjusted=os.path.join(_OTHER_MIXED_PLOTS_SUBDIR, _MIXED_UNADJUSTED_PLOT_FILENAME.format(species, f"anchors_recret")))

    if paranome_analysis and colinearity_analysis and reciprocal_retention_analysis:
        fcPlot.save_mixed_plot(fig_corr_para_col_rr, fig_uncorr_para_col_rr, ax_corr_para_col_rr, ax_uncorr_para_col_rr, species, correction_table_available, 
                            paranome=paranome_analysis, colinearity=reciprocal_retention_analysis, reciprocal_retention=reciprocal_retention_analysis,
                            filename_adjusted=os.path.join(_OTHER_MIXED_PLOTS_SUBDIR, _MIXED_ADJUSTED_PLOT_FILENAME.format(species, f"paranome_anchors_recret")),
                            filename_unadjusted=os.path.join(_OTHER_MIXED_PLOTS_SUBDIR, _MIXED_UNADJUSTED_PLOT_FILENAME.format(species, f"paranome_anchors_recret")))

    logging.info("")
    logging.info("All done")
