import pandas
import sys
import os
import logging
from ksrates.utils import init_logging
import ksrates.fc_plotting as fcPlot
import ksrates.fc_extract_ks_list as fc_extract_ks_list
import ksrates.fc_check_input as fcCheck
import ksrates.fc_configfile as fcConf
from ksrates.fc_rrt_correction import _ADJUSTMENT_TABLE

def plot_paralogs_distr(config_file, correction_table_file, paralog_tsv_file, anchors_ks_tsv_file):
    # INPUT
    config = fcConf.Configuration(config_file)
    init_logging("Generating mixed paralog and ortholog distributions", config.get_logging_level())
    logging.info("Loading parameters and input files")

    # GET PARAMETERS and INPUT FILES
    species = config.get_species()
    latin_names = config.get_latin_names()
    # Get analysis type
    paranome_analysis = config.get_paranome()
    colinearity_analysis = config.get_colinearity()
    if not colinearity_analysis and not paranome_analysis:
        logging.error("At least 'paranome_analysis' or 'colinearity_analysis' parameters in configuration file needs to be set to 'yes'.")
        logging.error("Exiting")
        sys.exit(1)

    # Get paralog and anchors TSV files
    # If a Ks file is not required, it will be equal to "None"
    # If the required input Ks file (paranome or anchor pairs or both) is missing, its path will be qual to an empty string ("") and the script will exit
    if paranome_analysis:
        default_path_paralog_tsv_file = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}.ks.tsv")
        paralog_tsv_file = fcCheck.get_argument_path(paralog_tsv_file, default_path_paralog_tsv_file, "Paralog Ks TSV file")
        if paralog_tsv_file == "":
            logging.error(f"Paralog Ks TSV file not found at default position [{default_path_paralog_tsv_file}].")
    if colinearity_analysis:
        default_path_anchors_tsv_file = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}.ks_anchors.tsv")
        anchors_ks_tsv_file = fcCheck.get_argument_path(anchors_ks_tsv_file, default_path_anchors_tsv_file, "Anchor pair Ks TSV file")
        if anchors_ks_tsv_file == "":
            logging.error(f"Anchor pair Ks TSV file not found at default position [{default_path_anchors_tsv_file}].")
    if paralog_tsv_file == "" or anchors_ks_tsv_file == "":
        logging.error("Exiting")
        sys.exit(1)

    # Get correction results TSV file
    # If correction_table is (still) missing, it will be equal to empty string (""), but the script will not exit
    default_path_correction_table_file = os.path.join("rate_adjustment", f"{species}", f"{_ADJUSTMENT_TABLE.format(species)}")
    correction_table_file = fcCheck.get_argument_path(correction_table_file, default_path_correction_table_file, "Correction table file")
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
    analysis_type = 'both'
    if not colinearity_analysis:
        analysis_type = 'paranome'
        logging.info(f"Plotting paranome Ks distribution for species [{species}]")
    elif not paranome_analysis:
        analysis_type = 'colinearity'
        logging.info(f"Plotting anchor pair Ks distribution for species [{species}]")
    else:
        logging.info(f"Plotting paranome and anchor pairs Ks distributions for species [{species}]")

    # PLOTTING THE BACKGROUND PARALOG AND/OR ANCHOR DISTRIBUTION
    fig_uncorr, ax_uncorr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, "un-corrected", correction_table_available, plot_correction_arrows)
    fig_corr, ax_corr = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_lim, "corrected", correction_table_available, plot_correction_arrows)

    if paranome_analysis:
        paranome_list, paranome_weights = fc_extract_ks_list.ks_list_from_tsv(paralog_tsv_file, max_ks_para, "paralogs")
        hist_paranome = fcPlot.plot_histogram("Whole-paranome (weighted)", ax_uncorr, paranome_list, bin_list, bin_width_para,
                            max_ks_para, kde_bandwidth_modifier, weight_list=paranome_weights)
        fcPlot.plot_histogram("Whole-paranome (weighted)", ax_corr, paranome_list, bin_list, bin_width_para,
                            max_ks_para, kde_bandwidth_modifier, weight_list=paranome_weights)

    if colinearity_analysis:
        anchors_list, anchors_weights = fc_extract_ks_list.ks_list_from_tsv(anchors_ks_tsv_file, max_ks_para, "anchor pairs")
        if len(anchors_list) == 0:
            logging.warning(f"No anchor pairs found! Maybe check your (gene) IDs between "
                            f"the [{species}.ks_anchors.tsv] file and the [{species}.ks.tsv] files.")
        # FOR NOW USE AN UNWEIGHTED HISTOGRAM
        hist_anchors = fcPlot.plot_histogram("Anchor pairs", ax_uncorr, anchors_list, bin_list, bin_width_para, max_ks_para,
                            kde_bandwidth_modifier, color=fcPlot.COLOR_ANCHOR_HISTOGRAM, weight_list=anchors_weights)
        fcPlot.plot_histogram("Anchor pairs", ax_corr, anchors_list, bin_list, bin_width_para, max_ks_para,
                            kde_bandwidth_modifier, color=fcPlot.COLOR_ANCHOR_HISTOGRAM, weight_list=anchors_weights)

    # Setting plot height based on tallest histogram (if paranome analysis is on, it will come from that distribution)
    if y_lim is None:
        if paranome_analysis:
            fcPlot.set_mixed_plot_height(ax_uncorr, y_lim, hist_paranome)
            fcPlot.set_mixed_plot_height(ax_corr, y_lim, hist_paranome)
        else:
            fcPlot.set_mixed_plot_height(ax_uncorr, y_lim, hist_anchors)
            fcPlot.set_mixed_plot_height(ax_corr, y_lim, hist_anchors)

    # PLOTTING THE ORTHOLOG DIVERGENCE LINES on the paralog distribution
    if correction_table_available:
        logging.info("Plotting ortholog divergence lines in the mixed plot")
        fcPlot.plot_divergences(correction_table, peak_stats, consensus_peak_for_multiple_outgroups, ax_uncorr, ax_corr, color_list, plot_correction_arrows)

    logging.info("")
    logging.info(f"Saving PDF figures of mixed plots [{fcPlot._MIXED_ADJUSTED_PLOT_FILENAME.format(species)}, {fcPlot._MIXED_UNADJUSTED_PLOT_FILENAME.format(species)}]")
    fcPlot.save_mixed_plot(fig_corr, fig_uncorr, ax_corr, ax_uncorr, species, paranome=paranome_analysis,
                        colinearity=colinearity_analysis)

    logging.info("")
    logging.info("All done")
