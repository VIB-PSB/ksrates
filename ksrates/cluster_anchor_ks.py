import os
import sys
import copy
import logging
import matplotlib.pyplot as plt
from numpy import array, arange
from pandas import read_csv
from statistics import mean, median
from scipy import stats
import ksrates.fc_cluster_anchors as fcCluster
import ksrates.fc_configfile as fcConf
import ksrates.fc_plotting as fcPlot
import ksrates.fc_check_input as fcCheck
import ksrates.fc_extract_ks_list as fc_extract_ks_list
from ksrates.utils import init_logging
from ksrates.fc_cluster_anchors import subfolder
from ksrates.fc_rrt_correction import _ADJUSTMENT_TABLE

def cluster_anchor_ks(config_file, correction_table_file, path_anchorpoints_txt, path_multiplicons_txt, path_segments_txt, path_list_elements_txt, path_ks_anchor_file, path_multiplicon_pair_txt):
    config = fcConf.Configuration(config_file)
    init_logging(f"Clustering anchorpoints Ks values to reconstruct recent WGD events", config.get_logging_level())
    logging.info("Loading parameters and input files")

    colinearity_analysis = config.get_colinearity()
    if not colinearity_analysis:
        logging.warning("Colinearity non required in configuration file; the anchor Ks clustering will be skipped.")
        sys.exit(0) # exit code 0 because no actual errors were thrown

    species = config.get_species()
    latin_names = config.get_latin_names()
    latin_name = latin_names[species].replace(' ', '\ ') # replace the space to have the correct format in titles
    max_ks_para = config.get_max_ks_para()
    bin_width = config.get_bin_width_para()
    bin_list = fcPlot.get_bins(max_ks_para, bin_width)
    x_max_lim = config.get_x_max_lim()
    y_max_lim = config.get_y_lim()
    color_list = config.get_color_list()
    plot_correction_arrows = config.plot_correction_arrows()

    max_EM_iterations = config.get_max_EM_iterations() # default 300
    num_EM_initializations = config.get_num_EM_initializations() # how many times the fitting with N given components is initialized 
    logging.info(f" - maximum EM iterations: {max_EM_iterations}")
    logging.info(f" - number of EM initializations: {num_EM_initializations}")

    # Getting the statistical measure for how to determine the representative value of an ortholog distribution
    peak_stats = config.get_peak_stats() # default is mode (other option, median)

    # Getting the choice on how to deal with the presence of multiple corrections for the same divergent pair
    consensus_peak_for_multiple_outgroups = config.get_consensus_peak_for_multiple_outgroups()

    # Checking user-defined path and / or default path for each required input file
    default_path_correction_table_file = os.path.join("rate_adjustment", f"{species}", f"{_ADJUSTMENT_TABLE.format(species)}")
    correction_table_file = fcCheck.get_argument_path(correction_table_file, default_path_correction_table_file, "Correction table file")
    if correction_table_file == "":
        logging.warning("Rate-adjustment data are not available yet, only anchor pair distribution will be plotted.")
        correction_table_available = False
    else:
        with open(correction_table_file, "r") as f:
            correction_table = read_csv(f, sep="\t")
            correction_table_available = True

    default_path_anchorpoints_txt = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}_i-adhore", "anchorpoints.txt")
    path_anchorpoints_txt = fcCheck.get_argument_path(path_anchorpoints_txt, default_path_anchorpoints_txt, "anchorpoints.txt file")
    if path_anchorpoints_txt == "":
        logging.error(f"anchorpoints.txt file not found at default position [{default_path_anchorpoints_txt}].")

    default_path_multiplicons_txt = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}_i-adhore", "multiplicons.txt")
    path_multiplicons_txt = fcCheck.get_argument_path(path_multiplicons_txt, default_path_multiplicons_txt, "multiplicons.txt file")
    if path_multiplicons_txt == "":
        logging.error(f"multiplicons.txt file not found at default position [{default_path_anchorpoints_txt}].")

    default_path_segments_txt = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}_i-adhore", "segments.txt")
    path_segments_txt = fcCheck.get_argument_path(path_segments_txt, default_path_segments_txt, "segments.txt file")
    if path_segments_txt == "":
        logging.error(f"segments.txt file not found at default position [{default_path_segments_txt}].")

    default_path_multiplicon_pair_txt = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}_i-adhore", "multiplicon_pairs.txt")
    path_multiplicon_pair_txt = fcCheck.get_argument_path(path_multiplicon_pair_txt, default_path_multiplicon_pair_txt, "multiplicon_pairs.txt file")
    if path_multiplicon_pair_txt == "":
        logging.error(f"multiplicon_pairs.txt file not found at default position [{default_path_multiplicon_pair_txt}].")

    default_path_list_elements_txt = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}_i-adhore", "list_elements.txt")
    path_list_elements_txt = fcCheck.get_argument_path(path_list_elements_txt, default_path_list_elements_txt, "list_elements.txt file")
    if path_list_elements_txt == "":
        logging.error(f"list_elements.txt file not found at default position [{default_path_list_elements_txt}].")

    default_path_ks_anchor_file = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}.ks_anchors.tsv")
    path_ks_anchor_file = fcCheck.get_argument_path(path_ks_anchor_file, default_path_ks_anchor_file, f"{species}.ks_anchors.tsv file")
    if path_ks_anchor_file == "":
        logging.error(f"{species}.ks_anchors.tsv file not found at default position [{default_path_ks_anchor_file}].")

    if path_anchorpoints_txt == "" or path_multiplicons_txt == "" or path_segments_txt == "" or path_multiplicon_pair_txt == "" or path_list_elements_txt == "" or path_ks_anchor_file == "":
        logging.error("Exiting")
        sys.exit(1)

    # Creating folder for secondary output files
    output = os.path.join(subfolder)
    if not os.path.isdir(os.path.join("rate_adjustment", species, output)):
        logging.info(f"Creating directory [rate_adjustment/{species}/{output}]")
        os.makedirs(os.path.join("rate_adjustment", species, output))

    # Parsing I-ADHoRe output files to get information about multiplicons, multiplicon levels, segments, anchorpoints, anchor pairs and Ks values. 
    segments_per_multip = fcCluster.parse_segments_file(path_segments_txt)

    segments_from_gene = fcCluster.parse_list_elements(path_list_elements_txt)

    ks_anchors = fcCluster.parse_ks_anchors_tsv_file(path_ks_anchor_file)

    multipl_per_level, level_of_each_multipl, max_level, level_list, level_list_filtered = fcCluster.parse_multiplicons_file(path_multiplicons_txt)

    anchors_per_multipl, levels_of_each_anchor = fcCluster.parse_multiplicon_pairs_file(path_multiplicon_pair_txt, level_of_each_multipl)

    anchorpoints_per_multipl, multipl_per_anchorpoint, levels_of_anchorpoints = fcCluster.parse_anchorpoints_file(path_anchorpoints_txt, level_of_each_multipl)

    anchor_ks_list = fc_extract_ks_list.ks_list_from_tsv(path_ks_anchor_file, max_ks_para, "anchor pairs") # Get complete anchor Ks list to be plotted in the background

    # -----------------------------------------------------------------------------

    # Getting non-redundant segment pairs based on group of anchorpoints
    # Filtering away segment pairs whose anchorpoints list is a subset of another segment pair Ks list
    logging.info("")
    logging.info(f"Obtaining a non-redundant list of anchorpoints Ks values")
    logging.info("")

    anchorpoints_ks_per_segment_pair = {} # Ks values per segment pair
    anchorpoints_per_segment_pair = {} # anchorpoint names per segment pair

    for multipl_id in level_of_each_multipl: # for each multiplicon
        anchorpoints_list = anchorpoints_per_multipl[multipl_id]
        segments_current_multip = segments_per_multip[multipl_id]
        for anchorpoints in anchorpoints_list:

            # Process only anchorpoints with acceptable Ks value
            if anchorpoints in ks_anchors:
                ks = ks_anchors[anchorpoints]
                if ks != "" and 0.05 <= float(ks) <= max_ks_para:
                    anchorpoint1, anchorpoint2 = anchorpoints[0], anchorpoints[1]
                    # Get the two segments where the two anchorpoints are placed in the current multiplicon
                    all_segments_anchorpoint1 = segments_from_gene[anchorpoint1]
                    all_segments_anchorpoint2 = segments_from_gene[anchorpoint2]
                    anchorpoint1_segment = list(set(all_segments_anchorpoint1) & set(segments_current_multip))[0]
                    anchorpoint2_segment = list(set(all_segments_anchorpoint2) & set(segments_current_multip))[0]
                    current_segment_pair_sorted = tuple(sorted([anchorpoint1_segment, anchorpoint2_segment]))
                    
                    if current_segment_pair_sorted not in anchorpoints_per_segment_pair:
                        anchorpoints_per_segment_pair[current_segment_pair_sorted] = [anchorpoints]
                    else:
                        anchorpoints_per_segment_pair[current_segment_pair_sorted].append(anchorpoints)  


    other_segment_pairs_dict = {} # links to each other the segment pairs that share anchorpoints
    # {key = segment pair; value = {key= other segment pair in which the anchorpoint is present; value= shared anchorpoints}
    # e.g. { (3, 5): { (1287, 1288): [shared anchorpoints] } }
    for segment_pair in anchorpoints_per_segment_pair:
        for anchorpoint in anchorpoints_per_segment_pair[segment_pair]:
            anchorpoint1, anchorpoint2 = anchorpoint[0], anchorpoint[1]
            # Get all multiplicons in which the current anchorpoint is found
            multipl_of_current_anchorpoint = multipl_per_anchorpoint[anchorpoint]

            for multipl_id in multipl_of_current_anchorpoint:
                # Let's take the segment pairs where this anchorpoints lies on
                segments_current_multip = segments_per_multip[multipl_id]
                all_segments_anchorpoint1 = segments_from_gene[anchorpoint1]
                all_segments_anchorpoint2 = segments_from_gene[anchorpoint2]
                anchorpoint1_segment = list(set(all_segments_anchorpoint1) & set(segments_current_multip))[0]
                anchorpoint2_segment = list(set(all_segments_anchorpoint2) & set(segments_current_multip))[0]
                current_segment_pair_sorted = tuple(sorted([anchorpoint1_segment, anchorpoint2_segment]))

                if segment_pair not in other_segment_pairs_dict:
                    other_segment_pairs_dict[segment_pair] = {}

                if current_segment_pair_sorted != segment_pair: # if the anchorpoint is also in another segment pair
                    if current_segment_pair_sorted not in other_segment_pairs_dict[segment_pair]:
                        other_segment_pairs_dict[segment_pair][current_segment_pair_sorted] = [anchorpoint]
                    else:
                        other_segment_pairs_dict[segment_pair][current_segment_pair_sorted].append(anchorpoint)

    clean_segment_pair_dictionary = copy.deepcopy(other_segment_pairs_dict)
    for segment_pair in other_segment_pairs_dict: # For each segment_pair that shares anchor points with other segments
        anchorpoints_current_segment_pair = anchorpoints_per_segment_pair[segment_pair]
        for other_segment_pair in other_segment_pairs_dict[segment_pair]:
            anchorpoints_other_segment_pair = anchorpoints_per_segment_pair[other_segment_pair]
            shortest_anchorpoints_list = min([anchorpoints_current_segment_pair, anchorpoints_other_segment_pair], key=len)
            intersection_current_segment_other_segment = list(set(anchorpoints_current_segment_pair) & set(anchorpoints_other_segment_pair))
            
            # If one list is a subset of the other list, let's remove the short list because it's redundant 
            # If the shorter list shares at least 1/3 of the anchorpoints with the longer list, it is removed (still considered redundant) 
            if set(intersection_current_segment_other_segment) == set(shortest_anchorpoints_list) or len(intersection_current_segment_other_segment) >= 1/3 * len(shortest_anchorpoints_list):   
                if shortest_anchorpoints_list == anchorpoints_current_segment_pair: # the shortest is the current segment pair
                    if segment_pair in clean_segment_pair_dictionary:
                        del clean_segment_pair_dictionary[segment_pair]
                elif shortest_anchorpoints_list == anchorpoints_other_segment_pair: # the shortest is the other segment pair
                    if other_segment_pair in clean_segment_pair_dictionary:
                        del clean_segment_pair_dictionary[other_segment_pair]


    # Getting the Ks lists of the non redundant segment pairs
    nonred_segment_pairs_segment_based = {}
    nonred_segment_pairs_segment_based_no_outliers = {}

    for segment_pair in clean_segment_pair_dictionary:
        anchorpoints_list_current_segment = anchorpoints_per_segment_pair[segment_pair]
        anchorpoints_ks_list_current_segment = []
        for anchorpoint in anchorpoints_list_current_segment:
                ks = ks_anchors[anchorpoint]
                if ks != "" and 0.05 <= float(ks) <= max_ks_para:
                    anchorpoints_ks_list_current_segment.append(float(ks))
        if len(anchorpoints_ks_list_current_segment) > 0: # if there is at least one Ks in the list, add the segment pair to the dictionaries
            # Updating the dictionary that takes into account also outliers
            nonred_segment_pairs_segment_based[segment_pair] = anchorpoints_ks_list_current_segment
            # Updating the dictionary that filters away the outliers by using MAD (median absolute deviation)
            if len(anchorpoints_ks_list_current_segment) <= 5: # if up to 5 Ks, too short for reliable median and MAD statistic
                nonred_segment_pairs_segment_based_no_outliers[segment_pair] = anchorpoints_ks_list_current_segment
            else: # if 6 or more: remove outliers
                median_ks = median(anchorpoints_ks_list_current_segment)
                mad = stats.median_absolute_deviation(anchorpoints_ks_list_current_segment)
                lower_bound, upper_bound = median_ks - mad, median_ks + mad
                ks_list_without_outliers = []
                for ks in anchorpoints_ks_list_current_segment:
                    if lower_bound <= ks <= upper_bound:
                        ks_list_without_outliers.append(ks)
                if len(ks_list_without_outliers) > 0:
                    nonred_segment_pairs_segment_based_no_outliers[segment_pair] = ks_list_without_outliers

    # -----------------------------------------------------------------------------

    # Taking median of each non-redundant Ks list per segment pair

    # There are two options for the Ks list source from which the medians are computed: with or without outliers
    ###chosen_segment_ks_list = nonred_segment_pairs_segment_based (OPTION ONE)
    chosen_segment_ks_list = nonred_segment_pairs_segment_based_no_outliers  # (OPTION TWO)

    all_segment_pairs_ks_median, all_segment_pairs_ks_median_dict = [], {}
    for segment_pair in chosen_segment_ks_list:
        median_ks = median(chosen_segment_ks_list[segment_pair])
        all_segment_pairs_ks_median.append([segment_pair, median_ks])
        all_segment_pairs_ks_median_dict[segment_pair] = median_ks

    # Preparing the list of medians to be given as input to the clustering function
    all_medians_list = []
    for segment_median in all_segment_pairs_ks_median:
        median_value = segment_median[1]
        all_medians_list.append(median_value)
    all_medians_list = array(all_medians_list) # must be a numpy array

    # -----------------------------------------------------------------------------

    # Choosing how many clusters to use, namely as many as the number of WGD events that explains the maximum multiplicon level
    # E.g. if the highest multiplicon level is 8, it is explained by 3 WGDs, then we set the number of clusters as 3
    # We don't consider the combination with WGTs
    num_wgd_per_level = {2: 1, 3: 2, 4: 2, 5: 3, 6: 3, 7: 3, 8: 3, 9: 4, 10: 4, 11: 4, 12: 4, 13: 4, 14: 4, 15: 4, 16: 4, 17: 4, 18: 4, 19: 4}

    if max_level in num_wgd_per_level:
        n_clusters = num_wgd_per_level[max_level]
    else: 
        n_clusters = 4

    logging.info(f"Highest multiplicon level in i-ADHoRe output files: {max_level}")
    logging.info(f"Number of WGDs inferred to explain the highest level: {n_clusters}")
    logging.info("")

    # -----------------------------------------------------------------------------

    # FIRST round of clustering (with all medians)

    # Clustering with GaussianMixtureModel, k-means or lognormal mixture modeling
    clustering_method = "Gaussian mixture modeling" # other options: "k-means" or "lognormal mixture modeling" 
    logging.info(f"Performing a first Ks clustering round with {n_clusters} clusters using {clustering_method}")
    if clustering_method == "Gaussian mixture modeling":
        gmm_clustered_medians = fcCluster.gmm(all_medians_list, n_clusters, max_EM_iterations, num_EM_initializations)
    elif clustering_method == "k-means":
        kmeans_clustered_medians = fcCluster.kmeans(all_medians_list, n_clusters)
    #elif clustering_method == "lognormal mixture modeling":
        # Temporary not available because pomegranate is not on midas:
        #lognormal_clustered_medians = fcCluster.lognormalmm(all_medians_list, n_clusters)

    # TODO: Decide whether to ignore Ks outliers in the segment pair Ks lists or decide instead to consider the complete Ks lists:
    ###chosen_nonred_segment_pair_ks_list = nonred_segment_pairs_segment_based
    chosen_nonred_segment_pair_ks_list = nonred_segment_pairs_segment_based_no_outliers

    # Get the clusters of medians and the resulting clusters of Ks
    cluster_of_ks, medians_per_cluster, segments_medians_per_cluster, cluster_color_letter_list = fcCluster.get_clusters_of_ks(gmm_clustered_medians, all_medians_list, all_segment_pairs_ks_median, chosen_nonred_segment_pair_ks_list, "first")


    logging.info(f"Saving the distribution of clustered medians [{subfolder}/{fcCluster._ANCHOR_CLUSTERS_MEDIANS.format(species)}]")
    # Plot the clusters of medians according to the used method
    fcCluster.plot_clusters_of_medians(medians_per_cluster, cluster_color_letter_list, x_max_lim, bin_list, species, latin_name, output)


    # Generate the plot for the mixed distribution with clusters of Ks
    fig_corr_first, ax_corr_first = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_max_lim, "corrected", correction_table_available, plot_correction_arrows)
    fig_uncorr_first, ax_uncorr_first = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_max_lim, "un-corrected", correction_table_available, plot_correction_arrows)

    # Plot the original complete anchor distribution in the background
    fcPlot.plot_histogram_for_anchor_clustering(ax_corr_first, anchor_ks_list, bin_list, y_max_lim)
    fcPlot.plot_histogram_for_anchor_clustering(ax_uncorr_first, anchor_ks_list, bin_list, y_max_lim)

    # Plot the clusters of anchor Ks and on top of them their KDEs
    clusters_sorted_by_median, cluster_color_letter_list = fcCluster.plot_clusters(ax_corr_first, cluster_of_ks, bin_width, max_ks_para, peak_stats, correction_table_available, plot_correction_arrows)
    fcCluster.plot_clusters(ax_uncorr_first, cluster_of_ks, bin_width, max_ks_para, peak_stats, correction_table_available, plot_correction_arrows)

    # Plotting the ortholog peaks coming from the correction, if available
    if correction_table_available:
        logging.info("Plotting divergence lines")
        fcPlot.plot_divergences(correction_table, peak_stats, consensus_peak_for_multiple_outgroups, ax_uncorr_first, ax_corr_first,
                                color_list, plot_correction_arrows)

    # -----------------------------------------------------------------------------

    # Removing the clusters that are poorly populated (Ks content <= 10%), very old (cluster median >= 3 Ks) or quite horizontally spread (IQR > 1.1 Ks)
    logging.info("")
    logging.info(f"Filtering away Ks clusters with unclear signal (poor Ks content, old Ks age or flat peak)...")
    clean_clusters_of_ks = fcCluster.filter_degenerated_clusters(cluster_of_ks, clusters_sorted_by_median, cluster_color_letter_list)

    # -----------------------------------------------------------------------------

    # SECOND round of clustering (only of medians coming from the good clusters)
    updated_n_clusters = len(clean_clusters_of_ks)
    if updated_n_clusters == n_clusters:
        logging.info("All clusters were retained")

        # Saving the first plot with the final name (leaving away "unfiltered")
        logging.info(f"Saving mixed Ks plot with anchor Ks clusters [{fcCluster._ANCHOR_CLUSTERS_FILTERED.format(species)}]")
        fcCluster.save_anchor_cluster_plot(fig_corr_first, fig_uncorr_first, ax_corr_first, ax_uncorr_first, species, latin_names, correction_table_available, cluster_of_ks, output, "second")
        logging.info("")

    elif updated_n_clusters == 0:
        logging.info("No clusters are left after filtering.")
        logging.info("")

    else: # if one ore more clusters were removed
        logging.info(f"Saving mixed Ks plot with unfiltered anchor Ks clusters [{fcCluster._ANCHOR_CLUSTERS_UNFILTERED.format(species)}]")
        fcCluster.save_anchor_cluster_plot(fig_corr_first, fig_uncorr_first, ax_corr_first, ax_uncorr_first, species, latin_names, correction_table_available, cluster_of_ks, output, "first")
        logging.info("")


        clean_medians_list = []
        clean_segment_medians_list = []

        for cluster in cluster_of_ks:
            if cluster in clean_clusters_of_ks: # if it is not too degenerated
                clean_medians_list.extend(medians_per_cluster[cluster])
                clean_segment_medians_list.extend(segments_medians_per_cluster[cluster])
        clean_medians_list = array(clean_medians_list)
        
        # Clustering with GaussianMixtureModel, k-means or lognormal mixture modeling
        logging.info("")
        logging.info(f"Performing a second Ks clustering round with {updated_n_clusters} clusters through {clustering_method} for the remaining Ks values")
        if clustering_method == "Gaussian mixture modeling":
            gmm_clustered_medians2 = fcCluster.gmm(clean_medians_list, updated_n_clusters, max_EM_iterations, num_EM_initializations)
        elif clustering_method == "k-means":
            kmeans_clustered_medians2 = fcCluster.kmeans(clean_medians_list, updated_n_clusters)
        #elif clustering_method == "lognormal mixture modeling":
            # Temporary not available because pomegranate is not on midas:
            #lognormal_clustered_medians2 = fcCluster.lognormalmm(clean_medians_list, updated_n_clusters)

        filtered_cluster_of_ks, __, __, filtered_cluster_color_list = fcCluster.get_clusters_of_ks(gmm_clustered_medians2, clean_medians_list, clean_segment_medians_list, chosen_nonred_segment_pair_ks_list, "second")

        # Generate the plot for the mixed distribution with clusters of Ks
        fig_corr_second, ax_corr_second = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_max_lim, "corrected", correction_table_available, plot_correction_arrows)
        fig_uncorr_second, ax_uncorr_second = fcPlot.generate_mixed_plot_figure(latin_names.get(species), x_max_lim, y_max_lim, "un-corrected", correction_table_available, plot_correction_arrows)
        
        # Plot the original complete anchor distribution in the background
        fcPlot.plot_histogram_for_anchor_clustering(ax_corr_second, anchor_ks_list, bin_list, y_max_lim)
        fcPlot.plot_histogram_for_anchor_clustering(ax_uncorr_second, anchor_ks_list, bin_list, y_max_lim)

        # Plot the clusters of anchor Ks and on top of them their KDEs
        fcCluster.plot_clusters(ax_corr_second, filtered_cluster_of_ks, bin_width, max_ks_para, peak_stats, correction_table_available, plot_correction_arrows)
        fcCluster.plot_clusters(ax_uncorr_second, filtered_cluster_of_ks, bin_width, max_ks_para, peak_stats, correction_table_available, plot_correction_arrows)

        # Plotting the ortholog peaks coming from the correction, if available
        if correction_table_available:
            logging.info("Plotting divergence lines")
            fcPlot.plot_divergences(correction_table, peak_stats, consensus_peak_for_multiple_outgroups, ax_uncorr_second,
                                    ax_corr_second, color_list, plot_correction_arrows)

        logging.info(f"Saving mixed Ks plot with filtered anchor Ks clusters [{fcCluster._ANCHOR_CLUSTERS_FILTERED.format(species)}]")
        logging.info("")
        fcCluster.save_anchor_cluster_plot(fig_corr_second, fig_uncorr_second, ax_corr_second, ax_uncorr_second, species, latin_names, correction_table_available, filtered_cluster_of_ks, output, "second")

    logging.info("All done")

# -----------------------------------------------------------------------------

# ALTERNATIVE STRATEGY THAT IS NOT USED IN THE CURRENT VERSION, but could be explored later on

# CHOOSING HOW MANY CLUSTERS FOR THE k-means: as many as the WGM events that explain the multiplicon levels
# num_wgd_per_level = {2: 1, 3: 2, 4: 2, 5: 3, 6: 3, 7: 3, 8: 3, 9: 4, 10: 4, 11: 4, 12: 4, 13: 4, 14: 4, 15: 4, 16: 4, 17: 4, 18: 4, 19: 4}

# # Bar plot of amount of anchor pairs per level
# anchors_distrib, axes = plt.subplots(nrows=1, ncols=1, figsize=(10,10))
# anchors_distrib.suptitle(f"Anchorpoints per level (tot. {total_anchorpoints_among_levels}) in ${latin_name}$")
# axes.set_xlabel("Multiplicon level")
# axes.set_ylabel("# anchorpoints")
# axes.set_xticks(sorted(list(multipl_per_level.keys()))) # it's the list of levels (e.g. [2,3,4,5,6])
# bar = axes.bar(sorted(list(multipl_per_level.keys())), list_num_anchorpoints_per_level, color="pink", ec="black")
# for rect in bar:
#     height = rect.get_height()
#     percentage = round(height / total_anchorpoints_among_levels * 100, 1)
#     axes.text(rect.get_x() + rect.get_width()/2.0, height, '%.1f' % percentage + "%", ha='center', va='bottom', color="darkred", fontweight="bold")
#     axes.text(rect.get_x() + rect.get_width()/2.0, height, '%d' % int(height), ha='center', va='top', fontweight="bold")
# anchors_distrib.savefig(os.path.join("level_distributions", f"anchorpoints_distribution_{species}"))
# plt.close()

# # STRATEGY 1:
# # LINKING THE MAXIMUM (significant) MULTIPLICON LEVEL TO THE NUMBER OF WGDs THAT CAN EXPLAIN IT [key=level, value=number of clusters]
# bin_percentages, group_percentages = {}, {}
# for rect in bar:
#     current_level = rect.get_x() + rect.get_width()/2.0
#     height = rect.get_height()
#     bin_percentages[current_level] = round(height / total_anchorpoints_among_levels * 100, 2)
#     current_group = num_wgd_per_level[current_level]
#     if current_group not in group_percentages:
#         group_percentages[current_group] = bin_percentages[current_level]
#     else:
#         group_percentages[current_group] = group_percentages[current_group] + bin_percentages[current_level]
# # Filtering away the groups with low percentage content
# group_percentages_filtered = {}
# threshold_percentage, threshold_number_anchors = get_threshold(total_anchorpoints_among_levels)
# print(f"Group threshold for the {total_anchorpoints_among_levels} anchors: {threshold_percentage} percent ({threshold_number_anchors} anchors)")
# for group in group_percentages:
#     if group_percentages[group] >= threshold_percentage:
#         group_percentages_filtered[group] = group_percentages[group]
# number_of_clusters = max(group_percentages_filtered)

# # STRATEGY 2 (not used for now):
# # TAKING INTO ACCOUNT ALSO THE POSSIBILITY OF a combination of WGDs and WGTs TO EXPLAIN THE MAXIMUM LEVEL: there are two options per level
# number_clusters_wgd_wgt = {2: [1], 3: [1, 2], 4: [2], 5: [2, 3], 6: [2, 3], 
#                         7: [3, 2], 8: [3, 2], 9: [2, 4], 10: [3, 4],
#                         11: [3, 4], 12: [3, 4], 13: [4, 3], 14: [4, 3],
#                         15: [4, 3], 16: [4, 3], 17: [3], 18: [3], 19: [3]}
# # Finding the maximum significant level and the two options it is explained by
# for level in sorted(bin_percentages.keys(), reverse=True): # the for loop goes from maximu level to level 2 (reversed count)
#     if bin_percentages[level] >= threshold_percentage: 
#         maximum_signifiant_level = int(level)
#         percentage_max_signif_level = bin_percentages[level]
#         number_of_clusters_list = number_clusters_wgd_wgt[level]
#         more_probable_n_clusters = number_of_clusters_list[0]
#         try:
#             less_probable_n_clusters = number_of_clusters_list[1]
#         except Exception: less_probable_n_clusters = None
#         break
