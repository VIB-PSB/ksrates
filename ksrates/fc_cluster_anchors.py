import os
from numpy import arange, percentile, where, linspace, log
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import logging
import seaborn
from statistics import median
from math import ceil
from scipy import stats
from sklearn.cluster import KMeans
# NOT USED, THE FUNCTION IS COMMENTED from scipy.cluster.vq import kmeans2
import ksrates.fc_kde_bootstrap as fcPeak
import ksrates.fc_plotting as fcPlot
from ksrates.fc_plotting import NEGATIVE_Y_FRACTION
from sklearn.mixture import GaussianMixture      
from matplotlib.legend import Legend
from matplotlib.colors import to_rgba
import matplotlib
# Use the Agg (Anti-Grain Geometry) backend to avoid needing a graphical user interface (X11 backend)
matplotlib.use('Agg')

Legend.update_default_handler_map({str: fcPlot.StringLegendHandler()})

plt.style.use(os.path.join(f"{os.path.dirname(os.path.abspath(__file__))}", "ks.mplstyle"))
# this should be the path to the directory of the current file "fc_plotting.py",
# and ks.mplstyle is supposed to be in the same directory

# Some constants
ALPHA_ANCHOR_CLUSTERS = 0.4
subfolder = "paralogs_analyses"

def parse_segments_file(path_segments_txt):
    """
    Gets from segments.txt all the segments associated to each multiplicon.

    :param path_segments_txt: path to the i-ADHoRe output file "segments.txt"
    :return segments_per_multip: dictionary assigning to each multiplicon all its segments
    """
    with open(path_segments_txt, "r+") as segments_file:
        segments_per_multip = {}
        for line in segments_file:
            line = line.rstrip().split("\t")
            if line and line[0].isdigit():
                multipl_id = int(line[1])
                if multipl_id in segments_per_multip:
                    segments_per_multip[multipl_id].append(int(line[0]))
                else:
                    segments_per_multip[multipl_id] = [int(line[0])]
    return segments_per_multip


def parse_list_elements(path_list_elements_txt):
    """
    Gets from list_elements.txt lists all the genes in the current multiplicon (either anchors or not) and associates with their segment(s).

    :param path_list_elements_txt: path to the i-ADHoRe output file 
    :return segments_from_gene: dictionary assigning to each gene the segment(s) where it lies on 
    """
    with open(path_list_elements_txt, "r+") as list_elements_file:
        segments_from_gene = {}
        for line in list_elements_file:
            line = line.rstrip().split("\t")
            if line and line[0].isdigit():
                gene_id = line[2]
                if gene_id in segments_from_gene:
                    segments_from_gene[gene_id].append(int(line[1]))
                else:
                    segments_from_gene[gene_id] = [int(line[1])]
    return segments_from_gene   


def parse_ks_anchors_tsv_file(path_ks_anchor_file):
    """
    Gets the Ks values in ks_anchors.tsv file.
    
    :param path_ks_anchor_file: path to the wgd anchor Ks output file (format: ".ks_anchors.tsv") 
    :return ks_anchors: dictionary assigning to each anchor pair its Ks value
    """
    with open(path_ks_anchor_file, "r+") as ks_anchors_file:
        ks_anchors = {} 
        next(ks_anchors_file) # to skip the headers in the first row
        for line in ks_anchors_file:
            line = line.rstrip().split("\t")
            anchor1 = line[11]
            anchor2 = line[12]
            anchor_pair_sorted = tuple(sorted((anchor1, anchor2)))
            if anchor_pair_sorted not in ks_anchors:
                ks_anchors[anchor_pair_sorted] = line[8]
    return ks_anchors


def parse_multiplicons_file(path_multiplicons_txt):
    """
    Gets from multiplicons.txt info about multiplicon levels.

    :param path_multiplicons_txt: path to the i-ADHoRe output file "multiplicons.txt"
    :return multipl_per_level: dictionary assigning to each level all the multiplicons with that level
    :return level_of_each_multipl: dictionary assigning to each multiplicon its level
    :return max_level: maximum multiplicon level found 
    :return level_list: list containing the level of each multiplicon
    :return level_list_filtered: list containing te level of long multiplicons (column 8 of input file)
    """
    with open(path_multiplicons_txt, "r+") as multiplicons_file:
        multipl_per_level = {}
        level_of_each_multipl = {} 
        level_list, level_list_filtered = [], []
        multipl_length, multipl_length_big = {}, {}
        # Get multiplicon levels and profile lengths
        next(multiplicons_file)
        for line in multiplicons_file:
            line = line.rstrip().split("\t")
            multipl_id = int(line[0])
            level = int(line[6])
            level_of_each_multipl[multipl_id] = level
            level_list.append(level)

            if level not in multipl_per_level.keys():
                multipl_per_level[level] = [multipl_id]
            else:
                multipl_per_level[level].append(multipl_id)

            length = int(line[8])
            if level in multipl_length.keys():
                multipl_length[level].append(length)
            else:
                multipl_length[level] = [length]

            if length >= 10:
                level_list_filtered.append(level)
                if level in multipl_length_big.keys():
                    multipl_length_big[level].append(length)
                else:
                    multipl_length_big[level] = [length]
        max_level = max(multipl_per_level.keys())
    return multipl_per_level, level_of_each_multipl, max_level, level_list, level_list_filtered


def parse_multiplicon_pairs_file(path_multiplicon_pair_txt, level_of_each_multipl):
    """
    Gets from multiplicon_pairs.txt info about anchor pairs (but not anchorpoints)

    :param path_multiplicon_pair_txt: path to the i-ADHoRe output file "multiplicon_pairs.txt"
    :param level_of_each_multipl: dictionary assigning to each multiplicon its level
    :return anchors_per_multipl: dictionary assigning to each multiplicon all its anchor pairs 
    :return levels_of_each_anchor: dictionary assigning to each anchor pairs the level of the multiplicon which it lies on
    """
    with open(path_multiplicon_pair_txt, "r+") as multiplicon_pairs_file:
        anchors_per_multipl = {}
        levels_of_each_anchor = {} # key = anchor pair, value = all levels in which it is found; let's plot it in the highest level bar.
        # Get from multiplicon_pairs.txt all the anchor pairs found in each multiplicon (complete set of anchors, right?)
        tot_anchor_pair_in_file = 0 # Same length as the anchor pair list in multiplicon_pairs.txt
        # Find unique anchors in multiplicons.txt and the levels in which they are present
        next(multiplicon_pairs_file)
        for line in multiplicon_pairs_file:
            line = line.rstrip().split("\t")
            multipl_id = int(line[1])
            # Sorting the anchor pair so that even if it is present as reversed repetition it is still counted as the same pair
            anchor_pair_sorted = tuple(sorted((line[2], line[3])))
            tot_anchor_pair_in_file += 1

            if multipl_id in anchors_per_multipl:
                anchors_per_multipl[multipl_id].append(anchor_pair_sorted)
            else:
                anchors_per_multipl[multipl_id] = [anchor_pair_sorted]
            
            # For each anchor pair, associates the level(s) in which it is present
            # E.g. levels_of_each_anchor[("GSVIVT01018989001", "GSVIVT01009691001")] = [2, 3, 4]
            # In multiplicons_pairs.txt, same anchor pairs are present more than once written in reversed order: in one line A-B and in another line B-A,
            # But this are repetitions of the same anchor pair that must be ignored
            level_of_current_multipl = level_of_each_multipl[multipl_id]
            if anchor_pair_sorted not in levels_of_each_anchor:
                levels_of_each_anchor[anchor_pair_sorted] = [level_of_current_multipl]
            else:
                # Let's not add the same level more than ones
                if level_of_current_multipl not in levels_of_each_anchor[anchor_pair_sorted]:
                    levels_of_each_anchor[anchor_pair_sorted].append(level_of_current_multipl)
            # Note: length of levels_of_each_anchor is the number of unique anchor pairs from the file (because of ignoring reversed repetitions)
    return anchors_per_multipl, levels_of_each_anchor


def parse_anchorpoints_file(path_anchorpoints_txt, level_of_each_multipl):
    """
    Gets from anchorpoints.txt info about the anchorpoints (multiplicons, levels)

    :param path_anchorpoints_txt: path to the i-adhore output file "anchorpoints.txt"
    :param level_of_each_multipl: dictionary assigning to each multiplicon its level
    :return anchorpoints_per_multipl: dictionary assigning to each multiplicons its anchorpoints
    :return multipl_per_anchorpoint: dictionary assigning to each anchorpoints the multiplicon(s) it lies in
    :return levels_of_anchorpoints: dictionary assigning to each anchorpoint its level
    """
    with open(path_anchorpoints_txt, "r+") as anchorpoints_file:
        anchorpoints_per_multipl = {}
        multipl_per_anchorpoint = {}

        next(anchorpoints_file)
        anchorpoints_lines = anchorpoints_file.readlines()
        levels_of_anchorpoints = {}
        for line in anchorpoints_lines:
            line = line.rstrip().split("\t")
            anchorpoint1, anchorpoint2 = line[3], line[4]
            sorted_anchorpoints = tuple(sorted([anchorpoint1, anchorpoint2]))
            multipl_id = int(line[1])
            if multipl_id not in anchorpoints_per_multipl:
                anchorpoints_per_multipl[multipl_id] = [sorted_anchorpoints]
            else:
                anchorpoints_per_multipl[multipl_id].append(sorted_anchorpoints)

            if sorted_anchorpoints not in multipl_per_anchorpoint:
                multipl_per_anchorpoint[sorted_anchorpoints] = [multipl_id]
            else:
                multipl_per_anchorpoint[sorted_anchorpoints].append(multipl_id)

            multipl_level = level_of_each_multipl[multipl_id]
            if sorted_anchorpoints not in levels_of_anchorpoints:
                levels_of_anchorpoints[sorted_anchorpoints] = [multipl_level]
            else:
                levels_of_anchorpoints[sorted_anchorpoints].append(multipl_level)
    return anchorpoints_per_multipl, multipl_per_anchorpoint, levels_of_anchorpoints


def get_threshold(tot_anchors):
    """
    Note: currently not in use, left for future improvements.
    Gets a threshold for considering a level significantly populated of multiplicons.
    The threshold depends on the total amount of anchorpoints and became stricter with
    increasing the number of anchorpoints.

    :param tot_anchors: total amount of anchorpoints
    :return threshold_percentage: minimum amount of multiplicons required (as percentage)
    :return threshold_number_anchors: minimum amount of multiplicons required (as absolute number)
    """
    if tot_anchors >= 4000:
        threshold_percentage = 0.5 # 0.5% of anchors at least in the bin group to be considered
        threshold_number_anchors = ceil(tot_anchors * 0.005)
        # threshold for 4000 is 20; threshold for 10000 is 50 (too stringent?)
    else:
        threshold_percentage = float(str(20/tot_anchors)[:6]) * 100
        threshold_number_anchors = 20
        # Threshold is the percentage associated to having at least 20 anchors in the bin group
        # Because less than 20 is too few to be considered as a cluster (?)
        # E.g. float(str(20/1500)[:6]) * 100 = 1.3 %
    return threshold_percentage, threshold_number_anchors



def gmm(all_medians_list, n_wgds, max_iter=100, n_init=1):
    """
    Log-transforms the segment pairs medians and fits a Gaussian mixture model.

    :param all_medians_list: list containing all segment pair medians
    :param n_wgds: number of putative WGDs
    :param max_iter: maximum number of iteration for the EM algorithm
    :param n_init: number of initializations
    :return labels: list of cluster IDs to which the medians have been assigned (keeps the order of the input medians)
    """
    # log-transform Ks medians (Ks data tend to have a lognormal shape, but GMM uses normal components)
    all_medians_list_logtransformed = log(all_medians_list)

    gmm = GaussianMixture(n_components=n_wgds, covariance_type="spherical", max_iter=max_iter, n_init=n_init)
    gmm.fit(all_medians_list_logtransformed.reshape(len(all_medians_list_logtransformed), 1))
    labels = gmm.predict(all_medians_list_logtransformed.reshape(len(all_medians_list_logtransformed), 1))
    return labels


def kmeans(all_medians_list, n_wgds):
    """
    Performs k-means to cluster the segment pair medians.

    :param all_medians_list: list containing all segment pair medians
    :param n_wgds: number of putative WGDs 
    :return kmeanModel.labels_: list of cluster IDs to which the medians have been assigned (keeps the order of the input medians)
    """
    kmeanModel = KMeans(n_clusters=n_wgds) # obtained from anchor pairs level distribution with percentage threshold
    kmeanModel.fit(all_medians_list.reshape(len(all_medians_list), 1))
    return kmeanModel.labels_


def get_clusters_of_ks(labels_of_clustered_medians, all_medians_list, all_segment_pairs_ks_median, nonredundant_ks_per_segment_pair, clustering_round):
    """
    Assigns to each cluster the Ks values that belong to that cluster.
    Previously another function (gmm(), Gaussian mixture model) clustered the medians,
    but what has to be actually plotted in the figure are the Ks values,
    not the medians coming from the segment pairs.

    Moreover, each cluster received a color based on its age, which is measured as the median value
    of the Ks list of the cluster. This way, when this plot is generated, 
    the order of colors is always the same (the first/youngest is blue, the next cluster is red and so on).

    :param labels_of_clustered_medians: list containing the cluster ID of each clustered median; 
           the element order is the same as the original median list
    :param all_medians_list: the list of medians
    :param all_segment_pairs_ks_median: a list that pairs each segment pairs with its median
    :param nonredundant_ks_per_segment_pair: dictionary that assigns to each segment pair its Ks list (with or without outliers)
    :param clustering_round: tag to state if it is the first or the second clustering round (allowed values: "first" or "second")
    :return cluster_of_ks: assigns to each cluster ID the Ks list belonging to it
    :return medians_per_cluster: assigns to each cluster ID the list of medians belonging to it
    :return segments_medians_per_cluster: assigns to each cluster ID the list of segment pairs and medians belonging to it
    :return cluster_color_letter_list: assigns to each cluster the right color and letter according to age
    """
    # For each cluster, get the list of medians and of Ks values that fall in it
    cluster_of_ks = {} # {cluster 0: [ks_list], cluster 1: [ks_list], ...}
    medians_per_cluster = {} # {cluster 0: [median_list], cluster 1: [medians_list], ...}
    segments_medians_per_cluster = {}

    for i in range(len(labels_of_clustered_medians)): # the labels of the cluster of each median look like [0 0 0 0 1 1 ...] and they respect the order of input median list
        cluster_current_median = labels_of_clustered_medians[i] # get the cluster ID of the current median (e.g. cluster 2)
        
        # If it's the first k-means round, get the list of medians per cluster and their also their segments for the second clustering
        if clustering_round == "first":
            if cluster_current_median not in medians_per_cluster:
                medians_per_cluster[cluster_current_median] = [all_medians_list[i]]
                segments_medians_per_cluster[cluster_current_median] = [all_segment_pairs_ks_median[i]]
            else:
                medians_per_cluster[cluster_current_median].append(all_medians_list[i])
                segments_medians_per_cluster[cluster_current_median].append(all_segment_pairs_ks_median[i])

        # Get the Ks associated to the median / segment pair
        segments_current_median = all_segment_pairs_ks_median[i][0] # segments which the current median comes from 
        if segments_current_median in nonredundant_ks_per_segment_pair:
            nonred_ks_list_current_segments = nonredundant_ks_per_segment_pair[segments_current_median].copy()

            # For each median, put its Ks list in the list of its cluster ID
            if cluster_current_median in cluster_of_ks:
                cluster_of_ks[cluster_current_median].extend(nonred_ks_list_current_segments)
            else:
                cluster_of_ks[cluster_current_median] = nonred_ks_list_current_segments

    _, cluster_color_letter_list = assign_cluster_colors(cluster_of_ks) # each cluster receives its specific color according to age

    return cluster_of_ks, medians_per_cluster, segments_medians_per_cluster, cluster_color_letter_list


def plot_clusters_of_medians(medians_per_cluster, cluster_color_letter_list, x_axis_max_limit_mixed_plot, bin_list, species, latin_name, output):
    """
    Generates a figure showing the clusters of segment pair medians. 

    :param medians_per_cluster: dictionary assigning to each cluster ID the list of medians belonging to it
    :param cluster_color_letter_list: dictionary assigning to each cluster the right color and letter according to age
    :param x_axis_max_limit_mixed_plot: upper range limit for the x axis
    :param bin_list: list of the edges of each bin (e.g. [0.0, 0.1, 0.2 ... ]); regulates how many bins are there per tick in the x axis 
    :param species: informal name of the species of interest
    :param latin_name: latin name of the species of interest
    :param output: name of the directory where the figure is saved in
    """
    # Plotting the clusters of the segment pair medians
    clusters_histo_medians, ax_clusters_medians = plt.subplots(1, 1, figsize=(12.0, 7.0))
    ax_clusters_medians.set_xlim(0, x_axis_max_limit_mixed_plot)
    ax_clusters_medians.set_xlabel("$K_\mathregular{S}$")
    ax_clusters_medians.set_ylabel('Number of segment pair medians')
    seaborn.despine(offset=10)
    for cluster_id in medians_per_cluster:
        ax_clusters_medians.hist(medians_per_cluster[cluster_id], linewidth=2, bins=bin_list, alpha=0.6,
                                 color=cluster_color_letter_list[cluster_id][0], label=f"Cluster {cluster_id}")
    ax_clusters_medians.legend()
    clusters_histo_medians.suptitle(f'Clustering of segment pair medians from ${latin_name}$')
    plt.setp(ax_clusters_medians.yaxis.get_majorticklabels(), rotation=90, verticalalignment='center')
    clusters_histo_medians.savefig(os.path.join("correction_analysis", f"{species}", output, 
            f"anchor_clusters_{species}_medians.pdf"), transparent=True, format="pdf")


def assign_cluster_colors(cluster_of_ks):
    """
    Assigns to each cluster its right color and letter according to age; 
    the proxy for the age is the median value of the Ks list of the cluster.

    :param cluster_of_ks: dictionary that associates to each cluster its anchor Ks list
    :return clusters_sorted_by_median: list of cluster IDs sorted according to their age
    :return cluster_color_letter_list: dictionary assigning to each cluster the right color and letter according to age
    """
    # Sorting the cluster IDs from youngest to oldest according to the median of their Ks list
    median_per_cluster = {}  # assigns to each cluster the median value of its Ks list; e.g. {0: 3.734, 1: 2.898, ... }
    for cluster_id in cluster_of_ks:
        median_per_cluster[cluster_id] = median(cluster_of_ks[cluster_id])
    clusters_sorted_by_median = sorted(median_per_cluster, key=median_per_cluster.get)

    # Assigning the right color to each cluster according to their median / age (youngest is dark blue, oldest is black)
    color_list = [("royalblue", "a"), ("red", "b"), ("green", "c"), ("black", "d"), ("yellow", "e")]

    cluster_color_letter_list = {}  # key=cluster, value=color
    for cluster_id in clusters_sorted_by_median:
        # take index of the current cluster ID in the list of sorted cluster IDs
        cluster_id_index = clusters_sorted_by_median.index(cluster_id)
        cluster_color_letter_list[cluster_id] = color_list[cluster_id_index]
    return clusters_sorted_by_median, cluster_color_letter_list


def plot_clusters(axis, cluster_of_ks, bin_width, max_ks_para, peak_stats, correction_table_available, plot_correction_arrows):
    """
    Plots the anchor Ks clusters, their KDE lines a marker to show their median and . 

    :param axis: axis object where to plot the clusters
    :param cluster_of_ks: dictionary that associates to each cluster its anchor Ks list
    :param bin_width: width of histogram bins in ortholog Ks plots
    :param max_ks_para: maximum paralog Ks value accepted for the analysis
    :param peak_stats: states whether the cluster peak is intended as median or mode
    :param correction_table_available: boolean tag to state if the correction table data are available or not
    :param plot_correction_arrows: boolean tag to state whether to plot the correction arrows at the bottom of the corrected plot (default: False)
    :return clusters_sorted_by_median: list of cluster IDs sorted according to their age
    :return cluster_color_letter_list: dictionary assigning to each cluster the right color and letter according to age
    """
    clusters_sorted_by_median, cluster_color_letter_list = assign_cluster_colors(cluster_of_ks)

    zorder_ID = 50
    for cluster_id in clusters_sorted_by_median:

        median_of_the_cluster = median(cluster_of_ks[cluster_id])
        
        __, kde_x, kde_y = fcPeak.compute_kde(cluster_of_ks[cluster_id], max_ks_para, bin_width)
        mode_of_the_cluster_KDE, __ = fcPeak.get_mode_from_kde(kde_x, kde_y)

        # transparently color the area under the KDE curve
        kde_area_xy = [(kde_x[0], 0), *zip(kde_x, kde_y), (kde_x[-1], 0)]
        polygon_color = to_rgba(cluster_color_letter_list[cluster_id][0], ALPHA_ANCHOR_CLUSTERS)

        if peak_stats == "mode":
            cluster_label = f"Cluster {cluster_color_letter_list[cluster_id][1]} (mode {round(mode_of_the_cluster_KDE, 2)})"
        elif peak_stats == "median":
            cluster_label = f"Cluster {cluster_color_letter_list[cluster_id][1]} (median {round(median_of_the_cluster, 2)})"

        polygon = mpatches.Polygon(kde_area_xy, facecolor=polygon_color, edgecolor='None', label=cluster_label,
                                   zorder=zorder_ID)
        axis.add_patch(polygon)

        # plot the KDE line
        axis.plot(kde_x, kde_y, color=cluster_color_letter_list[cluster_id][0], zorder=zorder_ID + 1)
        
        # plot a marker to highlight the median of each cluster
        plot_cluster_marker(axis, median_of_the_cluster, kde_x, kde_y, cluster_color_letter_list[cluster_id][0],
                                   cluster_color_letter_list[cluster_id][1], zorder_ID + 2, peak_stats, correction_table_available, plot_correction_arrows)

        zorder_ID += 10
    return clusters_sorted_by_median, cluster_color_letter_list


def plot_cluster_marker(axis, median_of_the_cluster, kde_x, kde_y, color, letter, zorder_ID, peak_stats, correction_table_available, plot_correction_arrows):
    """
    Draws a marker on the KDE curve of each cluster to show where the median 
    of the cluster data is. Label the median peak with the letter associated
    to the cluster.
    The function takes as input the median obtained from the Ks data, but then
    needs to compute the median of the KDE points, as the closest number to the
    real median.

    :param axis: axis object where to plot the marker
    :param median_of_the_cluster: median of the cluster obtained from the Ks data
    :param kde_x: x coordinates of the KDE line
    :param kde_y: y coordinates fo the KDE line
    :param color: color of the cluster whose median is calculated
    :param letter: letter of the cluster whose median is calculated
    :param zorder_ID: determines the order of objects along the z-axis (depth of the plot)
    :param peak_stats: states whether the cluster peak is intended as median or mode
    :param correction_table_available: boolean tag to state if the correction table data are available or not
    :param plot_correction_arrows: boolean tag to state whether to plot the correction arrows at the bottom of the corrected plot (default: False)
    """
    if peak_stats == "median":
        # Get the point of the KDE line that is as closest as possible to the median of the raw data
        x_value = min(kde_x, key=lambda x: abs(x-median_of_the_cluster))
        index_median = where(kde_x == x_value)
        y_value = kde_y[index_median][0]

    elif peak_stats == "mode":
        x_value, y_value = fcPeak.get_mode_from_kde(kde_x, kde_y)
    
    # Plot a vertical line from the median x coordinate up until a fraction of y above the KDE curve
    y_height = axis.get_ylim()[1]
    marker_y_height_fraction = 0.15

    if correction_table_available and plot_correction_arrows:
        ymin = NEGATIVE_Y_FRACTION
    else:
        ymin = 0

    axis.axvline(x=x_value, ymin=ymin, ymax=y_value / y_height + marker_y_height_fraction, color=color,
                 linestyle=(0, (3, 3)), linewidth=1.1, solid_capstyle='butt', solid_joinstyle='miter',
                 zorder=zorder_ID + 1)

    # Plot a filled colored circular label with a white letter to flag it
    # white circle as background to prevent transparency
    axis.scatter(x=x_value, y=y_value + (y_height * marker_y_height_fraction), s=285,
                 facecolor="w", edgecolor="w",
                 linewidth=fcPlot.LINEWIDTH_CIRCLE_LABEL_BORDER, zorder=zorder_ID + 2)
    # circle with same color as background anchor histogram on top of previous
    axis.scatter(x=x_value, y=y_value + (y_height * marker_y_height_fraction), s=285,
                 facecolor=fcPlot.COLOR_ANCHOR_HISTOGRAM, edgecolor="w",
                 linewidth=fcPlot.LINEWIDTH_CIRCLE_LABEL_BORDER, zorder=zorder_ID + 3)
    # circle with transparent cluster color on top of previous
    axis.scatter(x=x_value, y=y_value + (y_height * marker_y_height_fraction), s=285,
                 facecolor=to_rgba(color, ALPHA_ANCHOR_CLUSTERS), edgecolor="w",
                 linewidth=fcPlot.LINEWIDTH_CIRCLE_LABEL_BORDER, zorder=zorder_ID + 4)
    # white letter
    axis.text(x=x_value, y=y_value + (y_height * marker_y_height_fraction), s=letter, fontsize=10,
              horizontalalignment='center', verticalalignment='center', clip_on=True, zorder=zorder_ID + 5, color='w')


def filter_degenerated_clusters(cluster_of_ks, clusters_sorted_by_median, cluster_color_letter_list):
    """
    Removes the clusters that are likely to be just noise or that are too degenerated for the visualization.
    Criteria for degenerated clusters are: Ks content less than 10% of the total (poorly populated),
    median greater than 3 (very old signal) and inter-quartile range greater than 1.1 (flat signal).

    :param cluster_of_ks: dictionary that associates to each cluster its anchor Ks list
    :param clusters_sorted_by_median: list of cluster IDs sorted according to their age
    :param cluster_color_letter_list: dictionary assigning to each cluster the right color and letter according to age
    :return clean_clusters_of_ks: clusters_of_ks without the clusters with bad signal 
    """
    tot_ks_in_clusters = 0
    tot_ks_list_in_clusters = []
    for cluster in cluster_of_ks:
        tot_ks_in_clusters += len(cluster_of_ks[cluster])
        tot_ks_list_in_clusters.extend(cluster_of_ks[cluster])

    clean_clusters_of_ks = {}
    for cluster in clusters_sorted_by_median:
        ks_list = cluster_of_ks[cluster]
        # Get median, first quartile range, third quartile and interquartile range
        median_ks = percentile(ks_list, 50, interpolation="midpoint")
        q1 = percentile(ks_list, 25, interpolation="midpoint") 
        q3 = percentile(ks_list, 75, interpolation="midpoint") 
        iqr = q3 - q1

        evaluation = []
        percentage_ks_content = len(ks_list) / tot_ks_in_clusters
        if percentage_ks_content <= 0.10: # Discard clusters that are poorly populated ( less than 10% of the total Ks)
            evaluation.append(f"neglectable Ks content (< {round(10)}% of data)")
        if median_ks >= 3: # Discard clusters whose median is very old
            evaluation.append("high median (> 3 Ks)")
        if iqr >= 1.1: # Discard clusters that tend to be flat, with no clear peak
            evaluation.append("highly spread data (IQR > 1.1)")
        
        # If the clusters has no negative features, let's add it in the "cleaned list", otherwise just explain why is discarded
        cluster_letter = cluster_color_letter_list[cluster][1]
        if len(evaluation) == 0:
            clean_clusters_of_ks[cluster] = ks_list
            logging.info(f"- Cluster {cluster_letter}: Ks content {round(percentage_ks_content * 100, 2)}%, median {round(median_ks, 2)}, IQR {round(iqr, 2)} --> Accepted")
        elif len(evaluation) == 1:
            logging.info(f"- Cluster {cluster_letter}: Ks content {round(percentage_ks_content * 100, 2)}%, median {round(median_ks, 2)}, IQR {round(iqr, 2)} --> Discarded due to {evaluation[0]}")
        elif len(evaluation) == 2:
            logging.info(f"- Cluster {cluster_letter}: Ks content {round(percentage_ks_content * 100, 2)}%, median {round(median_ks, 2)}, IQR {round(iqr, 2)} --> Discarded due to {evaluation[0]} and {evaluation[1]}")   
        elif len(evaluation) == 3:
            logging.info(f"- Cluster {cluster_letter}: Ks content {round(percentage_ks_content * 100, 2)}%, median {round(median_ks, 2)}, IQR {round(iqr, 2)} --> Discarded due to {evaluation[0]}, {evaluation[1]} and {evaluation[2]}")   
    return clean_clusters_of_ks


# TODO: why update the figure title later, set correct in the first place
# --> Not possible since at the very beginning we don't know the number of clusters
def update_figure_title_cluster_anchors(fig, ax, corrected_or_not, species, latin_names, correction_table_available, cluster_of_ks, round_number):
    """
    Updates the title of the figure showing anchor Ks clusters.
    In case the correction data are not available (yet), the figure title 
    will not mention the corrected divergence lines.

    :param fig: mixed plot figure object
    :param ax: axis object in the figure 
    :param corrected_or_not: flag that states whether the figure shows "rate-adjusted" or "un-corrected" divergence lines 
    :param species: informal name of the species of interest
    :param latin_names: dictionary of scientific names
    :param correction_table_available: tag to state if the correction table data are available or not (allowed values: True or False)
    :param cluster_of_ks: dictionary that associates to each cluster its anchor Ks list
    :param round_number: tag to state if it is the first or the second clustering round (allowed values: "first" or "second")
    :return: the figure suptitle object
    """
    latinSpecies = latin_names[species]
    species_escape_whitespace = latinSpecies.replace(' ', '\ ')
    if round_number == "first":
        filtered_or_not = "Non-filtered anchor "
    elif round_number == "second":
        filtered_or_not = "Anchor "

    if not correction_table_available:
        sup = fig.suptitle(filtered_or_not + "$K_\mathregular{S}$ " + f"clusters for ${species_escape_whitespace}$", y=0.98)
    else:
        sup = fig.suptitle(filtered_or_not + "$K_\mathregular{S}$ " + f"clusters with {corrected_or_not} divergences for ${species_escape_whitespace}$", y=0.98)
    return sup


def create_legend(axis, legend_size):
    """
    Places the legend elements associated to the total anchor Ks histogram and to the clusters 
    at the beginning of the legend (by default they are placed at the end of the legend)

    :param axis: the axis object from which the legend is taken
    :param legend_size: size of the legend box as a tuple of format "(x, y, width, height)"
    :return: the updated legend object
    """
    handles, labels = axis.get_legend_handles_labels()
    # empty patch used as spacer between histograms and divergence line legend entries
    empty_rect = mpatches.Patch(fill=False, edgecolor='none', visible=False)
    
    background_distr_label_index = labels.index("All anchor pairs")
    total_anchors_and_clusters_labels = labels[background_distr_label_index:]
    del labels[background_distr_label_index:]
    total_anchors_and_clusters_handles = handles[background_distr_label_index:]
    del handles[background_distr_label_index:]

    sorted_handles = total_anchors_and_clusters_handles + [empty_rect, "Divergence with:"] + handles
    sorted_labels = total_anchors_and_clusters_labels + ["", ""] + labels

    lgd = axis.legend(sorted_handles, sorted_labels, handlelength=1.5, mode="expand", loc="upper left",
                      bbox_to_anchor=legend_size)
    return lgd


def save_anchor_cluster_plot(fig_corr, fig_uncorr, ax_corr, ax_uncorr, species, latin_names,
                             correction_table_available, cluster_of_ks, output, round_number):
    """
    This function must be called to save the figure of anchor clusters in order to adjust the figure layout:
    the plot area is shrunk to the left and some reasonable space is left on the right side for the legend.

    :param fig_corr: figure object of the corrected mixed distribution
    :param fig_uncorr: figure object of the un-corrected mixed distribution
    :param ax_corr: axis object of the corrected mixed distribution
    :param ax_uncorr: axis object of the un-corrected mixed distribution
    :param species: species of interest
    :param latin_names: dictionary of scientific names
    :param correction_table_available: boolean to state if the correction table data are available or not
    :param cluster_of_ks: dictionary that associates to each cluster its anchor Ks list
    :param output: output subfolder name where to store all pictures associated to paralogs analyses
    :param round_number: tag to state if it is the first or the second clustering round
           (allowed values: "first" or "second")
    """
    legend_size = fcPlot.define_legend_size(ax_corr)
    chart_box = ax_uncorr.get_position()

    if correction_table_available:
        ax_corr.set_position([chart_box.x0, chart_box.y0, chart_box.width*0.65, chart_box.height])
        lgd = create_legend(ax_corr, legend_size)
    else:
        ax_corr.set_position([chart_box.x0, chart_box.y0, chart_box.width*0.9, chart_box.height])
        lgd = ax_corr.legend(handlelength=1.5, mode="expand", loc="upper left", bbox_to_anchor=(0.63, 0.0, 0.75, 1))

    update_figure_title_cluster_anchors(fig_corr, ax_corr, "corrected", species, latin_names,     
                                         correction_table_available, cluster_of_ks, round_number)

    if round_number == "first": # unfiltered
        figure_file_path = os.path.join("correction_analysis", f"{species}", output, f"mixed_{species}_anchor_clusters_unfiltered.pdf")
    elif round_number == "second": # filtered, only significant clusters
        figure_file_path = os.path.join("correction_analysis", f"{species}", f"mixed_{species}_anchor_clusters.pdf")

    fig_corr.savefig(figure_file_path, bbox_extra_artists=(ax_corr, lgd, fig_corr._suptitle), bbox_inches="tight",
                     transparent=True, format="pdf")
