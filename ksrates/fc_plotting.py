import os
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerBase
from matplotlib.text import Text
from matplotlib.legend import Legend
from matplotlib.colors import to_rgba
import random
import seaborn
import numpy as np
import pandas
from pandas import Series
from scipy import stats
import ksrates.fc_kde_bootstrap as fcPeak
import matplotlib
# Use the Agg (Anti-Grain Geometry) backend to avoid needing a graphical user interface (X11 backend)
matplotlib.use('Agg')

plt.style.use(os.path.join(f"{os.path.dirname(os.path.abspath(__file__))}", "ks.mplstyle"))
# this should be the path to the directory of the current file "fc_plotting.py",
# and ks.mplstyle is supposed to be in the same directory

# A legend handler to add simple left-aligned text strings to the legend, e.g. (sub)titles
class StringLegendHandler(HandlerBase):
    def create_artists(self, legend, text ,xdescent, ydescent, width, height, fontsize, trans):
        txt = Text(-1, 1, text, fontsize='x-large', ha="left", va="baseline")
        txt.set_transform(trans)
        return [txt]

Legend.update_default_handler_map({str: StringLegendHandler()})

# Some constants to use in plotting
ALPHA_PARANOME_HISTOGRAM = 0.75
ALPHA_ANCHOR_HISTOGRAM = 0.75
ALPHA_DIVERGENCE_RECT = 0.3
ALPHA_DIVERGENCE_LINE = 0.6

COLOR_PARANOME_HISTOGRAM = to_rgba("0.79", ALPHA_PARANOME_HISTOGRAM)
COLOR_ANCHOR_HISTOGRAM = to_rgba("0.64", ALPHA_ANCHOR_HISTOGRAM)
COLOR_PARANOME_KDE = "0.6"
COLOR_ANCHOR_KDE = "0.4"

SIZE_CIRCLE_LABEL_WHITE = 220
LINEWIDTH_CIRCLE_LABEL_BORDER = 1

# Fraction of total y range by which to extend the plot below if correction arrows are drawn
NEGATIVE_Y_FRACTION = 0.08


def generate_mixed_plot_figure(species, x_max_lim, y_max_lim, corrected_or_not, correction_table_available, plot_correction_arrows):
    """
    Initializes a figure with a single empty plot for the mixed distribution.

    :param species: name of the species of interest in the current analysis
    :param x_max_lim: upper limit of the Ks range in the x axis of the plot
    :param y_max_lim: upper limit of the Ks range in the y axis of the plot
    :param corrected_or_not: string to specify if the plot will be the corrected or the un-corrected mixed plot ["corrected", "un-corrected"]
    :param correction_table_available: boolean flag to state whether the correction table is available yet or not
    :param plot_correction_arrows: boolean stating whether there will be plotted correction arrows or not
    :return: figure and axis objects
    """
    species_escape_whitespaces = species.replace(' ', '\ ')

    # Increase figure size in height if extra space for the correction arrows is required and
    # any divergence line will be actually plotted (thus, the "and" condition) 
    if plot_correction_arrows and correction_table_available:
        fig, ax = plt.subplots(1, 1, figsize=(14.0, 7.6))
    else:
        fig, ax = plt.subplots(1, 1, figsize=(14.0, 7.0))

    if correction_table_available:
        if corrected_or_not == "corrected":
            fig.suptitle("Rate-adjusted mixed " + "$K_\mathregular{S}$" + f" distribution for ${species_escape_whitespaces}$", y=0.98)
        elif corrected_or_not == "un-corrected":
            fig.suptitle(f"Mixed " + "$K_\mathregular{S}$" + f" distribution for ${species_escape_whitespaces}$", y=0.98)
    else:
        fig.suptitle("$K_\mathregular{S}$" + f" distribution for ${species_escape_whitespaces}$", y=0.98)

    seaborn.despine(offset=10)
    ax.set_xlabel("$K_\mathregular{S}$")
    ax.set_ylabel("Number of retained duplicates")

    ax.set_xlim(0, x_max_lim)
    if isinstance(y_max_lim, float):
        y_max_lim = float(y_max_lim)
        ax.set_ylim(0, y_max_lim)
    plt.setp(ax.yaxis.get_majorticklabels(), rotation=90, verticalalignment='center')
    return fig, ax


def reflect_KDE(ks_list, weight_list, max_ks):
    """
    Reflects the KDE on the left boundary (values from 0 to 0.5) and on the right boundary
    (values from max_ks to max_ks - 0.5). This removes the KDE boundary effect that may produce
    artefact peaks.

    :param ks_list: list of paralog Ks values to be plotted
    :param weight_list: (optional, default None) list containing the weights for Ks values
    :param max_ks: maximum Ks value accepted for this analysis (e.g. all Ks > 5 are ignored)
    :return ks_list_reflected: a list containing the Ks list between 0 and max_ks and its reflected parts
    :return weight_list_reflected: a list containing the weights of the actual list and of its reflected parts
    :return reflection_threshold: range of Ks values (0.5) reflected at the boarders
    """
    ks_array = np.array(ks_list)
    left_boundary_refl_ks, right_boundary_refl_ks = np.array([]), np.array([])
    left_boundary_refl_weights, right_boundary_refl_weights = np.array([]), np.array([])
    # Set the threshold for the reflection (0.50 Ks left and 0.50 Ks right)
    reflection_threshold = 0.5

    for i in range(len(ks_array)):
        if ks_array[i] <= 0 + reflection_threshold: # (0 + 0.5) will reflect to left the data until Ks of 0.5"
            left_boundary_refl_ks = np.append(left_boundary_refl_ks, ks_array[i] * -1) # reflection around 0
            if weight_list is not None:
                left_boundary_refl_weights = np.append(left_boundary_refl_weights, weight_list[i])
        elif ks_array[i] >= max_ks - reflection_threshold: # e.g. max_ks = 5; (5 - 0.5) will reflect to right the data from Ks of 4.5"
            # TODO: here there should be a "if ks_array[i] != max_ks", to avoid to duplicate also the upper limit number, it should be present once
            right_boundary_refl_ks = np.append(right_boundary_refl_ks, ks_array[i] * -1 + 2 * max_ks) # reflection around max_ks
            if weight_list is not None:
                right_boundary_refl_weights = np.append(right_boundary_refl_weights, weight_list[i])

    ks_list_reflected = np.hstack([left_boundary_refl_ks, ks_array, right_boundary_refl_ks])
    if weight_list is None:
        weight_list_reflected = None
    else:
        weight_list_reflected = np.hstack([left_boundary_refl_weights, np.array(weight_list),
                                            right_boundary_refl_weights])
    return ks_list_reflected, weight_list_reflected, reflection_threshold


def plot_histogram(legend_label, axis, ks_list, bin_list, bin_width, max_ks, bw_modifier, weight_list=None,
                   color=COLOR_PARANOME_HISTOGRAM, plot_kde=True, density=False):
    """
    Generates a histogram from a Ks value list, optionally weighted. Can plot a KDE on top of the Ks histogram.
    
    :param legend_label: string value to set as legend label
    :param axis: matplotlib axis on which the distribution will be plotted
    :param ks_list: list of paralog Ks values to be plotted
    :param bin_list: list of the edges of each bin (e.g. [0.0, 0.1, 0.2 ... ]); regulates how many bins are there
                     per tick in the x axis
    :param bin_width: width of histogram bins (needs to match the values in bin_list)
    :param max_ks: maximum Ks value accepted for this analysis (e.g. all Ks > 5 are ignored)
    :param bw_modifier: modifier to reduce the Scott's factor and improve the KDE fitting on paranome or anchor distributions (default 0.4)
    :param weight_list: (optional, default None) list containing the weights for Ks values
    :param color: (optional, default a dark gray) color of histogram bars
    :param plot_kde: (optional, default True) calculate plot a KDE line on top of the histogram
    :param density: boolean flag stating whether to plot a density plot or not (default: False)
    :return hist: the histogram object
    """
    hist = axis.hist(ks_list, bins=bin_list, weights=weight_list, histtype='stepfilled', color=color,
                     label=legend_label, density=density)

    if plot_kde:
        if bin_width != bin_list[1] - bin_list[0]:
            logging.warning(f"plot_histogram: value of bin_width [{bin_width}] does not match "
                            f"bins in bin_list [{bin_list[0]}, {bin_list[1]}, {bin_list[2]}, ...], "
                            f"the KDE line may not be correctly scaled")

        ks_list_reflected, weight_list_reflected, reflection_threshold = reflect_KDE(ks_list, weight_list, max_ks)

        kde = stats.gaussian_kde(ks_list_reflected, bw_method="scott", weights=weight_list_reflected)
        # change bandwidth so that the line follows more tightly the underlying distribution
        kde.set_bandwidth(kde.factor * bw_modifier)

        kde_x = np.linspace(0 - reflection_threshold, max_ks + reflection_threshold, num=512)
        kde_y = kde(kde_x)
        if weight_list is None:
            axis.plot(kde_x, kde_y * len(ks_list_reflected) * bin_width, color=COLOR_ANCHOR_KDE, linewidth=1.3)
        else:
            axis.plot(kde_x, kde_y * len(ks_list_reflected) * bin_width * np.mean(weight_list_reflected),
                      color=COLOR_PARANOME_KDE, linewidth=1.3)
    return hist


def plot_histogram_for_anchor_clustering(axis, anchor_ks_list, bin_list, y_max_lim):
    """
    Plots the histogram distribution of the anchor pair Ks list in the mixed plot showing anchor Ks clusters.

    :param axis: matplotlib axis on which the distribution will be plotted
    :param anchor_ks_list: list of anchor pair Ks values to be plotted
    :param bin_list: list of the edges of each bin (e.g. [0.0, 0.1, 0.2 ... ]); regulates how many bins are there
                     per tick in the x axis 
    """
    hist = axis.hist(anchor_ks_list, bins=bin_list, histtype='stepfilled', color=COLOR_ANCHOR_HISTOGRAM,
                     label="All anchor pairs", zorder=-20)
    if y_max_lim is None:
        set_mixed_plot_height(axis, y_max_lim, hist)


def set_mixed_plot_height(axis, y_max_lim, hist):
    """
    Sets the height of the plot based on the tallest histogram bin (either of the
    paranome distribution or of the anchor pair distribution).
    """
    tallest_bin = max(hist[0])
    axis.set_ylim(0, tallest_bin * 1.25)


def get_bins(max_ks, bin_width):
    """
    Function to set the number and the width of histogram bins in paralog Ks distribution plot.

    :param max_ks: maximum Ks value accepted for this analysis (e.g. all Ks > 5 are ignored)
    :param bin_width: width of histogram bins in paralogs Ks plots
                      (default: 0.1; there are 10 bins per unit in the x axis)
    :return: bin_list, a list of numbers that are the starting and ending points of each histogram bin:
             [0, 0.1, 0.2, 0.3 ... 9.8, 9.9, 10]
    """
    # "max_ks + bin_width" is used to get also the last number in the list (e.g. 10),
    # otherwise bin_list stops at e.g. 9.9
    bin_list = np.arange(0.0, max_ks + bin_width, bin_width)
    return bin_list


def plot_divergences(correction_table, peak_stats, consensus_peak_for_multiple_outgroups, ax_uncorr, ax_corr, color_list,
                     plot_correction_arrows):
    """
    Plots as vertical lines all the speciation events between the species of interest
    and the sister species (one or more).

    :param correction_table: file with correction data (previously generated by correct.py)
    :param consensus_peak_for_multiple_outgroups: choice of how to deal with multiple corrections
           for the same divergence (expected values: either "mean among outgroups" or "best outgroup")
    :param ax_uncorr: axis object where un-corrected divergence lines will be drawn
    :param ax_corr: axis object where rate-corrected divergence lines will be drawn
    :param color_list: list of colors that will be used to distinguish all divergence lines
    :param plot_correction_arrows: boolean to indicate whether divergence correction arrows should be plotted
    """
    n_divergence_lines = len(correction_table)

    if plot_correction_arrows:
        # extend y-axes into negative range to add space below plot for divergence correction arrows
        __, y_max_corr = ax_corr.get_ylim()
        ax_corr.set_ylim(y_max_corr * -NEGATIVE_Y_FRACTION, y_max_corr)
        ax_corr.spines['left'].set_bounds(0, y_max_corr)
        __, y_max_uncorr = ax_uncorr.get_ylim()
        ax_uncorr.set_ylim(y_max_uncorr * -NEGATIVE_Y_FRACTION, y_max_uncorr)
        ax_uncorr.spines['left'].set_bounds(0, y_max_uncorr)

    y_min_corr, y_max_corr = ax_corr.get_ylim()
    __, y_max_uncorr = ax_uncorr.get_ylim()

    arrow_y_gap = y_min_corr / n_divergence_lines

    zorder_ID = -10  # initialize the zorder_ID to take care of the object order along the z-axis (depth)
    for i, row in correction_table.iterrows():
        zorder_ID -= 2  # it is always decreased by 2 because the zorder involves pairs of objects
        
        node = row["Node"]
        divergence_id = str(node)

        lines_per_node = correction_table.loc[correction_table["Node"] == node]
        if isinstance(lines_per_node, Series):
            n_lines_per_node = 1
        else:
            n_lines_per_node = len(lines_per_node)

        species, sister = row['Species'], row['Sister_Species']

        if peak_stats == "mode": # choosing the MODE as peak Ks for the WGD event
            peak, sd = 'Original_Mode', 'Original_Mode_SD'
        elif peak_stats == "median": # choosing the MEDIAN as peak Ks for the WGD event
            peak, sd = 'Original_Median', 'Original_Median_SD'
        original_peak, original_peak_sd = row[peak], row[sd]
        
        peak_mean_out, sd_mean_out = row['Peak_MeanOut'], row['Peak_MeanOut_SD']
        peak_best_out, sd_best_out = row['Peak_BestOut'], row['Peak_BestOut_SD']
        if consensus_peak_for_multiple_outgroups == "mean among outgroups":
            corrected_peak, corrected_peak_sd = peak_mean_out, sd_mean_out
        elif consensus_peak_for_multiple_outgroups == "best outgroup": 
            corrected_peak, corrected_peak_sd = peak_best_out, sd_best_out
        line_color = color_list[node-1]

        sister_escape_whitespaces = sister.replace(' ', '\ ')

        # Figure with uncorrected divergences
        mean_label = str(round(original_peak, 2))
        ortho_label = f"({divergence_id}) ${sister_escape_whitespaces}$ ({mean_label})"
        # transparent rectangle with width mean +/- sd and dashed mean line
        plot_divergence_line(ax_uncorr, original_peak, original_peak_sd, line_color, ortho_label,
                             plot_correction_arrows, zorder_ID=-300)
        # divergence id number in circle on top of mean line (with jitter in y direction to reduce overlap)
        plot_divergence_id(ax_uncorr, divergence_id, original_peak, y_max_uncorr, n_lines_per_node, line_color,
                           zorder_ID)

        # Figure with corrected divergences
        if original_peak < corrected_peak:
            means_label = str(round(original_peak, 2)) + r'$\rightarrow$' + str(round(corrected_peak, 2))
        else:
            means_label = str(round(corrected_peak, 2)) + r'$\leftarrow$' + str(round(original_peak, 2))
        ortho_label = f"({divergence_id}) ${sister_escape_whitespaces}$ ({means_label})"
        # transparent rectangle with width mean +/- sd and dashed mean line
        plot_divergence_line(ax_corr, corrected_peak, corrected_peak_sd, line_color, ortho_label,
                             plot_correction_arrows, zorder_ID=-300)
        # arrow from original uncorrected to corrected divergence mean below plot
        if plot_correction_arrows:
            plot_divergence_correction_arrow(ax_corr, original_peak, corrected_peak, (i + 1) * arrow_y_gap, line_color,
                                             zorder_ID=-300)
        # divergence id number in circle on top of mean line (with jitter in y direction to reduce overlap)
        plot_divergence_id(ax_corr, divergence_id, corrected_peak, y_max_corr, n_lines_per_node, line_color, zorder_ID)


def plot_divergence_line(axis, peak, sd, line_color, label, plot_correction_arrows, zorder_ID, height=0.88):
    """
    Plots a single divergence line in the mixed plot.

    For more info about zorder_ID see docstring of function "plot_divergence_id". In this function it makes sure that
    the newly added divergence line is in the background of those already present in the plot; makes sure also that
    the rectangular shape is always in the background of the line.

    :param axis: axis object where the divergence line will be drawn on
    :param peak: the Ks coordinate of the peak of the ortholog Ks distribution of the two species
    :param sd: standard deviation associated to the peak
    :param line_color: HTML color of the divergence line
    :param label: legend label for the divergence line
    :param plot_correction_arrows: boolean to indicate whether divergence correction arrows are plotted to determine
                                   correct ymin value of divergence lines
    :param zorder_ID: an integer for the object order along the z axis (the depth of the axis)
    :param height: (default 0.80) height of the line, considering the total height of y axis being equal to 1
    """
    if plot_correction_arrows:
        y_min = NEGATIVE_Y_FRACTION / (1 + NEGATIVE_Y_FRACTION)
    else:
        y_min = 0
    # transparent rectangle with width peak +/- sd
    axis.axvspan(xmin=peak - sd, xmax=peak + sd, ymin=y_min, ymax=height, facecolor=line_color,
                 alpha=ALPHA_DIVERGENCE_RECT, capstyle='projecting', joinstyle='bevel', zorder=zorder_ID - 2)
    # dashed line in correspondence of the distribution peak
    axis.axvline(x=peak, ymin=y_min, ymax=height, color=line_color, alpha=ALPHA_DIVERGENCE_LINE,
                 linestyle=(0, (6, 6)), linewidth=1.1, solid_capstyle='butt', solid_joinstyle='miter',
                 zorder=zorder_ID - 1, label=label)
    # additional transparent white dashed line on top of previous to lighten line on plot but not the lines in legend
    axis.axvline(x=peak, ymin=y_min, ymax=height, color='w', alpha=0.3,
                 linestyle=(0, (6, 6)), linewidth=1.1, solid_capstyle='butt', solid_joinstyle='miter',
                 zorder=zorder_ID)


def plot_divergence_correction_arrow(axis, uncorrected_peak, corrected_peak, arrow_y, line_color, zorder_ID):
    """
    Plots an arrow below the plot pointing from the original uncorrected peak to the (rate-)corrected peak
    of a species divergence.

    :param axis: axis object where the arrow will be drawn on
    :param uncorrected_peak: the Ks value of the original uncorrected peak of the ortholog Ks distribution
                             of two species
    :param corrected_peak: the Ks value of the (rate-)corrected peak of the ortholog Ks distribution
                           of two species
    :param arrow_y: the y location of the arrow
    :param line_color: color of the arrow
    :param zorder_ID: an integer for the object order along the z axis (the depth of the axis)
    """
    if abs(corrected_peak - uncorrected_peak) > 0.02:
        if np.sign(corrected_peak - uncorrected_peak) > 0:
            arrow_connectionstyle = "arc3,rad=0.015"
        else:
            arrow_connectionstyle = "arc3,rad=-0.015"
        axis.annotate("", xy=(corrected_peak, arrow_y), xytext=(uncorrected_peak, arrow_y),
                      arrowprops=dict(arrowstyle="->,head_length=0.3,head_width=0.18", shrinkA=0, shrinkB=0,
                                      edgecolor=line_color, alpha=ALPHA_DIVERGENCE_LINE,
                                      connectionstyle=arrow_connectionstyle, linewidth=0.9))
    else:
        axis.annotate("", xy=(corrected_peak, arrow_y), xytext=(uncorrected_peak, arrow_y),
                      arrowprops=dict(arrowstyle="->,head_length=0.08,head_width=0.15", shrinkA=0, shrinkB=0,
                                      edgecolor=line_color, linewidth=0.9))


# TODO: fix zorder_ID documentation to match implementation
def plot_divergence_id(axis, divergence_id, x, y_max, n_ids, circle_color, zorder_ID):
    """
    Plots the given divergence ID number in a circle above the divergence line
    (with jitter in y direction to reduce overlap).

    zorder_ID is used to organize the visualization order of the different number-circle labels along the z-axis (depth).
    zorder_ID is different every time this function is called, because it is externally increased by 2.
    Example: the first time this function is called, zorder_ID is 1, the number has zorder=-1 and the circle has zorder=-2,
    therefore the circle is in the background of the number. The second time this function is called, zorder_ID =3,
    the number zorder=-3 and the circle zorder=-4, therefore again the circle is in the background of the number,
    but also these two new objects are in the background of the first number-circle pair.
    This way, when a new number/circle pair is added it leis in the background of all the already existing labels.

    :param axis: axis object where the divergence line will be drawn on
    :param divergence_id: identification number for the speciation event in the species of interest evolutionary history;
        the older the speciation, the higher the number; starts from 1
    :param x: x-coordinate of the corrected peak
    :param y_max: height of the y axis
    :param n_ids: total number of speciation events involving the species of interest in the current divergence node in the tree
    :param circle_color: color of the circle around the ID label; it's the same color of the correspondent divergence line
    :param zorder_ID: an integer for the object order along the z axis (the depth of the axis)
    """
    # n_ids is used to spread divergence ID circles proportionally to the amount of lines expected in the same Ks range,
    # so to reduce the chance that circles overlap.
    y_jitter_max = 0.008
    if n_ids > 1:
        y_jitter = random.uniform(-y_jitter_max * n_ids, y_jitter_max * n_ids)
    else:
        y_jitter = 0

    # white circle with colored border and number inside
    y = y_max * (0.922 + y_jitter + (n_ids - 1) * y_jitter_max)
    axis.scatter(x=x, y=y, s=SIZE_CIRCLE_LABEL_WHITE, facecolor='w', edgecolor=circle_color,
                 linewidth=LINEWIDTH_CIRCLE_LABEL_BORDER, zorder=zorder_ID - 1)
    y = y_max * (0.920 + y_jitter + (n_ids - 1) * y_jitter_max)
    axis.text(x=x, y=y, s=divergence_id, fontsize=10, horizontalalignment='center', verticalalignment='center',
              clip_on=True, zorder=zorder_ID, color='k')


def define_legend_size(axis):
    """
    Sets the legend box width for the mixed distribution plots.
    The legend box width is increased in case of very long name labels.

    :param axis: matplotlib axis
    :return: the size of the legend box in tuple format "(x, y, width, height)"
    """
    __, labels = axis.get_legend_handles_labels()
    final_legend_size = (1.01, 0.0, 0.7, 1)
    # the third number is the legend width, initialized as 0.7, which is the
    # basic length to assure a symmetric title in the final figure
    for label in labels:
        label = label.replace("\\", "")
        # Get a legend width that can host the longest latin name of the list 
        if len(label) <= 59:
            current_legend_width = (len(label) + 17) / 100
        elif 60 <= len(label) <= 65:
            current_legend_width = (len(label) + 23) / 100
        elif 66 <= len(label) <= 70:
            current_legend_width = (len(label) + 25) / 100
        elif len(label) >= 71:
            current_legend_width = (len(label) + 30) / 100

        if current_legend_width > final_legend_size[2]:
            final_legend_size = (1.01, 0.0, current_legend_width, 1)
    return final_legend_size


def create_legend(axis, paranome, colinearity, legend_size):
    """
    Places the legend elements associated to the histograms at the beginning of the legend,\\
    while by default they are placed at the end.

    :param axis: matplotlib axis object from which the legend is taken
    :param paranome: the config file field that states if the whole-paranome has to be plotted [yes/no]
    :param colinearity: the config file field that states if the anchor pairs have to be plotted [yes/no]
    :param legend_size: size of the legend box as a tuple of format "(x, y, width, height)"
    :return: the updated legend object
    """
    handles, labels = axis.get_legend_handles_labels()
    sorted_handles, sorted_labels = handles.copy(), labels.copy()
    paranome_rect = Patch(facecolor=COLOR_PARANOME_HISTOGRAM, edgecolor="w")
    anchors_rect = Patch(facecolor=COLOR_ANCHOR_HISTOGRAM, edgecolor="w")
    # empty patch used as spacer between histograms and divergence line legend entries
    empty_rect = Patch(fill=False, edgecolor='none', visible=False)

    if paranome and not colinearity:
        sorted_handles = [paranome_rect, empty_rect, "Divergence with:"] + sorted_handles[:-1]
        sorted_labels = [sorted_labels[-1]] + ["", ""] + sorted_labels[:-1]
    elif not paranome and colinearity:
        sorted_handles = [anchors_rect, empty_rect, "Divergence with:"] + sorted_handles[:-1]
        sorted_labels = [sorted_labels[-1]] + ["", ""] + sorted_labels[:-1]
    elif paranome and colinearity:
        sorted_handles = [paranome_rect, anchors_rect, empty_rect, "Divergence with:"] + sorted_handles[:-2]
        sorted_labels = sorted_labels[-2:] + ["", ""] + sorted_labels[:-2]

    lgd = axis.legend(sorted_handles, sorted_labels, handlelength=1.5, mode="expand", loc="upper left",
                      bbox_to_anchor=legend_size)
    return lgd


def save_mixed_plot(fig_corr, fig_uncorr, ax_corr, ax_uncorr, species, paranome, colinearity):
    """
    This function must be called to save the mixed distribution figure in order to adjust the figure layout:
    the plot area is shrunk to the left and some reasonable space is left on the right side for the legend.

    :param fig_corr: figure object of the corrected mixed distribution
    :param fig_uncorr: figure object of the un-corrected mixed distribution
    :param ax_corr: axis object of the corrected mixed distribution
    :param ax_uncorr: axis object of the un-corrected mixed distribution
    :param species: species of interest
    :param paranome: the config file field that states if the whole-paranome has to be plotted [True/False]
    :param colinearity: the config file field that states if the anchor pairs have to be plotted [True/False]
    """
    legend_size = define_legend_size(ax_corr)

    chart_box = ax_uncorr.get_position()
    # For the un-corrected plot:
    ax_uncorr.set_position([chart_box.x0, chart_box.y0, chart_box.width*0.65, chart_box.height])
    lgd = create_legend(ax_uncorr, paranome, colinearity, legend_size)
    fig_uncorr.savefig(os.path.join("correction_analysis", f"{species}", f"mixed_{species}_uncorrected.pdf"),
                       bbox_extra_artists=(lgd, fig_uncorr._suptitle), bbox_inches="tight", transparent=True, format="pdf")
    # Same thing for the corrected plot:
    ax_corr.set_position([chart_box.x0, chart_box.y0, chart_box.width*0.65, chart_box.height])
    lgd = create_legend(ax_corr, paranome, colinearity, legend_size)
    fig_corr.savefig(os.path.join("correction_analysis", f"{species}", f"mixed_{species}_corrected.pdf"),
                     bbox_extra_artists=(lgd, fig_corr._suptitle), bbox_inches="tight", transparent=True, format="pdf")


def generate_orthologs_figure(species, sister_species, outgroup_species, x_lim):
    """
    Initializes a figure with 2x3 empty subplots where ortholog distributions will be plotted:
    - three columns for the three ortholog Ks distributions between species-sister, species-outspecies and sister-outspecies.
    - two rows: the upper one shows the first 20 KDE lines drawn on the histogram distribution, while the lower one shows the 
      KDE line obtained with the thumb-up Scott's rule and the KDE line obtained by multiplying Scott's factor by 0.7.

    :param species: name of the species of interest
    :param sister_species: name of the sister species of the species of interest
    :param outgroup_species: name of the outspecies for the species of interest and its sister
    :param x_lim: upper limit of the Ks range in the x axis of the subplots
    :return: figure, axes
    """
    species_escape_whitespaces = species.replace(' ', '\ ')
    sister_escape_whitespaces = sister_species.replace(' ', '\ ')
    outspecies_escape_whitespaces = outgroup_species.replace(' ', '\ ')

    fig, axes = plt.subplots(2, 3, sharex=True, sharey="row", figsize=(20.0, 14.0))
    fig.suptitle(f"Ortholog distributions\n\n(${species_escape_whitespaces}$, ${sister_escape_whitespaces}$, "
                 f"outgroup = ${outspecies_escape_whitespaces}$)", fontsize=19, y=0.99)
    axes[0, 0].set_xlim(0, x_lim)
    axes[0, 0].set_title(f"${species_escape_whitespaces}$ — ${sister_escape_whitespaces}$", fontsize="large")
    axes[0, 0].set_ylabel("Number of retained orthologs")

    axes[0, 1].set_xlim(0, x_lim)
    axes[0, 1].tick_params("y", reset=True)
    axes[0, 1].set_title(f"${species_escape_whitespaces}$ — ${outspecies_escape_whitespaces}$", fontsize="large")

    axes[0, 2].set_xlim(0, x_lim)
    axes[0, 2].tick_params("y", reset=True)
    axes[0, 2].set_title(f"${sister_escape_whitespaces}$ — ${outspecies_escape_whitespaces}$", fontsize="large")

    axes[1, 0].set_xlabel("$K_\mathregular{S}$")
    axes[1, 0].set_ylabel("Number of retained orthologs")
    axes[1, 0].tick_params("y", reset=True)
    axes[1, 1].set_xlabel("$K_\mathregular{S}$")
    axes[1, 1].tick_params("y", reset=True)
    axes[1, 2].set_xlabel("$K_\mathregular{S}$")
    axes[1, 2].tick_params("y", reset=True)
    return fig, axes


def plot_orthologs_histogram_kdes(ks_list, bin_list_ortho, bin_width_ortho, ax_boot, ax_kde, x_max_lim_ortho,
                                  bootstrap_kde):
    """
    Plots in the 6-panel figure the histograms of ortholog Ks distributions and their related KDE lines.
    
    :param ks_list: list of ortholog Ks values to be plotted
    :param bin_list_ortho: list of intervals used to delimit histogram bins (e.g. [0, 0.1, 0.2, 0.3 ... 9.8, 9.9, 10])
    :param bin_width_ortho: bin width in orthologs Ks histograms
    :param ax_boot: coordinate to select the desired plot in the upper row of the figure (e.g. axes[0, 0])
    :param ax_kde: coordinate to select the desired plot in the lower row of the figure (e.g. axes[1, 0])
    :param x_max_lim_ortho: upper limit of the Ks range in the ortholog distribution plots (default: 5)
    :param bootstrap_kde: list of data points of the first 20 KDE lines
    """
    # plotting into first row panels
    ax_boot.hist(ks_list, bins=bin_list_ortho, alpha=0.5)
    for kde_data in bootstrap_kde:  # kde_data = [kde_x,kde_y]
        ax_boot.plot(kde_data[0], kde_data[1])

    # plotting into second row panels
    kde, kde_x, kde_y = fcPeak.compute_kde(ks_list, x_max_lim_ortho, bin_width_ortho)

    # extracting the bandwidths
    used_bw = kde.factor
    scotts_rule_bw = kde.scotts_factor()
    silverman_rule_bw = kde.silverman_factor()

    ax_kde.hist(ks_list, bins=bin_list_ortho, alpha=0.5)
    # KDE computed with the automatic Scott bw
    ax_kde.plot(kde_x, kde_y, label=f"Scott's factor: {round(scotts_rule_bw, 3)}")

    kde_avg, kde_x_avg, kde_y_avg = fcPeak.compute_kde(ks_list, x_max_lim_ortho, bin_width_ortho,
                                                       bandwidth_scaling="bootstrap_average")
    # KDE computed with the average bw
    ax_kde.plot(kde_x_avg, kde_y_avg, label=f"Avg bootstrap factor: {round(scotts_rule_bw * 0.7, 3)}")
    ax_kde.legend(fontsize=11)


def plot_orthologs_peak_lines(ortholog_db, tag, ax_num):
    """
    Plots on an ortholog Ks distribution two vertical lines in correspondence of the mode and median
    estimated through bootstrap, plus two background rectangles to visualize their standard deviation.
    
    :param ortholog_db: database of ortholog Ks distribution peaks
    :param tag: string used to search the species pair in the database
                (format example: "Arabidopsis thaliana_Brassica rapa")
    :param ax_num: coordinate of the current plot in the 6-panel figure (e.g. axes[0, 0] for the upper-left)
    """
    peak = ortholog_db.at[tag, 'Ortholog_Mode']
    sd_peak = ortholog_db.at[tag, 'Ortholog_Mode_SD']
    median_value = ortholog_db.at[tag, 'Ortholog_Median']
    sd_median = ortholog_db.at[tag, 'Ortholog_Median_SD']
    # Using rectangle with dashed line:
    plot_divergence_line(ax_num, peak, sd_peak, "k", f"Mean mode: {round(peak, 2)} ± {round(sd_peak, 2)}",
                         plot_correction_arrows=False, zorder_ID=-1, height=1)
    plot_divergence_line(ax_num, median_value, sd_median, "navy",
                         f"Mean median: {round(median_value,2)} ± {round(sd_median,2)}", plot_correction_arrows=False,
                         zorder_ID=-2, height=1)
    ax_num.legend(fontsize=11)
