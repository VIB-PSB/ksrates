from numpy import log, exp, array, arange, random, linspace, histogram, append, hstack, concatenate, flip, argmax, round, repeat, mean, float64
import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
from scipy.stats import norm, lognorm, expon, gaussian_kde
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy.interpolate import UnivariateSpline
from math import floor, sqrt
import logging
import os
import ksrates.fc_plotting as fcPlot
import ksrates.fc_extract_ks_list as fc_extract_ks_list
from ksrates.fc_cluster_anchors import ALPHA_ANCHOR_CLUSTERS
from ksrates.fc_plotting import NEGATIVE_Y_FRACTION
from ksrates.fc_cluster_anchors import subfolder


# TODO: when the script will be moved to the final directory, do we have to change this path?
plt.style.use(os.path.join(f"{os.path.dirname(os.path.abspath(__file__))}", "ks.mplstyle"))
# this should be the path to the directory of the current file,
# and ks.mplstyle is supposed to be in the same directory

def generate_peak_model_figure(species_escape_whitespace, x_max_lim):
  """
  Generates the figure hosting the plots for the mixture models
  obtained from Ks data and from data plus a random lognormal.

  :param species_escape_whitespace: scientific name of the focal species with space escape
  :param x_max_lim: upper Ks limit in the x-axis range
  :return: figure, axes and figure title
  """
  fig_peaks, [[ax_peaks_ks, ax_peaks_logks], [ax_peaks2_ks, ax_peaks2_logks]] = plt.subplots(nrows=2, ncols=2, figsize=(20, 16), sharey="row")
  sup_peaks = fig_peaks.suptitle(f"Exponential-Lognormal mixture model on ${species_escape_whitespace}$ " +"$K_\mathregular{S}$ " + f"paranome\n\nInitialized from data")
  ax_peaks2_ks.set_xlabel("$K_\mathregular{S}$")
  ax_peaks2_logks.set_xlabel("ln $K_\mathregular{S}$")
  for ax_pair in [[ax_peaks_ks, ax_peaks_logks], [ax_peaks2_ks, ax_peaks2_logks]]:
    ax_pair[0].set_xlim(0, x_max_lim)
    ax_pair[1].set_xlim(-5, 2)
  return fig_peaks, ax_peaks_ks, ax_peaks_logks, ax_peaks2_ks, ax_peaks2_logks, sup_peaks


def generate_random_model_figure(species_escape_whitespace, min_num_comp, max_num_comp, x_max_lim):
  """
  Generates the figure hosting the plots for the mixture models
  obtained from random initialization.

  :param species_escape_whitespace: scientific name of the focal species with space escape
  :param min_num_comp: minimum number of components used in the mixture model
  :param max_num_comp: maximum number of components used in the mixture model
  :param x_max_lim: upper Ks limit in the x-axis range
  :return: figure, axes and figure title
  """
  fig_random, axes_random = plt.subplots(nrows=(max_num_comp-min_num_comp+1), ncols=2, figsize=(20, 8*(max_num_comp-min_num_comp+1)), sharey="row")
  sup_random = fig_random.suptitle(f"Exponential-Lognormal mixture model on ${species_escape_whitespace}$ " + "$K_\mathregular{S}$ " + f"paranome\n\nRandom initialization")
  if len(axes_random.shape) == 1: # there is only one row, because the max_num_comp is 3
    axes_random = [axes_random]   # convert into a list of list just to mimic a multi-row matrix
  axes_random[-1][0].set_xlabel("$K_\mathregular{S}$")
  axes_random[-1][1].set_xlabel("ln $K_\mathregular{S}$")
  for ax_pair in list(axes_random):
    ax_pair[0].set_xlim(0, x_max_lim)
    ax_pair[1].set_xlim(-5, 2)
  return fig_random, axes_random, sup_random


def generate_best_model_figure(latin_species, x_max_lim, y_max_lim, correction_table_available, plot_correction_arrows):
  """
  Generates the figure hosting the best mixture model
  (according to lowest BIC score).

  :param latin_species: scientific name of the focal species
  :param x_max_lim: upper Ks limit in the x-axis range
  :param y_max_lim: upper Ks limit in the y-axis range
  :param correction_table_available: boolean stating whether the correction results are available or not
  :param plot_correction_arrows: boolean stating whether there will be plotted correction arrows or not
  :return figure, axes and figure title
  """
  fig_best_model, ax_best_ks = fcPlot.generate_mixed_plot_figure(latin_species, x_max_lim, y_max_lim, "corrected", correction_table_available, plot_correction_arrows)
  ax_best_ks.set_ylabel("Number of retained duplicates")
  return fig_best_model, ax_best_ks


def init_parameters_randomly(model_id, num_comp, ax_ks, max_ks_mixture):
  """
  Initializes parameters for the given number of components by randomly drawing
  from suitable ranges for Ks values.
  Appends an extra buffer-lognormal component targeted to the right-most Ks area
  to avoid that the other components stretch towards right in the forced attempt to cover that area.
  
  :param model_id: current model identification number
  :param num_comp: number of components (1 exponential distribution + N lognormal components)
  :param ax_ks: axis object showing Ks paranome in background 
  :param max_ks_mixture: upper Ks limit considered when fitting the mixture models 
  :return init_means: list of randomly initialized means for the Gaussian components
  :return init_stdevs: list of randomly initialized standard deviations for the Gaussian components
  :return init_lambd: randomly initialized rate parameter for the exponential component
  :return init_weights: initial component weights (they all have even weights)
  """
  ax_ks.set_title(f"Model {model_id}")
  init_means, init_stdevs, init_weights, init_lambd = [], [], [1/num_comp] * num_comp, 0

  init_lambd = round(random.choice(arange(0.2, 1, 0.1)), 2)
  for n in range(num_comp-2): # the range doesn't consider the exponential and the buffer-gaussian components
    mean = round(random.choice(arange(-0.5, 1, 0.1)), 1)
    sigma = round(random.choice(arange(0.3, 0.9, 0.1)), 1)
    init_means.append(mean)
    init_stdevs.append(sigma)
  
  # Add an extra "buffer" lognormal that is meant to cover the right-most area;
  # e.g. to have a lognormal around 5, the mean of the relative normal must be set
  # around ln(5) = 1.6 Ks.
  init_means.append(log(max_ks_mixture))
  init_stdevs.append(0.3)
  return init_means, init_stdevs, init_lambd, init_weights


def init_parameters_from_data(ax_ks, species, ks_data_log, ks_weights_log, ks_data, ks_weights, species_escape_whitespace, output, max_ks_mixture):
    """
    Initializes the component parameters from data: the exponential rate is guessed from the
    height of the first histogram bin (width 0.1), while the normal means and standard 
    deviations are guessed by detecting peaks in the log-transfromed paranome data.
    If too many lognormals are found, four are randomly drawn and the others ignored (it should not happen often).
    Appends an extra buffer-lognormal component targeted to the 4-5 Ks area
    to avoid that the other components stretch towards right in the forced attempt to cover that area.
    Prints the initial parameters on screen.

    :param ax_ks: axis object showing Ks paranome in background 
    :param species: informal name of the focal species
    :param ks_data_log: log-transformed paranome Ks values
    :param ks_weights_log: weights associated to the log-transformed paranome Ks values
    :param ks_data: paranome Ks values of the focal species
    :param ks_weights: weights associated to the paranome Ks values
    :param species_escape_whitespace: species name escaping the white space
    :param output: output folder
    :param max_ks_mixture: upper Ks limit considered when fitting the mixture models 
    :return init_lambd: rate parameter of the exponential component initialized as the height of the first histogram bin (width 0.1)
    :return init_means: list of means for the Gaussian components initialized by peak detection in log-transformed data
    :return init_stdevs: list of standard deviations for the Gaussian components initialized by peak detection in log-transformed data
    :return init_weights: initial component weights (they all have even weights)
    """
    # NOTE: This histogram has by convention bin_width equal to 0.1,
    # the height of its first bin is used to guess the rate of the exponential component
    bin_list_for_lambda = arange(0, 5.1, 0.1)
    hist_data = histogram(ks_data, weights=ks_weights, bins=bin_list_for_lambda, density=True) 
    init_lambd = round(hist_data[0][0], 2) 

    # Get the spline of the log-transformed Ks paranome
    spl_x, spl_y = get_spline(ks_data_log, ks_weights_log, species, species_escape_whitespace, output)

    # Initialize the lognormals based on peak detection in the spline
    init_means, init_stdevs = find_peak_init_parameters(spl_x, spl_y, species, species_escape_whitespace, output)

    # Add an extra "buffer" lognormal that is meant to cover the area around 4-5 Ks
    reduced_gaussians = False
    if len(init_means) > 4:
        # if too many lognormals were already found, retain four of them by chance
        init_means = list(random.choice(init_means, size=4, replace=False))
        init_stdevs = list(random.choice(init_stdevs, size=4, replace=False))
        reduced_gaussians = True
    init_means.append(log(max_ks_mixture)) # e.g. the means of the normal is at 1.6 Ks to have a lognormal around 5 Ks
    init_stdevs.append(0.3)

    num_comp = len(init_means)+1 # counts also the buffer lognormal
    init_weights = [1/num_comp] * num_comp
    return init_lambd, init_means, init_stdevs, init_weights, reduced_gaussians


def reflect_logT_data(ks_data_log, ks_weights_log, max_logks):
    """
    Reflects log-transformed data around the upper Ks limit by 0.7 Ks.
    The reflection helps in buffering the boundary effect when plotting the KDE.
    It is not performed on the left because the shape of the log-transformed
    generally has a thin and long left tail but is truncated
    at the right boundary because of stopping the original Ks distribution at 5 Ks
    (or at any customized cutoff from the expert configuration file).

    :param ks_data_log: log-transformed paranome Ks values
    :param ks_weights_log: weights associated to the log-transformed paranome Ks values
    :param max_logks: maximum value in the log-transformed paranome
    :return ks_list_reflected: Ks paranome including extra Ks values obtained by reflection
    :return weight_list_reflected: weights of the Ks paranome including extra weights associated to the reflected part
    :return reflection_threshold: Ks range (0.7 Ks) that is reflected beyond the upper boundary
    """
    ks_array = array(ks_data_log)
    right_boundary_refl_ks = array([]) # array for the reflected values around the right boundary
    right_boundary_refl_weights = array([]) # array for the weights of the reflected values around the right boundary
    reflection_threshold = 0.7

    for i in range(len(ks_array)):
      if ks_array[i] >= max_logks - reflection_threshold: # e.g. max_logks = 5; (5 - 0.7) will reflect to right the data from Ks of 4.3"
        if ks_array[i] != max_logks: # avoid repeating the values located precisely at the right boundaries of reflection
          right_boundary_refl_ks = append(right_boundary_refl_ks, ks_array[i] * -1 + 2 * max_logks) # reflection around max_logks
          right_boundary_refl_weights = append(right_boundary_refl_weights, ks_weights_log[i])

    ks_list_reflected = hstack([ks_array, right_boundary_refl_ks])
    weight_list_reflected = hstack([array(ks_weights_log), right_boundary_refl_weights])
    return ks_list_reflected, weight_list_reflected, reflection_threshold


def get_spline(ks_data_log, ks_weights_log, species, species_escape_whitespace, output):
  """
  Obtains the KDE from the log-transformed paranome Ks, then obtains the spline
  of such KDE to smooth out the noise (little bumps) in the curve (the noise tricks 
  the successive peak detection step).
  Generates a figure for comparison purposes between the KDE and the spline.

  :param ks_data_log: log-transformed paranome Ks values
  :param ks_weights_log: weights associated to the log-transformed paranome Ks values
  :param species: informal name of the focal species
  :param species_escape_whitespace: species latin name escaping the white space
  :param output: output folder
  :return spl_x: x-axis coordinates of the spline
  :return spl_y: y-axis coordinates of the spline
  """
  kde_spl, axes_kde_spl = plt.subplots(nrows=1, ncols=3, figsize=(25, 7), sharey=True, sharex=True)
  axes_kde_spl[0].set_xlim(-5,2)
  kde_spl.suptitle("KDE and spline of log-transformed $K_\mathregular{S}$ paranome of " + f"${species_escape_whitespace}$")

  for ax_spl in axes_kde_spl:
    ax_spl.set_xlabel("ln $K_\mathregular{S}$")
  axes_kde_spl[0].set_ylabel("Density of retained duplicates")

  # Get KDE
  max_ks = ks_data_log.max()
  min_ks = ks_data_log.min()
  ks_list_reflected, weight_list_reflected, reflection_threshold = reflect_logT_data(ks_data_log, ks_weights_log, max_ks)

  kde = gaussian_kde(ks_list_reflected, bw_method="scott", weights=weight_list_reflected)
  bw_modifier = 0.4  # change bandwidth so that the line follows more tightly the underlying distribution
  kde.set_bandwidth(kde.factor * bw_modifier)
  # Get the x-coordinates points for plotting the KDE
  # (the linspace is set with 100 points per x axis unit for consistency with later in the code) 
  kde_x = linspace(min_ks-reflection_threshold, max_ks+reflection_threshold, num=int(round((abs(min_ks) + max_ks) * 100 + (reflection_threshold*2*100))))
  kde_y = kde(kde_x)
  for w in [0, 1]:
    axes_kde_spl[w].plot(kde_x, kde_y, color="k", lw=1, label="KDE")
  
  # Get spline 
  spl = UnivariateSpline(kde_x, kde_y) # here it also covers the reflected part
  spl.set_smoothing_factor(0.01)
  for w in [0, 2]:
    axes_kde_spl[w].plot(kde_x, spl(kde_x), 'r', lw=1, label="Spline on KDE")
  # Limit the spline to the range of interest (min Ks value - max Ks value),
  # but make it go a bit beyond the max Ks to include a possible peak placed at the right boundary
  # NOTE: there must be 100 data points per x-axis unit for later being able to compute the peak width (e.g. 500 points within 0-5 Ks range)
  spl_x = linspace(min_ks, max_ks+0.1, num=int(round((abs(min_ks) + (max_ks+0.1)) *100)))
  spl_y = spl(spl_x)

  # Plot histogram of reflected log-transformed paranome and a vertical line on the reflection boundary
  for ax_kde_spl in axes_kde_spl:
    ax_kde_spl.hist(ks_list_reflected, weights=weight_list_reflected, bins=arange(-5, 5, 0.1), color="r", alpha=0.2, density=True)
    ax_kde_spl.axvline(x=max_ks, color="0.5", linestyle="--", lw=1.1, label=f"Reflection boundary")
    ax_kde_spl.legend(loc="upper left")

  axes_kde_spl[0].set_ylim(0, axes_kde_spl[0].get_ylim()[1] * 1.15) # resize to not overlap with legend
  kde_spl.savefig(os.path.join("rate_adjustment", f"{species}", output, f"elmm_{species}_kde_spline.pdf"), bbox_inches="tight")
  plt.close(kde_spl)
  return spl_x, spl_y


def find_peak_init_parameters(spl_x, spl_y, species, species_escape_whitespace, output):
  """
  Finds peaks in the spline of the log-transformed paranome and interprets them as 
  traces of Gaussian distributions coming from WGM events; guesses then their mean
  and standard deviations. More details below.

  First, Scipy is used to detect peaks. The prominence (height) of a peak can be informative of whether
  the peak is just a small noise bump or if it has an actual WGM signal.
  However, due to having possible (close) overlapping peak signals, WGM peak signals do not always have a 
  tall rising-decreasing pattern and the peak prominence can be also quite small.
  This brings the problem on how to distinguish noise from potential good peaks based on prominences.

  An strategy to face this problem is to reflect the distribution around each peak in both directions
  and re-evaluate the prominence of the peak: if at least one reflection gives a significant prominence,
  the peak is retained.
  Note that if a peak is just a small noise bump, it will be also filtered away after the reflections.

  After this filtering, each peak x-coordinate becomes the mean of a Gaussian component and the width
  of the peak becomes its standard deviation. If the width is too large, a default intermediate value
  is set to avoid having too flat initial Gaussians.

  :param spl_x: x-axis coordinates of the spline
  :param spl_y: y-axis coordinates of the spline
  :param species: informal name of the focal species
  :param species_escape_whitespace: species name escaping the white space
  :param output: output folder
  :return init_means: list of randomly initialized means for the lognormal components
  :return init_stdevs: list of randomly initialized standard deviations for the lognormal components
  """
  # Find peaks and their prominences (heights) through Scipy on the spline
  peaks, __ = find_peaks(spl_y)
  prominences = peak_prominences(spl_y, peaks)[0]

  fig_refl_RL, axes_refl = plt.subplots(nrows=len(peaks)+1, ncols=2, figsize=(14, 7*(len(peaks)+1)), sharey=True)
  fig_refl_RL.suptitle("Peak detection in log-transformed $K_\mathregular{S}$ paranome of " + f"${species_escape_whitespace}$")
  for w in range(len(peaks)+1):
    axes_refl[w][0].set_ylabel("Density of retained duplicates")
  for w in [0,1]:
    axes_refl[len(peaks)+1-1][w].set_xlabel("ln $K_\mathregular{S}$")
  # Plotting peaks and prominences found by Scipy
  for ax_refl in axes_refl:
    ax_refl[0].plot(spl_x, spl_y, color="0.5", linewidth=1)
    ax_refl[1].plot(spl_x, spl_y, color="0.5", linewidth=1)
  axes_refl[0,0].scatter(spl_x[peaks], spl_y[peaks], marker="x", c="r", label="peaks") 
  axes_refl[0,0].vlines(x=spl_x[peaks], ymin=spl_y[peaks]-prominences, ymax=spl_y[peaks], color="r", label="prominences")
  axes_refl[0,1].scatter(spl_x[peaks], spl_y[peaks], marker="x", c="r", label="peaks") 
  axes_refl[0,1].vlines(x=spl_x[peaks], ymin=spl_y[peaks]-prominences, ymax=spl_y[peaks], color="r", label="prominences")

  # Filtering away non-significant peaks through reflection
  peak_threshold = 0.06 # arbitrary threshold based on empirical data

  # 1) reflect L<--R
  reflect_prominences_R1, reflect_widths_R1 = array([]), array([])

  for i in range(len(peaks)):
    peak_index = peaks[i]
    # Reflect and plot the spline around the peak
    spl_peak_refl_y = concatenate((flip(spl_y[peak_index+1:]), spl_y[peak_index:]))
    spl_peak_refl_x = concatenate((flip(spl_x[peak_index+1:] * -1 + 2 * spl_x[peak_index]), spl_x[peak_index:]))
    axes_refl[i+1,0].plot(spl_peak_refl_x, spl_peak_refl_y)

    current_peak_index = floor(len(spl_peak_refl_y)/2) 
    # Get the new peak prominence and width after reflection
    peak_prominences_after_reflect = round(peak_prominences(spl_peak_refl_y, [current_peak_index])[0][0], 2)
    peak_width_after_reflect, height, __, __ = peak_widths(spl_peak_refl_y, [current_peak_index], rel_height=0.4)
    if peak_width_after_reflect > 150: # probably the width is horizontally stretched due to an overlapping signal in the reflected spline
      peak_width_after_reflect = 150
    else:
      peak_width_after_reflect = round(peak_width_after_reflect[0], 2)
    reflect_prominences_R1 = append(reflect_prominences_R1, peak_prominences_after_reflect)
    reflect_widths_R1 = append(reflect_widths_R1, peak_width_after_reflect)
    if peak_prominences_after_reflect >= peak_threshold:
      color = "r"
    else:
      color = "0.3"
    # Plot the prominence after reflection
    axes_refl[i+1,0].vlines(x=spl_x[peak_index], ymin=spl_y[peak_index] - peak_prominences_after_reflect, ymax = spl_y[peak_index], color=color, label=peak_prominences_after_reflect)
    if peak_prominences_after_reflect >= peak_threshold:
      # Plot the width after reflection for significant prominences
      axes_refl[i+1,0].hlines(y=height[0], xmin=spl_x[peak_index], xmax=spl_x[peak_index]+(peak_width_after_reflect/2)/100, linestyles="-", color="darkred", lw=1, label=round((peak_width_after_reflect / 2 / 100), 2))
    axes_refl[i+1,0].legend(frameon=False)

  significant_peaks_R1 = reflect_prominences_R1  >= peak_threshold # boolean list stating which prominences are significant

  # 2) reflect L-->R
  reflect_prominences_R2, reflect_widths_R2 = array([]), array([])

  for i in range(len(peaks)):
    peak_index = peaks[i]
    # Reflect and plot the spline around the peak
    spl_peak_refl_y = concatenate((spl_y[:peak_index+1], flip(spl_y[:peak_index])))
    spl_peak_refl_x = concatenate((spl_x[:peak_index+1], flip(spl_x[:peak_index]) * -1 + 2*spl_x[peak_index]))
    axes_refl[i+1,1].plot(spl_peak_refl_x, spl_peak_refl_y)

    current_peak_index = floor(len(spl_peak_refl_y)/2)
    # Get the new peak prominence and width after reflection
    peak_prominences_after_reflect = round(peak_prominences(spl_peak_refl_y, [current_peak_index])[0][0], 2)
    peak_width_after_reflect, height, __, __ = peak_widths(spl_peak_refl_y, [current_peak_index], rel_height=0.4)
    if peak_width_after_reflect > 150: # probably the width is horizontally stretched due to an overlapping signal in the reflected spline
      peak_width_after_reflect = 150
    else:
      peak_width_after_reflect = round(peak_width_after_reflect[0], 2)
    reflect_prominences_R2 = append(reflect_prominences_R2, peak_prominences_after_reflect)
    reflect_widths_R2 = append(reflect_widths_R2, peak_width_after_reflect)
    if peak_prominences_after_reflect >= peak_threshold:
      color = "r"
    else:
      color = "0.3"
    # Plot the prominence after reflection
    axes_refl[i+1,1].vlines(x=spl_x[peak_index], ymin=spl_y[peak_index] - peak_prominences_after_reflect, ymax = spl_y[peak_index], color=color, label=peak_prominences_after_reflect)
    if peak_prominences_after_reflect >= peak_threshold:
      # Plot the width after reflection for significant prominences
      axes_refl[i+1,1].hlines(y=height[0], xmin=spl_x[peak_index], xmax=spl_x[peak_index]+(peak_width_after_reflect/2)/100, linestyles="-", color="darkred", lw=1, label=round((peak_width_after_reflect / 2 / 100), 2))
    axes_refl[i+1,1].legend(frameon=False)

  significant_peaks_R2 = reflect_prominences_R2  >= peak_threshold # boolean list stating which prominences are significant

  axes_refl[0,0].set_title('Reflection L <-- R')
  axes_refl[0,0].legend()
  axes_refl[0,1].set_title('Reflection L --> R')
  axes_refl[0,1].legend()
  original_y_lim = axes_refl[0][0].get_ylim()[1]
  for w in range(len(axes_refl)):
    for z in [0,1]:
      axes_refl[w][z].set_ylim(0, original_y_lim * 1.2) # resize to not overlap with legend
  fig_refl_RL.savefig(os.path.join("rate_adjustment", f"{species}", output, f"elmm_{species}_peaks.pdf"), bbox_inches="tight")
  plt.close(fig_refl_RL)

  # Guessing the component means and stdevs based on reflected peaks
  final_prominences, init_means, init_stdevs = [], [], []
  for i in range(len(peaks)):
    if significant_peaks_R1[i] == True or significant_peaks_R2[i] == True:
      mean = round(spl_x[peaks[i]], 2)
      init_means.append(mean)
      # The stdev is guess bases on the reflection that generates the tallest prominence for the peak
      index_best_reflect = argmax((reflect_prominences_R1[i], reflect_prominences_R2[i])) # get the "best" reflection (highest prominence)
      prominence = (reflect_prominences_R1[i], reflect_prominences_R2[i])[index_best_reflect] # get prominence of best reflection
      final_prominences.append(prominence)
      width = (reflect_widths_R1[i], reflect_widths_R2[i])[index_best_reflect] # get width of best reflection
      # Assuming that the peak width around half prominence is twice the standard deviation of the underlying Gaussian signal;
      # With the trick of using 100 points per x-axis unit when building the spline, the length between the two spline points 
      # used to get the peak width is in the same scale as the x-axis (i.e. a length of 50 between the two points of the spline
      # is 0.5 Ks and will give a standard deviation of 0.25)
      proxy_stdev = round((width / 2 / 100), 2)
      init_stdevs.append(proxy_stdev)
  return init_means, init_stdevs


def deconvolute_data(tsv_file, max_ks, data_type):
  """
  Generates an artificial proxy dataset for the weighted paranome Ks that doesn't need weights.
  If the dataset is being generated for exp-log mixture model, adds a tail of 1 Ks data beyond 
  the maximum Ks value to avoid fitting EM on a truncated distribution. The data tail comes from
  real data if the maximum Ks is 4 or below, and is made from artificial data otherwise.

  The function obtain the paranome histogram with a very thin bin width (0.01, so 100 bins per Ks unit),
  then it generates the new proxy dataset by taking as values the middle points of each bin and counting them
  according to the histogram bin count. Having narrow bins helps in more reliably representing the actual Ks distribution.
  E.g. the first bin goes from 0 to 0.01 and contains 50 Ks values; the new dataset will contain 50 times the 0.005 value;
  and so on for the remaining bins. 
  The reason for this choice is that the EM algorithm is implemented in a way that is not able to handle the Ks weights.

  :param tsv_file: wgd output file containing either paranome or anchor pairs Ks values (suffix formats: ".ks.tsv", "ks_anchors.tsv")
  :param max_ks: maximum Ks value consideref for the mixture model algorithm
  :param data_type: flag stating whether input Ks are "paralogs" or "anchor pairs"
  :return deconvoluted_data: artificial dataset that produces the same histogram shape, used to avoid having weights in EM
  """
  tail_length = 0.5 # tail spans for 0.5 extra Ks range

  if data_type == "paralogs" or data_type == "anchor pairs":
    ks_data, ks_weights = fc_extract_ks_list.ks_list_from_tsv(tsv_file, max_ks, data_type)
  elif data_type == "orthologs":
    ks_data = fc_extract_ks_list.ks_list_from_tsv(tsv_file, max_ks, data_type)
    ks_weights = [1] * len(ks_data) # dummy weights all equal to 1

  if max_ks <= 4.5:
    # Avoid having a truncated distribution for the EM fitting in the exp-log mixture model
    # If we have real data on the right boundary, let's use an extra Ks as right tail
    # Real data reach only 5 Ks because after that they are not weighted anymore
    max_ks_tail = max_ks + tail_length # the max Ks reached with the tail (0.5 Ks extra than the maximum Ks)
    with open(tsv_file, "r") as tsv_file:
        tsv = read_csv(tsv_file, sep="\t")
        filtered_tsv = tsv.loc[(tsv["Ks"].dtypes == float64) & (tsv["Ks"] >= max_ks) & (tsv["Ks"] <= max_ks_tail)]
        tail_ks = filtered_tsv["Ks"].to_list()
        tail_weights = filtered_tsv["WeightOutliersExcluded"].to_list()
        ks_data.extend(tail_ks)
        ks_weights.extend(tail_weights)
    bin_list_for_deconvoluted_data = arange(0, max_ks_tail + 0.01, 0.01) # must have 0.01 bin width
  else:
    # Ks dataset remains up to max_ks_para without adding tail
    bin_list_for_deconvoluted_data = arange(0, max_ks + 0.01, 0.01) # must have 0.01 bin width
  
  hist_data = histogram(ks_data, bins=bin_list_for_deconvoluted_data, weights=ks_weights)
  deconvoluted_data = array([])
  for i in range(1, len(hist_data[0]) + 1):
      midpoint = round((hist_data[1][i] - 0.01 / 2), 3) # the subtracting factor is half of the bin width (0.01 / 2)
      deconvoluted_data = append(deconvoluted_data, repeat(midpoint, round(hist_data[0][i-1])))
  # NOTE: histogram bin counts are float numbers, while here the artificial bin counts are integers and
  # the total amount of data is not exactly the same (e.g. for A.thaliana is 12124 instead of about 12119)
  
  if max_ks > 4.5:
    # If instead there aren't enough real data to extend for the Ks tail, then add artificial tail
    deconvoluted_data = add_right_tail(deconvoluted_data, hist_data, max_ks, tail_length, tsv_file)
  return deconvoluted_data


def add_right_tail(deconvoluted_data, hist_data, max_ks, tail_length, tsv_file):
  """
  Adds extra right tail data for a range of 1 Ks to avoid the problem of
  fitting the EM on a truncated distribution (truncated at the maximum accepted Ks,
  default 5). The data are added so that they form 100 new histogram bins
  of 0.01 width and of the same height/count as the mean count of the last 
  50 bins of the distribution (corresponding to a range of 0.5 Ks).
  This way, the buffer gaussian can cover the higher Ks without
  being biased by the truncation.

  :param deconvoluted_data: artificial dataset that produces the same histogram shape, used to avoid having weights in EM
  :param hist_data: count data and bins of the real weighted Ks distribution
  :param max_ks: maximum accepted Ks taken from TSV file (default 5)
  :param tail_length: how much does the tail span beyond the maximum Ks (default 0.5 Ks)
  :param tsv_file: wgd output file containing either paranome or anchor pairs Ks values (suffix formats: ".ks.tsv", "ks_anchors.tsv")
  :return deconvoluted_data: same artificial set with the right tail added
  """
  max_ks_tail = max_ks + tail_length # the max Ks reached with the tail (0.5 Ks extra than the maximum Ks)
  for i in arange(max_ks+0.01, max_ks_tail + 0.01, 0.01):
    midpoint = round((i - 0.01 / 2), 3)
    deconvoluted_data = append(deconvoluted_data, repeat(midpoint, round(mean(hist_data[0][-51:]))))
  return deconvoluted_data


def logtransformation(tsv_file, max_ks):
  """
  Performs log-transformation of the paranome Ks values.

  :param tsv_file: wgd output file containing either paranome or anchor pairs Ks values (suffix formats: ".ks.tsv", "ks_anchors.tsv")
  :param max_ks: maximum Ks value to be accepted for the analysis
  :return ks_data_log: log-transformed paranome Ks values
  :return ks_weights_clean: weights associated to the log-transformed paranome Ks values
  """
  ks_data, ks_weights = fc_extract_ks_list.ks_list_from_tsv(tsv_file, max_ks, "paralogs")
  ks_data_clean, ks_weights_clean = remove_ks_zeros(ks_data, ks_weights)
  ks_data_log = log(ks_data_clean)
  return ks_data_log, ks_weights_clean


def remove_ks_zeros(ks_data, ks_weights):
  """
  Removes the zeroes from the original Ks paranome dataset and the weights
  associated to those zeroes, preparing the data for the log-transformation.

  :param ks_data: paranome Ks values of the focal species
  :param ks_weights: weights associated to the paranome Ks values
  :return ks_data_clean: paranome Ks values without zero values
  :return ks_weights_clean: weights associated to non-zero paranome Ks values
  """
  ks_data_clean, ks_weights_clean = array(ks_data), array(ks_weights)
  zero_boolean = ks_data_clean > 0
  ks_data_clean = ks_data_clean[zero_boolean]
  ks_weights_clean = ks_weights_clean[zero_boolean]
  return ks_data_clean, ks_weights_clean


def e_step(num_comp, data, mean, stdev, weights, lambd): 
  """
  Performs the expectation step of the EM algorithm.

  :param num_comp: number of components (1 exponential distribution + N lognormal components)
  :param data: input Ks data
  :param mean: Normal distribution mean vector (one mean per Normal distribution)
  :param stdev: Normal distribution standard deviation vector (one stdev per Normal distribution)
  :param weights: list of component weights
  :param lambd: lambda parameter of the exponential distribution
  :return loglikelihood: measure of how well the model (i.e. the parameters) explains the Ks dataset
  :return posteriors: posteriors or responsabilities
  """
  products = []
  for k in range(0, num_comp):
    if k == 0: # exp component is always number one
      prod = weights[k] * expon.pdf(data, scale=1/lambd)   # expon function wants parameters scale as 1/lambda
    else:
      prod = weights[k] * lognorm.pdf(data, scale=exp(mean[k-1]), s=stdev[k-1])
    products.append(prod)
  sum_comp = sum(products)

  log_sum_comp = log(sum_comp)
  loglikelihood = sum(log_sum_comp)

  posteriors = [] # posteriors
  for k in range(0, num_comp):
    post = products[k] / sum_comp
    posteriors.append(post)

  return loglikelihood, posteriors


def m_step(num_comp, data, posteriors):
  """
  Performs the maximization step of the EM algorithm.  

  :param num_comp: number of components (1 exponential distribution + N lognormal components)
  :param data: input Ks data
  :param posteriors: posterior probabilities (responsibility that component k takes for ‘explaining’ the observation x)
  :return new_means: updated list of Gaussian means after maximization step
  :return new_stdevs: updated list of Gaussian standard deviations after maximization step
  :return new_weights: updated list of component weights after maximization step
  :return new_lambda: updated exponential rate after maximization step
  """
  # Update lambda of the exponential component
  new_lambda = sum(posteriors[0]) / sum(posteriors[0] * data)

  # Computes number of points assigned to each cluster
  points_per_k = []
  for k in range(0, num_comp):
    points_curr_k = sum(posteriors[k])
    points_per_k.append(points_curr_k)

  # Computes updated weights of each component
  new_weights = []
  for k in range(0, num_comp):
    weight = points_per_k[k] / len(data)
    new_weights.append(round(weight, 2))

  # Computes updated means and stdevs of the lognormal components
  new_means, new_stdevs = [], []
  for k in range(0, num_comp-1):
    mean = sum(posteriors[k+1] * log(data)) / points_per_k[k+1]
    new_means.append(mean)

    stdev = sqrt( sum(posteriors[k+1] * pow(log(data) - new_means[k], 2)) / points_per_k[k+1] )
    new_stdevs.append(stdev)
  return new_means, new_stdevs, new_weights, new_lambda


def em(num_comp, max_iter, ks_data_proxy, init_lambd, init_means, init_stdevs, init_weights, model_id,
        max_model_iteration, max_num_comp, parameter_table, outfile, reduced_gaussians_flag=None, model_iteration=None,
        EM_data=False, EM_data_random=False, EM_random=False):
  """
  Performs the EM algorithm with a given number of components, their parameters on the proxy dataset
  for the weighted Ks paranome (the deconvoluted data). If convergence of the algorithm based on
  threshold of 1e-6 is not reached in 300 iterations, it prints a warning.

  :param num_comp: number of components (1 exponential distribution + N lognormal components)
  :param max_iter: maximum number of iterations for the EM steps
  :param ks_data_proxy: artificial dataset (deconvoluted data) used as proxy for the weighted paranome Ks values
  :param init_lambd: initial rate parameter for the exponential component
  :param init_means: list of initial means for the Gaussians related to the lognormal components
  :param init_stdevs: list of initial standard deviations for the Gaussians related to the lognormal components
  :param init_weights: initial component weights (they all have even weights)
  :param model_id: current model identification number
  :param max_model_iteration: number of iterations for a model with a given number of components (e.g. the model with N components is performed 5 times)
  :param max_num_comp: maximum number of components used in the mixture model
  :param parameter_table: list collecting the component parameters
  :param outfile: output file listing the component parameters
  :param reduced_gaussians_flag: boolean flag stating whether the number of detected gaussian was reduced because they were too many
  :param model_iteration: current iteration for the model (the number of iterations is given by max_model_iteration)
  :param EM_data: boolean flag stating whether the current EM algorithm was initializated with "data-driven" strategy
  :param EM_data_random: boolean flag stating whether the current EM algorithm was initializated with "data-driven plus one random component" strategy 
  :param EM_random: boolean flag stating whether the current EM algorithm was initializated with "random components" strategy
  :return bic: BIC score for fitting evaluation
  :return means: list of means of the Gaussians related to the lognormal components after fitting
  :return stdevs: list of standard deviation of the Gaussians related to the lognormal components after fitting
  :return lambd: rate parameter of the exponential component after fitting
  :return weights: list of component weights after fitting
  """
  convergence = False
  ks_data_proxy, __ = remove_ks_zeros(ks_data_proxy, [1]*len(ks_data_proxy))
  for i in range(0, max_iter):
    if i == 0: # initialization
      curr_loglik, posteriors = e_step(num_comp, ks_data_proxy, init_means, init_stdevs, init_weights, init_lambd)
      means, stdevs, weights, lambd = m_step(num_comp, ks_data_proxy, posteriors)

    else: # other iterations
      new_loglik, posteriors = e_step(num_comp, ks_data_proxy, means, stdevs, weights, lambd)
      means, stdevs, weights, lambd = m_step(num_comp, ks_data_proxy, posteriors)

      # Check if the algorithm is converging
      loglik_diff = abs(curr_loglik - new_loglik)
      if loglik_diff < 1e-6:
        convergence = True
        needed_iterations = i
        break
      else:
        curr_loglik = new_loglik
  final_loglik = new_loglik

  bic = compute_bic(final_loglik, num_comp, ks_data_proxy)

  print_details_model(model_id, model_iteration, max_model_iteration, max_num_comp, max_iter, convergence, i, num_comp, reduced_gaussians_flag, init_means, init_stdevs, 
                        init_lambd, init_weights, means, stdevs, lambd, weights, final_loglik, loglik_diff, bic, parameter_table, outfile,
                        EM_data=EM_data, EM_data_random=EM_data_random, EM_random=EM_random)
  return bic, means, stdevs, lambd, weights


def print_details_model(model_id, model_iteration, max_model_iteration, max_num_comp, max_em_iteration, convergence_flag, convergence_iteration, num_comp, reduced_gaussians_flag,
                        init_means, init_stdevs, init_lambd, init_weights, fitted_means, fitted_stdevs, fitted_lambd, fitted_weights, 
                        final_loglik, loglik_diff, bic, parameter_table, outfile, EM_data=False, EM_data_random=False, EM_random=False):
  """
  Prints details about mixture models and adds new line to the parameter
  table which will becaome the TSV output file listing the parameters coming
  from all the computed models, i.e. from each iteration.
  :param model_id: current model identification number
  :param model_iteration: current iteration for the model (the number of iterations is given by max_model_iteration)
  :param max_model_iteration: number of iterations for a model with a given number of components (e.g. the model with N components is performed 5 times)
  :param max_num_comp: maximum number of components used in the mixture model
  :param max_em_iteration: maximum number of iterations allowed for EM algorithm to reach convergence
  :param convergence_flag: boolean flag stating whether the EM algorithm reached convergence or not
  :param convergence_iteration: number of iterations needed by the EM algorithm to reach convergence
  :param num_comp: number of components (1 exponential distribution + N lognormal components)
  :param reduced_gaussians_flag: boolean flag stating whether the number of detected gaussian was reduced because they were too many
  :param init_means: list of initial means for the Gaussians related to the lognormal components
  :param init_stdevs: list of initial standard deviations for the Gaussians related to the lognormal components
  :param init_lambd: initial rate parameter for the exponential component
  :param init_weights: initial component weights (they all have even weights)
  :param fitted_means: list of fitted means for the Gaussians related to the lognormal components
  :param fitted_stdevs: list of fitted standard deviations for the Gaussians related to the lognormal components
  :param fitted_lambd: fitted rate parameter for the exponential component
  :param fitted_weights: fitted component weights
  :param final_loglik: log-likelihood of the fitted model
  :param loglik_diff: log-likelihood difference; if the model didn't convergence, this number is greather than the threshold (1e-6)
  :param bic: BIC score of the fitted model
  :param parameter_table: list collecting the component parameters
  :param outfile: output file listing the component parameters
  :param EM_data: boolean flag stating whether the current EM algorithm was initializated with "data-driven" strategy
  :param EM_data_random: boolean flag stating whether the current EM algorithm was initializated with "data-driven plus one random component" strategy 
  :param EM_random: boolean flag stating whether the current EM algorithm was initializated with "random components" strategy
  """
  # For the more human friendly output format
  outfile.write(f"Model {model_id}\n")
  if not EM_random and not EM_data_random:
    outfile.write("\n")
  if EM_random:
    outfile.write(f"Number of components: {num_comp}\n")
  if EM_data_random or EM_random:
    outfile.write(f"Iteration n. {model_iteration}\n\n")
  outfile.write("")

  if reduced_gaussians_flag:
    outfile.write(f"Reducing the number of lognormal components from {len(init_means)} to 4\n") # todo: 5 ?

  print_parameters("starting", init_means, init_stdevs, init_lambd, init_weights, outfile)

  outfile.write("\n")
  if convergence_flag:
    outfile.write(f"The EM algorithm has reached convergence after {convergence_iteration} iterations\n")
  else:
    outfile.write(f"The EM algorithm didn't reach convergence after {max_em_iteration} iterations (diff is ~{round(loglik_diff, 2)})\n")
  outfile.write("\n")

  print_parameters("fitted", round(fitted_means, 2), round(fitted_stdevs, 2), round(fitted_lambd, 2), round(fitted_weights, 2), outfile)
  outfile.write(f"Log-likelihood: {round(final_loglik, 3)}\n")
  outfile.write(f"BIC: {round(bic, 3)}\n")

  # if EM_data_random or EM_random:
  if (model_id == 2 or model_id == 5) and model_iteration == max_model_iteration:
    outfile.write("\n-----------------------------------------------------\n\n")
  elif model_id == 1 or (model_id != max_num_comp and model_iteration == max_model_iteration):
    outfile.write("\n---------------------\n\n")
  else:
    outfile.write("\n")

  # Add lines for the tabular output format
  if not convergence_flag:
    convergence_iteration = "NA"
  if not model_iteration:
    model_iteration = "NA"
  else:
    model_iteration = int(model_iteration)
  for mean, stdev, weight in zip(fitted_means, fitted_stdevs, fitted_weights[1:]):
    parameter_table.append([model_id, model_iteration, round(bic, 3), round(final_loglik, 3), convergence_iteration,
                      round(fitted_lambd, 2), fitted_weights[0], round(mean, 2), round(stdev, 2), weight])


def make_parameter_table_file(parameter_table, species):
  """
  Generates the text output file with the dataframe containing all component parameters.

  :param parameter_table: list collecting the component parameters
  :param species: informal name of the focal species
  """
  headers = ["Model", "Iteration", "BIC", "Loglikelihood", "Convergence", 
            "Exponential_Rate", "Exponential_Weight", "Normal_Mean", "Normal_SD", "Normal_Weight"]
  parameter_df = DataFrame.from_records(array(parameter_table), columns=headers)
  with open (os.path.join("rate_adjustment", f"{species}", subfolder, f"elmm_{species}_parameters.tsv"), "w+") as outfile:
    outfile.write(parameter_df.to_csv(sep="\t", index=False))


def print_parameters(starting_or_fitted, means, stdevs, lambd, weights, outfile):
  """
  Prints the component (initial or fitted) parameters on screen.

  :param starting_or_fitted: string label to flag if the parameters are the initial or the fitted ones
  :param means: list of means of the Gaussians related to the lognormal components
  :param stdevs: list of standard deviation of the Gaussians related to the lognormal components
  :param lambd: rate parameter of exponential component
  :param weights: list of component weights
  :param outfile: output file listing the component parameters
  """
  if starting_or_fitted == "starting":
    outfile.write("STARTING PARAMETERS:\n")
  elif starting_or_fitted == "fitted":
    outfile.write("FITTED PARAMETERS:\n")
    
  outfile.write(f"  EXP   : {lambd}\n")
  for n in range(len(means)):
      outfile.write(f"  NORM {n+1}: {means[n]} +- {stdevs[n]}\n")
  rounded_weights = []
  for w in weights:
      w = round(w, 2)
      rounded_weights.append(w)
  outfile.write(f"  WEIGHT: {rounded_weights}\n")


def plot_init_comp(ax_ks, ax_logks, means, stdevs, lambd, weights, plot_logtranformed=True):
  """
  Plots the mixture components with initial parameters.
  On the left axis it plots the lognormal and the exponential components,
  on the right axis, if required, the Gaussians related to the lognormals. 

  :param ax_ks: axis object showing Ks paranome in background (left axis)
  :param ax_logks: axis object showing log-transformed Ks paranome in background (right axis)
  :param means: list of initial means of the Gaussians related to the lognormal components
  :param stdevs: list of initial standard deviations of the Gaussians related to the lognormal components
  :param lambd: initial rate parameter of the exponential component
  :param weights: list of component weights after initialization
  :param plot_logtranformed: boolean to flag whether to plot or not also the Gaussians on the right plot with log-transformed data
  """
  x_points = linspace(-5, 10, 15*100)
  ax_ks.plot(x_points, weights[0] * expon.pdf(x_points, scale=1/lambd),'g:', lw=2, alpha=0.5)
  colors = ["b:", "r:", "c:", "m:", "k:"][:len(means)-1] + ["y:"]
  for comp, style in zip(range(0, len(means)), colors):
    ax_ks.plot(x_points, weights[comp+1] * lognorm.pdf(x_points, scale=exp(means[comp]), s=stdevs[comp]), style, lw=1.5, alpha=0.4)
    if plot_logtranformed:
      ax_logks.plot(x_points, weights[comp+1] * norm.pdf(x_points, means[comp], stdevs[comp]), style, lw=1.5, alpha=0.4)


def plot_fitted_comp(ax_ks, ax_logks, means, stdevs, lambd, weights, max_x_axis_lim, peak_stats, correction_table_available, plot_correction_arrows, scaling=1, plot_peak_markers=False, plot_logtranformed=True):
  """
  Plots the mixture components with fitted parameters plus their overall PDF curve.
  On the left axis it plots the lognormal and the exponential components,
  on the right axis, if required, the Gaussians related to the lognormals. 

  :param ax_ks: axis object showing Ks paranome in background  (left axis)
  :param ax_logks: axis object showing log-transformed Ks paranome in background (right axis)
  :param means: list of fitted means of the Gaussians related to the lognormal components
  :param stdevs: list of fitted standard deviations of the Gaussians related to the lognormal components
  :param lambd: fitted rate parameter of the exponential component
  :param weights: list of component weights after fitting
  :param max_x_axis_lim: upper limit in the x-axis 
  :param peak_stats: states whether the cluster peak is intended as median or mode
  :param correction_table_available: boolean stating whether the correction results are available or not
  :param plot_correction_arrows: boolean stating whether there will be plotted correction arrows or not
  :param scaling: proportional factor used to convert density plot into real count plot (in case this is not needed, scaling is set to 1)
  :param plot_peak_markers: boolean flag stating whether to plot markers in correspondence of the peaks of the lognormal components
  :param plot_logtranformed: boolean to flag whether to plot or not also the Gaussians on the right plot with log-transformed data
  """
  x_points = linspace(-5, max_x_axis_lim, int((5 + max_x_axis_lim) *100)) # hundred points per Ks unit
  x_points_strictly_positive = linspace(0, max_x_axis_lim, int(max_x_axis_lim * 100))
  total_pdf_log = 0
  total_pdf = weights[0] * expon.pdf(x_points_strictly_positive, scale=1/lambd)

  if not plot_peak_markers:
    color = "g"
    linestyle = "-"
    alpha = 0.8
  else:  # True only for the final plot, namely for the best model
    color = "dimgrey"  # components in the final plot will be grey, but in the "all_models" plots will be colored
    linestyle = ":"
    alpha = 1

  ax_ks.plot(x_points_strictly_positive, scaling * weights[0] * expon.pdf(x_points_strictly_positive, scale=1/lambd), c=color, ls=linestyle, lw=1.5, alpha=alpha, label='Exponential')
    
  # Getting the representative value of the lognormal peak and sort lognormals accordingly
  lognormal_peaks = {} # key: lognormal component ID, value: peak
  for comp in range(len(means)):
    if peak_stats == "mode":
      peak_coord = round(exp(means[comp] -  pow(stdevs[comp], 2)), 2)
    elif peak_stats == "median":
      peak_coord = round(exp(means[comp]), 2)
    lognormal_peaks[comp] = peak_coord
  lognormals_sorted_by_peak = sorted(lognormal_peaks, key=lognormal_peaks.get) # sort by values (peaks)
  letter_dict = dict(zip(lognormals_sorted_by_peak, [ "a", "b", "c", "d", "e", "f", "g"][:len(lognormals_sorted_by_peak)])) # assign progressive letter

  colors = ["b-", "r-", "c-", "m-", "k-"][:len(means)-1] + ["y-"]

  for comp, style in zip(lognormals_sorted_by_peak, colors):
    color = style[0]
    linestyle = style[1:]
    if plot_peak_markers: # True only for the final plot
      color = "dimgrey"
      linestyle = "--"
      plot_marker(ax_ks, lognormal_peaks[comp], scaling * weights[comp+1] * lognorm.pdf(lognormal_peaks[comp], scale=exp(means[comp]), s=stdevs[comp]), "k", letter_dict[comp], correction_table_available, plot_correction_arrows)
    ax_ks.plot(x_points_strictly_positive, scaling * weights[comp+1] * lognorm.pdf(x_points_strictly_positive, scale=exp(means[comp]), s=stdevs[comp]), c=color, ls=linestyle, lw=1.5, alpha=alpha, label=f'Lognormal {letter_dict[comp]} ({peak_stats} {lognormal_peaks[comp]})')

    if plot_logtranformed:
      ax_logks.plot(x_points, scaling * weights[comp+1] * norm.pdf(x_points, means[comp], stdevs[comp]), c=color, ls=linestyle, lw=1.5, alpha=alpha, label=f'Norm {letter_dict[comp]}')

    total_pdf_log += weights[comp+1] * norm.pdf(x_points, means[comp], stdevs[comp])
    total_pdf += weights[comp+1] * lognorm.pdf(x_points_strictly_positive, scale=exp(means[comp]), s=stdevs[comp])
  total_pdf = total_pdf * scaling

  ax_ks.plot(x_points_strictly_positive, total_pdf, "k-", lw=1.5, label=f'Exp-lognormal mixture model')
  ax_ks.legend(loc="upper right")
  if plot_logtranformed:
    ax_logks.plot(x_points, total_pdf_log, "k-", lw=1.5, label=f'Total PDF')
    ax_logks.legend(loc="upper left")


def plot_histograms_mixture(ax_ks, ax_logks, ks_data, ks_weights, ks_data_log, ks_weights_log, bin_list, bin_width, y_lim, best_model=False):
  """
  Plots the weighted Ks paranome of the focal species on the left axis object,
  where the fitting of exponential and lognormal will be performed.
  Plots the weighted log-transformed Ks paranome on the right axis object, where the
  Gaussians related to the lognormal components will be plotted.

  :param ax_ks: axis object showing Ks paranome in background 
  :param ax_logks: axis object showing log-transformed Ks paranome in background  
  :param ks_data: paranome Ks values of the focal species
  :param ks_weights: weights associated to the paranome Ks values
  :param ks_data_log: log-transformed paranome Ks values
  :param ks_weights_log: weights associated to the log-transformed paranome Ks values
  :param bin_list: list of the edges of each bin (e.g. [0.0, 0.1, 0.2 ... ]); regulates how many bins are there per tick in the x axis 
  :param bin_width: width of the histogram bins
  :param best_model: boolean to flag if the result to be plotted is coming from the best model or not
  """
  fcPlot.plot_histogram("Whole-paranome (weighted)", ax_ks, ks_data, bin_list, None,
                      None, None, weight_list=ks_weights, plot_kde=False, density=True)
  ax_logks.hist(ks_data_log, weights=ks_weights_log, bins=arange(-10, 10 + bin_width, bin_width), histtype='stepfilled', density=True, color="r", alpha=0.3, label=f"Log-transformed paranome")


def compute_bic(loglik, num_comp, ks_data_proxy):
  """
  Computes the BIC score of a fitted model.
  Note that it keeps using the proxy dataset for the weighted paranome Ks,
  previously obtained by deconvoluting the weighted histogram by using bin of
  0.01 width.

  :params loglik: list of log-likelihoods (computed starting from data and by random initialization)
  :param num_comp: number of components (1 exponential distribution + N lognormal components)
  :param ks_data_proxy: artificial dataset (deconvoluted data) used as proxy for the weighted paranome Ks values
  :return bic: BIC score of the fitted model
  """
  bic = -2 * loglik + num_comp * log(len(ks_data_proxy))
  return bic


def eval_best_model(bic_dict, outfile):
  """
  Select the best model according to lowest BIC score.
  Then evaluates the other models by comparing their score to the best one.

  :param bic_dict: dictionary assogining to each model ID its BIC socre
  :param outfile: output file listing the component parameters
  :return best_model_id: ID number of the model with lowest BIC score
  """
  best_bic, best_model_id = float("inf"), 0
  for m in bic_dict:
    if bic_dict[m] < best_bic:
      best_bic = bic_dict[m]
      best_model_id = m

  outfile.write("Model evaluation through BIC score:\n\n")
  outfile.write(f"Number {best_model_id} is the best model (lowest BIC)\n")
  outfile.write(f"BIC list:\n")
  for m in bic_dict:
    outfile.write(f"   Model {m}: {bic_dict[m]}\n")

  l = [
    "0 to  2:   Very weak",
    "2 to  6:    Positive",
    "6 to 10:      Strong",
    "    >10: Very Strong"
  ]
  outfile.write("\n")
  
  if len(bic_dict) != 1:
    outfile.write("Delta BIC assessment: \n")
    for m in bic_dict:
      delta_bic = bic_dict[m] - bic_dict[best_model_id]
      j = 0
      if delta_bic > 2: j = 1
      if delta_bic > 6: j = 2
      if delta_bic > 10: j = 3
      outfile.write(f"   Model {m}: delta(BIC) = {round(delta_bic, 2):>8} ({l[j]})\n")
    outfile.write("\n")
  return best_model_id


def create_legend_mixture_model(axis, legend_size, num_mixture_model_lines):
  """
  Places the legend elements associated to the histograms at the beginning of the legend,\\
  while by default they are placed at the end.

  :param axis: matplotlib axis object from which the legend is taken
  :param legend_size: size of the legend box as a tuple of format "(x, y, width, height)"
  :param num_mixture_model_lines: number of lines generated by mixture models (components + total PDF)
  :return: the updated legend object
  """
  handles, labels = axis.get_legend_handles_labels()
  sorted_handles, sorted_labels = handles.copy(), labels.copy()
  paranome_rect = Patch(facecolor=fcPlot.COLOR_PARANOME_HISTOGRAM, edgecolor="w")
  # empty patch used as spacer between histograms and divergence line legend entries
  empty_rect = Patch(fill=False, edgecolor='none', visible=False)

  sorted_handles = [paranome_rect, empty_rect, sorted_handles[-2]] + sorted_handles[-1-num_mixture_model_lines:-2] + [empty_rect, "Divergence with:"] + sorted_handles[:-1-num_mixture_model_lines]
  sorted_labels = [sorted_labels[-1], "", sorted_labels[-2]] + sorted_labels[-1-num_mixture_model_lines:-2] + ["", ""] + sorted_labels[:-1-num_mixture_model_lines]
  lgd = axis.legend(sorted_handles, sorted_labels, handlelength=1.5, mode="expand", loc="upper left",
                    bbox_to_anchor=legend_size)
  return lgd


def plot_marker(axis, x_value, y_value, color, comp_ID, correction_table_available, plot_correction_arrows, zorder_ID=0):
  """
  Draws a marker on a lognormal component to show where the peak is.
  The peak can be either the mode or the median. 
  Labels the peak with the number associated to the component.

  :param axis: the axis object of the plot
  :param x_value: the peak Ks value
  :param y_value: the peak y-axis value
  :param color: the color associated to the lognormal curve
  :param comp_ID: the letter ID associated to the lognormal component
  :param correction_table_available: boolean stating whether the correction results are available or not
  :param plot_correction_arrows: boolean stating whether there will be plotted correction arrows or not
  :param zorder_ID: takes care of the z-axis (depth) order of the plotted items
  """
  # Plot a vertical line from the median x coordinate up until a fraction of y above the KDE curve
  y_height = axis.get_ylim()[1]
  marker_y_height_fraction = 0.15
  if correction_table_available and plot_correction_arrows:
    ymin = NEGATIVE_Y_FRACTION
  else:
    ymin = 0

  axis.axvline(x=x_value, ymin=ymin, ymax=y_value / y_height + marker_y_height_fraction, color=color,
              linestyle=(0, (3, 3)), alpha=0.4, linewidth=1.1, solid_capstyle='butt', solid_joinstyle='miter',
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
  # white component with letter associated to the lognormal component
  axis.text(x=x_value, y=y_value + (y_height * marker_y_height_fraction), s=comp_ID, fontsize=10,
            horizontalalignment='center', verticalalignment='center', clip_on=True, zorder=zorder_ID + 5, color='w')


def plot_best_model(fig_best_model, ax_best_model, species, ks_data, ks_weights, bin_list, bin_width, max_x_axis_lim,
                    y_lim, best_model_id, all_models_init_parameters, all_models_fitted_parameters,
                    correction_table, correction_table_available, consensus_peak_for_multiple_outgroups,
                    peak_stats, color_list, plot_correction_arrows, deconvoluted_data, max_ks_EM):
  """
  Plots the backgound paranome histogram distribution with the divergence lines,
  then plots the initial components as dotted curves and the fitted components as solid curves. 

  :param fig_best_model: figure object that will show the best mixture model
  :param ax_best_model: axis object where the best model result will be plotted
  :param species: informal name of the focal species
  :param ks_data: paranome Ks values of the focal species
  :param ks_weights: weights associated to the paranome Ks values
  :param bin_list: list of the edges of each bin (e.g. [0.0, 0.1, 0.2 ... ]); regulates how many bins are there per tick in the x axis 
  :param bin_width: width of the histogram bins
  :param max_x_axis_lim: upper limit in the x-axis 
  :param best_model_id: ID number of the model with lowest BIC score
  :param all_models_init_parameters: dictionary assigning to each model ID its initial parameters
  :param all_models_fitted_parameters: dictionary assigning to each model ID its fitted parameters
  :param correction_table: correction results in DataFrame format (contains both possible types of consensus strategy for how to deal with multiple outgroups)
  :param correction_table_available: boolean stating whether the correction results are available or not
  :param consensus_peak_for_multiple_outgroups: choice of how to deal with multiple corrections
          for the same divergence (expected values: either "mean among outgroups" or "best outgroup")
  :param peak_stats: states whether the cluster peak is intended as median or mode
  :param color_list: list of colors that will be used to distinguish all divergence lines
  :param plot_correction_arrows: boolean stating whether there will be plotted correction arrows or not
  :param deconvoluted_data: artificial dataset that produces the same histogram shape, used to avoid having weights in EM
  :param max_ks_EM: upper limit of the Ks range considered for the mixture modeling fitting
  """
  hist_paranome = fcPlot.plot_histogram("Whole-paranome (weighted)", ax_best_model, ks_data, bin_list, None, None, None, ks_weights, plot_kde=False, density=False)

  # Set the height of the distribution according to heighest histogram bin
  if y_lim is None:
    fcPlot.set_mixed_plot_height(ax_best_model, y_lim, hist_paranome)

  # PLOTTING THE ORTHOLOG DIVERGENCE LINES on the paralog distribution
  if correction_table_available:
    dummy_fig, dummy_axis = plt.subplots()
    fcPlot.plot_divergences(correction_table, peak_stats, consensus_peak_for_multiple_outgroups, dummy_axis, ax_best_model, color_list, plot_correction_arrows)

  init_means, init_stdevs, init_lambd, init_weights = all_models_init_parameters[best_model_id]

  final_means, final_stdevs, final_lambd, final_weights = all_models_fitted_parameters[best_model_id]

  scaling = bin_width * len(deconvoluted_data[deconvoluted_data <= max_ks_EM])
  plot_fitted_comp(ax_best_model, None, final_means, final_stdevs, final_lambd, final_weights, max_x_axis_lim, peak_stats, correction_table_available, plot_correction_arrows, scaling=scaling, plot_peak_markers=True, plot_logtranformed=False)

  legend_size = fcPlot.define_legend_size(ax_best_model)
  chart_box = ax_best_model.get_position()
  ax_best_model.set_position([chart_box.x0, chart_box.y0, chart_box.width*0.65, chart_box.height])
  lgd = create_legend_mixture_model(ax_best_model, legend_size, len(init_means)+2) # number of plotted lines is: exp + lognormals + total PDF
  fig_best_model.savefig(os.path.join("rate_adjustment", f"{species}", f"mixed_{species}_elmm.pdf"),
                    bbox_extra_artists=(lgd, fig_best_model._suptitle), bbox_inches="tight", transparent=True, format="pdf")

  # # TEMPORARY FOR A FIGURE PLOT WITH DENSITY FOR COMPARISON AFTER SCALING
  # fig, ax = generate_best_model_figure("elaeis", "Elaeis guineensis", 3, None, True, True)
  # fcPlot.plot_histogram("Whole-paranome (weighted)", ax, ks_data, bin_list, None, None, None, ks_weights, plot_kde=False, density=True)
  # if correction_table_available:
  #   fcPlot.plot_divergences(correction_table, consensus_peak_for_multiple_outgroups, dummy_axis, ax, color_list, plot_correction_arrows)
  # plot_fitted_comp(ax, None, final_means, final_stdevs, final_lambd, final_weights, peak_stats, scaling=1, plot_peak_markers=True, plot_logtranformed=False)

  # legend_size = fcPlot.define_legend_size(ax)
  # chart_box = ax.get_position()
  # ax.set_position([chart_box.x0, chart_box.y0, chart_box.width*0.65, chart_box.height])
  # lgd = create_legend_mixture_model(ax, legend_size, len(init_means)+2) # number of plotted lines is: exp + lognormals + total PDF
  # fig.savefig(os.path.join("rate_adjustment", f"{species}", f"mixed_density_elmm.pdf"),
  #                   bbox_extra_artists=(lgd, fig._suptitle), bbox_inches="tight", transparent=True, format="pdf")


