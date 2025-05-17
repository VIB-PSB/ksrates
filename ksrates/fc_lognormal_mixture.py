import os
import logging
from pandas import DataFrame
import matplotlib.pyplot as plt
import seaborn as sns
import ksrates.fc_plotting as fcPlot
from ksrates.fc_exp_log_mixture import deconvolute_data
from numpy import exp, sqrt, linspace, arange, log, argmin, mean, around, array
from wgd_ksrates.modeling import mixture
import scipy.stats as ss
from wgd_ksrates.modeling import plot_aic_bic
from wgd_ksrates.modeling import plot_mixture, plot_probs
import ksrates.fc_exp_log_mixture as fcEM
from matplotlib.patches import Patch
from ksrates.fc_cluster_anchors import subfolder

def plot_all_models_gmm(models, data, l, u, bins, out_file):
    """
    Modified from wgd.
    Plot a bunch of GMMs. Increased fig height to fully host
    a legend with 6 or 7 components.

    :param models: list of GMM model objects
    :param data: Ks array
    :param l: lower Ks limit
    :param u: upper Ks limit
    :param bins: number of histogram bins
    :param out_file: output file
    :return: nada
    """
    fig, axes = plt.subplots(len(models), 3, figsize=(15, 3.5 * len(models)))
    for i, model in enumerate(models):
        plot_mixture(model, data, axes[i, 0], l, u, bins=bins)
        plot_mixture(model, data, axes[i, 1], log=True, l=log(l + 0.0001),
                     u=log(u), bins=bins)
        plot_probs(model, axes[i, 2], l, u)
    sns.despine(offset=5)
    fig.tight_layout()
    fig.savefig(out_file)


def inspect_aic(aic, outfile):
    """
    Evaluate the Akaike information criterion results for mixture models.
    Modified from wgd.

    :param aic: AIC values
    :param outfile: output file where to add the AIC scores
    """
    im = argmin(aic)
    outfile.write("\n")
    outfile.write("AIC assessment:\n")
    outfile.write("min(AIC) = {:.2f} for model {}\n".format(aic[im], im + 1))
    outfile.write("Relative probabilities compared to model {}:\n".format(im + 1))
    outfile.write("   /                          \\\n")
    outfile.write("   |      (min(AIC) - AICi)/2 |\n")
    outfile.write("   | p = e                    |\n")
    outfile.write("   \                          /\n")
    for i, aic_i in enumerate(aic):
        p_i = exp((aic[im] - aic_i) / 2)
        outfile.write(".. model{:>4}: p = {:.4f}\n".format(i + 1, p_i))
    outfile.write("\n")


def inspect_bic(bic, outfile):
    """
    Evaluate the BIC values
    Modified from wgd.

    :param bic: BIC values
    :param outfile: output file where to add the BIC scores
    """
    im = argmin(bic)
    l = [
        "0 to  2:   Very weak",
        "2 to  6:    Positive",
        "6 to 10:      Strong",
        "    >10: Very Strong"
    ]
    outfile.write("\n")
    outfile.write("Delta BIC assessment: \n")
    outfile.write("min(BIC) = {:.2f} for model {}\n".format(bic[im], im + 1))
    for i, bic_i in enumerate(bic):
        dbic = bic_i - bic[im]
        j = 0
        if dbic > 2: j = 1
        if dbic > 6: j = 2
        if dbic > 10: j = 3
        outfile.write(".. model{:>4}: delta(BIC) = {:>8.2f} ({})\n".format(
                i + 1, dbic, l[j]))
    outfile.write("\n")


def log_components(X, model_id, m, outfile, parameter_table, max_iter):
    """
    Modified from wgd.

    :param X: data frame (log transformed Ks values)
    :param model_id: ID number of the model 
    :param m: model
    :param outfile: output file listing the component parameters
    :param parameter_table: list collecting the component parameters
    """
    outfile.write(f"Model {model_id}\n")
    outfile.write(f"Number of components: {len(m.means_)}\n")
    if m.n_iter_ == max_iter:
        convergence = "NA"
        outfile.write(f"The EM algorithm didn't reach convergence after {max_iter} iterations\n")
    else:
        convergence = m.n_iter_
        outfile.write(f"The EM algorithm has reached convergence after {m.n_iter_} iterations\n\n")

    outfile.write("FITTED PARAMETERS:\n")
    for j in range(len(m.means_)):
        outfile.write(f"  NORM {j+1}: {m.means_[j][0]} +- {sqrt(m.covariances_[j][0][0])}\n")
        
        parameter_table.append([model_id, m.bic(X), m.lower_bound_, convergence, 
                                m.means_[j][0], m.covariances_[j][0][0], m.weights_[j]])

    weight_list = []
    for w in m.weights_:
        weight_list.append(w)
    outfile.write(f"  WEIGHT: {weight_list}\n")
    outfile.write(f"Log-likelihood: {m.lower_bound_}\n")
    outfile.write(f"BIC: {m.bic(X)}\n\n")


def fit_gmm(X, n1, n2, outfile, parameter_table, max_iter=600, n_init=1, **kwargs):
    """
    Modified from wgd.
    Compute Gaussian mixtures for different numbers of components

    :param X: data frame (log transformed Ks values)
    :param n1: minimum number of components
    :param n2: maximum number of components
    :param outfile: output file listing the component parameters
    :param parameter_table: list collecting the component parameters
    :param max_iter: maximum number of iterations
    :param n_init: number of k-means initializations
    :param kwargs: other keyword args for `GaussianMixture`
    :return: models, bic, aic, best model
    """
    # fit models with 1 to n components
    N = arange(n1, n2 + 1)
    models = [None for i in range(len(N))]
    for i in range(len(N)):
        # Apply LMM only if the number of components is smaller than or equal to the number of data points
        if N[i] <= len(X):
            models[i] = mixture.GaussianMixture(
                    n_components=N[i], covariance_type='full', max_iter=max_iter,
                    n_init=n_init, tol=1e-6, **kwargs
            ).fit(X)
            log_components(X, i+1, models[i], outfile, parameter_table, max_iter)
        else:
            logging.warning(f"Lognormal mixture model with {N[i]} or more components is skipped due to too few input Ks data points")
            break
    # Remove any residual None from the list (None's are left in the list if there are
    # too many components for the number of data points, e.g. only 2 anchor points)
    models = [m for m in models if m]

    # compute the AIC and the BIC
    aic = [m.aic(X) for m in models]
    bic = [m.bic(X) for m in models]
    best = models[argmin(bic)]

    return models, bic, aic, best


def plot_mixture_model(model, data, max_x_axis_lim, ax, bin_width, scaling, peak_stats, correction_table_available, plot_correction_arrows, l=0, u=5, color='black', alpha=0.2,
                 log=False, bins=25):
    """
    Plot a mixture model. Assumes a log-transformed model and data
    and will back-transform.
    Modified from wgd.

    Note from scipy docs:
    ---------------------
    A common parametrization for a lognormal random variable Y
    is in terms of the mean, mu, and standard deviation, sigma,
    of the unique normally distributed random variable X such
    that exp(X) = Y. This parametrization corresponds to setting
    s = sigma and scale = exp(mu).

    So since we fit a normal mixture on the logscale, we can get
    the lognormal by setting scale to np.exp(mean[k]) and s to
    np.sqrt(variance[k]) for component k.

    :param model: best model
    :param data: data array
    :param max_x_axis_lim: upper limit in the x-axis 
    :param ax: figure ax
    :param bin_width: width of the histogram bins
    :param scaling: proportional factor used to convert density plot into real count plot (in case this is not needed, scaling is set to 1)
    :param peak_stats: states whether the cluster peak is intended as median or mode
    :param correction_table_available: boolean stating whether the adjustment results are available or not
    :param plot_correction_arrows: boolean stating whether there will be plotted adjustment arrows or not
    :param l: lower Ks limit
    :param u: upper Ks limit
    :param color: color for histogram
    :param alpha: alpha value
    :return: ax
    """
    x = linspace(l, max_x_axis_lim, 1000).reshape((-1, 1))
    if not log:
        data = exp(data)

    means = model.means_
    varcs = model.covariances_
    weights = model.weights_
    total_pdf = None
    first = True

    lognormal_peaks = {} # key: lognormal component ID, value: peak

    for comp in range(len(means)):
        # Getting the lognormal peak and sort lognormals accordingly
        if peak_stats == "mode":
            peak_coord = around(exp(means[comp][0] -  pow(sqrt(varcs[comp][0][0]), 2)), 2)
        elif peak_stats == "median":
            peak_coord = around(exp(means[comp]), 2)
        lognormal_peaks[comp] = peak_coord
    lognormals_sorted_by_peak = sorted(lognormal_peaks, key=lognormal_peaks.get) # sort by values (peaks)
    letter_dict = dict(zip(lognormals_sorted_by_peak, [ "a", "b", "c", "d", "e", "f", "g"][:len(lognormals_sorted_by_peak)])) # assign progressive letter

    for comp in lognormals_sorted_by_peak:
        fcEM.plot_marker(ax, lognormal_peaks[comp], scaling * weights[comp] * ss.lognorm.pdf(lognormal_peaks[comp], scale=exp(means[comp][0]), s=sqrt(varcs[comp][0][0])), "k", letter_dict[comp], correction_table_available, plot_correction_arrows)

    for comp in lognormals_sorted_by_peak:
        curve = scaling * weights[comp] * ss.lognorm.pdf(x, scale=exp(means[comp]), s=sqrt(varcs[comp]))
        ax.plot(x, curve, c="dimgrey", ls='--', lw=1.5, label=f"Lognormal {letter_dict[comp]} ({peak_stats} {lognormal_peaks[comp]})")

        if first:
            total_pdf = curve
            first = False
        else:
            total_pdf += curve
    ax.plot(x, total_pdf, '-k', lw=1.5, label="Lognormal mixture model")
    return ax


def lmm(
    fig, max_x_axis_lim, data_type, tsv_file, species, axis, ks_range, min_ks_anchors,
    components, bins, bin_width_para, max_iter, n_init,
    output_dir, outfile, parameter_table, datatype_tag, peak_stats, correction_table_available, plot_correction_arrows):
    """
    Modified from wgd code.

    Mixture modeling tools.
    Note that histogram weighting is done after applying specified filters. Also
    note that mixture models are fitted to node-averaged (not weighted)
    histograms. Please interpret mixture model results with caution.
    :param fig: figure object
    :param max_x_axis_lim: upper limit in the x-axis 
    :param data_type: strings stating whether the data are "paralogs" or "anchor pairs"
    :param tsv_file: wgd output file containing either paranome or anchor pairs Ks values (suffix formats: ".ks.tsv", "ks_anchors.tsv")
    :param species: informal name of the focal species
    :param axis: axis object
    :param ks_range: Ks range used for models
    :param min_ks_anchors: minimum anchor Ks value to be modelled 
    :param components: number of components to use (tuple: (min, max))
    :param bins: number histogram bins for visualization
    :param bin_width_para: bin width of paralog Ks histogram
    :param max_iter: number of iterations
    :param n_init: number of k-means initializations (best is kept)
    :param output_dir: output folder
    :param outfile: output file for logging details
    :param parameter_table: list collecting the component parameters
    :param datatype_tag: string for figure title stating if data is paranome or anchors
    :param peak_stats: states whether the cluster peak is intended as median or mode
    :param correction_table_available: boolean stating whether the adjustment results are available or not
    :param plot_correction_arrows: boolean stating whether there will be plotted adjustment arrows or not
    """
    # Generating artificial dataset with same shape as the Ks histogram with 0.01 bin width
    deconvoluted_data = deconvolute_data(tsv_file, ks_range[1], data_type, min_ks_anchors)
    deconvoluted_data = log(deconvoluted_data).reshape(-1, 1) # log-transform and reshape data for fit_gmm

    models, bic, aic, best = fit_gmm(
            deconvoluted_data, components[0], components[1], outfile, parameter_table,
            max_iter=max_iter, n_init=n_init)

    inspect_aic(aic, outfile)
    inspect_bic(bic, outfile)
    # For now, let's not generate the AIC and BIC PDF plots
    # logging.info("Plotting AIC and BIC")
    # plot_aic_bic(aic, bic, components[0], components[1],
    #                  os.path.join(output_dir, f"lmm_{species}_aic_bic_{datatype_tag}.pdf"))

    logging.info("Plotting mixtures")
    plot_all_models_gmm(models, deconvoluted_data, ks_range[0], ks_range[1], bins=bins,
                        out_file=os.path.join(output_dir, f"lmm_{species}_all_models_{datatype_tag}.pdf"))
    
    # Plotting the components of the best model on the final picture;
    # Components are scaled up to the size of actual count data and not to density histogram
    scaling = bin_width_para * len(deconvoluted_data)
    plot_mixture_model(best, deconvoluted_data, max_x_axis_lim, axis, bin_width_para, scaling, peak_stats,
                       correction_table_available, plot_correction_arrows, ks_range[0], ks_range[1], bins=bins)
    return best


def create_legend_mixture_model(axis, legend_size, num_mixture_model_lines, datatype):
    """
    Places the legend elements associated to the histograms at the beginning of the legend,\\
    while by default they are placed at the end.

    :param axis: matplotlib axis object from which the legend is taken
    :param legend_size: size of the legend box as a tuple of format "(x, y, width, height)"
    :param num_mixture_model_lines: number of lines generated by mixture models (components + total PDF)
    :param datatype: string for figure title stating if data is "paranome" or comes from "colinearity" analysis
    :return: the updated legend object
    """
    handles, labels = axis.get_legend_handles_labels()
    sorted_handles, sorted_labels = handles.copy(), labels.copy()
    if datatype == "paranome":
        paralog_rect = Patch(facecolor=fcPlot.COLOR_PARANOME_HISTOGRAM, edgecolor="w")
    elif datatype == "anchors":
        paralog_rect = Patch(facecolor=fcPlot.COLOR_ANCHOR_HISTOGRAM, edgecolor="w")
    # empty patch used as spacer between histograms and divergence line legend entries
    empty_rect = Patch(fill=False, edgecolor='none', visible=False)

    sorted_handles = [paralog_rect, empty_rect, sorted_handles[num_mixture_model_lines-1]] + sorted_handles[:num_mixture_model_lines-1] + [empty_rect, "Divergence with:"] + sorted_handles[num_mixture_model_lines:-1]
    sorted_labels = [sorted_labels[-1], "", sorted_labels[num_mixture_model_lines-1]] + sorted_labels[:num_mixture_model_lines-1] + ["", ""] + sorted_labels[num_mixture_model_lines:-1]

    lgd = axis.legend(sorted_handles, sorted_labels, handlelength=1.5, mode="expand", loc="upper left",
                        bbox_to_anchor=legend_size)
    return lgd


def save_lmm(fig, axis, species, best_model, datatype, correction_table_available):
    """
    This function must be called to save the mixed plot figure with lognormal mixture model
    in order to adjust the figure layout: the plot area is shrunk to the left and some reasonable 
    space is left on the right side for the legend.

    :param fig: figure object of the adjusted mixed distribution
    :param axis: axis object of the adjusted mixed distribution
    :param species: focal species
    :param datatype: string for figure title stating whether data is paranome or comes from "colinearity" analysis
    """
    if correction_table_available:
        num_mixture_model_lines = len(best_model.means_) + 1 # components + total PDF

        legend_size = fcPlot.define_legend_size(axis)
        chart_box = axis.get_position()

        axis.set_position([chart_box.x0, chart_box.y0, chart_box.width*0.65, chart_box.height])
        lgd = create_legend_mixture_model(axis, legend_size, num_mixture_model_lines, datatype)

        fig.savefig(os.path.join("rate_adjustment", f"{species}", f"mixed_{species}_lmm_{datatype}.pdf"),
                    bbox_extra_artists=(axis, lgd, fig._suptitle), bbox_inches="tight", transparent=True, format="pdf")
    else:
        # if not correction_table_available use a simpler layout with the legend
        # inside the plot and no right margin
        lgd = axis.legend(handlelength=1.5, loc="upper right")
        fig.savefig(os.path.join("rate_adjustment", f"{species}", f"mixed_{species}_lmm_{datatype}.pdf"),
                    transparent=True, format="pdf")


def make_parameter_table_file(parameter_table, species, datatype):
  """
  Generates the text output file with the dataframe containing all component parameters.

  :param parameter_table: list collecting the component parameters
  :param species: informal name of the focal species
  :param datatype: string for figure title stating whether data is paranome or comes from "colinearity" analysis
  """
  #headers = ["Model", "Iteration", "BIC", "Loglikelihood", "Convergence", 
  #          "Exponential_Rate", "Exponential_Weight", "Normal_Mean", "Normal_SD", "Normal_Weight"]
  headers = ["Model", "BIC", "Loglikelihood", "Convergence", "Normal_Mean", "Normal_SD", "Normal_Weight"]
  parameter_df = DataFrame.from_records(array(parameter_table), columns=headers)
  with open (os.path.join("rate_adjustment", f"{species}", subfolder, f"lmm_{species}_parameters_{datatype}.tsv"), "w+") as outfile:
    outfile.write(parameter_df.to_csv(sep="\t", index=False))
