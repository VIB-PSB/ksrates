import os
from pandas import read_csv
import matplotlib.pyplot as plt
from ksrates.utils import init_logging
import logging
import sys
from numpy import argmin, arange, random
import ksrates.fc_configfile as fcConf
import ksrates.fc_check_input as fcCheck
import ksrates.fc_plotting as fcPlot
import ksrates.fc_extract_ks_list as fc_extract_ks_list
from ksrates.fc_plotting import COLOR_ANCHOR_HISTOGRAM
import ksrates.fc_exp_log_mixture as fcEM
from ksrates.fc_cluster_anchors import subfolder
from ksrates.fc_rrt_correction import _ADJUSTMENT_TABLE


def exp_log_mixture(config_file, expert_config_file, paralog_tsv_file, correction_table_file):
  # INPUT
  config = fcConf.Configuration(config_file, expert_config_file)
  init_logging("Exponential-Lognormal mixture model on Ks paranome", config.get_logging_level())
  logging.info("Loading parameters and input files")

  # GET PARAMETERS and INPUT FILES
  species = config.get_species()
  latin_names = config.get_latin_names()
  latinSpecies = latin_names[species]
  species_escape_whitespace = latinSpecies.replace(' ', '\ ')
  # NOTE: only for this script, the max accepted Ks value is constrained to 5,
  # instead of being taken from configfile, due to the fact that the code uses
  # a built-in "buffer" lognormal to cover the area of high Ks values (4-5 Ks). 
  max_ks_para = config.get_max_ks_para()
  bin_width = config.get_bin_width_para()
  bin_list = fcPlot.get_bins(max_ks_para, bin_width)
  x_max_lim = config.get_x_max_lim()
  y_max_lim = config.get_y_lim()
  color_list = config.get_color_list()
  plot_correction_arrows = config.plot_correction_arrows()
  paranome_analysis = config.get_paranome()
  # Getting the statistical measure for how to determine the representative value of an ortholog distribution
  peak_stats = config.get_peak_stats() # default is mode (other option, median)
  # Getting the choice on how to deal with the presence of multiple adjustments for the same divergent pair
  consensus_peak_for_multiple_outgroups = config.get_consensus_peak_for_multiple_outgroups()

  # Parameters used during the mixture modeling
  max_ks_EM = config.get_max_ks_for_mixture_model(max_ks_para) # upper Ks limit considered for the mixture model fitting
  max_EM_iterations = config.get_max_EM_iterations() # default 300
  num_EM_initializations = config.get_num_EM_initializations() # how many times the fitting with N given components is initialized 
                                                               # in the random method and in the "peak + random" method). default 10
  max_num_comp = config.get_max_mixture_model_components() # max number of components used in the fitting with random components (exp + buffer lognormal + lognormal)
  min_num_comp = 2 # there are always at least two components (the exponential and the buffer lognormal)

  logging.info(f" - maximum EM iterations: {max_EM_iterations}")
  logging.info(f" - number of EM initializations: {num_EM_initializations}")
  logging.info(f" - maximum number of components: {max_num_comp}")
  if max_ks_EM != max_ks_para:
    logging.info(f" - Ks range considered for the mixture modeling: up to {max_ks_EM} Ks.")
  logging.info("")

  if paranome_analysis:
    default_path_paralog_tsv_file = os.path.join("paralog_distributions", f"wgd_{species}", f"{species}.ks.tsv")
    paralog_tsv_file = fcCheck.get_argument_path(paralog_tsv_file, default_path_paralog_tsv_file, "Paralog Ks TSV file")
    if paralog_tsv_file == "":
      logging.error(f"Paralog Ks TSV file not found at default position [{default_path_paralog_tsv_file}].")
      logging.error("Exiting")
      sys.exit(1)
  else:
    logging.error("Mixture modeling is not performed since paranome analysis is not required in configuration file")
    logging.error("Exiting")
    sys.exit(0) # exit code 0 because no actual errors were thrown

  ks_data, ks_weights = fc_extract_ks_list.ks_list_from_tsv(paralog_tsv_file, max_ks_para, "paralogs")

  # Get adjustment results TSV file
  # If correction_table is (still) missing, it will be equal to empty string (""), but the script will not exit
  default_path_correction_table_file = os.path.join("rate_adjustment", f"{species}", f"{_ADJUSTMENT_TABLE.format(species)}")
  correction_table_file = fcCheck.get_argument_path(correction_table_file, default_path_correction_table_file, "Rate-adjustment table file")
  if correction_table_file == "":
      logging.warning("Rate-adjustment data are not available yet, only Ks paranome distribution will be plotted.")
      correction_table = None
      correction_table_available = False
  else:
      with open(correction_table_file, "r") as f:
          correction_table = read_csv(f, sep="\t")
          correction_table_available = True

  # Creating folder for secondary output files
  output = os.path.join(subfolder)
  if not os.path.isdir(os.path.join("rate_adjustment", species, output)):
      logging.info(f"Creating directory [rate_adjustment/{species}/{output}]")
      logging.info("")
      os.makedirs(os.path.join("rate_adjustment", species, output))

  # Generating figures for the mixture models
  fig_peaks, ax_peaks_ks, ax_peaks_logks, ax_peaks2_ks, ax_peaks2_logks, sup_peaks = fcEM.generate_peak_model_figure(
                                                                                species_escape_whitespace, x_max_lim)
  fig_random, axes_random, sup_random = fcEM.generate_random_model_figure(species_escape_whitespace, 
                                                                          min_num_comp, max_num_comp, x_max_lim)
  fig_best_model, ax_best_ks = fcEM.generate_best_model_figure(latinSpecies, x_max_lim, y_max_lim, 
                                                               correction_table_available, plot_correction_arrows)

  # Generating a proxy dataset for the weighted Ks paranome (deconvoluting the histogram)
  deconvoluted_data = fcEM.deconvolute_data(paralog_tsv_file, max_ks_EM, "paralogs")
  # Log-transformation of Ks paranome
  ks_data_log, ks_weights_log = fcEM.logtransformation(paralog_tsv_file, max_ks_EM)

  bic_dict, parameters_list = {}, {} # will contain BIC scores and parameters of all models 
  all_models_init_parameters = {} # will contain the initial parameters of all models (for plotting purpose)
  all_models_fitted_parameters = {} # will contain the fitted parameters of all models (for plotting purpose)
  parameter_table = [] # will contain parameters for every model iteration for tabular output text file

  # -----------------------------------------------------------------------------

  with open (os.path.join("rate_adjustment", f"{species}", subfolder, f"elmm_{species}_parameters.txt"), "w+") as outfile:
    
    # Performing EM algorithm multiple times with different types of initializations

    logging.info("Performing EM algorithm with initialization from Ks paranome data")
    model_id = 1
    ax_peaks_ks.set_title(f"Model {model_id}\n")
    # Plotting background Ks paranome histograms for the figures (original and log-transformed)
    fcEM.plot_histograms_mixture(ax_peaks_ks, ax_peaks_logks, ks_data, ks_weights, ks_data_log, ks_weights_log, bin_list, bin_width, y_max_lim)

    # Initializing component parameters and plotting them
    init_lambd, init_means, init_stdevs, init_weights, reduced_gaussians = fcEM.init_parameters_from_data(ax_peaks_ks, species, ks_data_log,
                                                                                                          ks_weights_log, ks_data, ks_weights, 
                                                                                                          species_escape_whitespace, output, max_ks_EM)
    all_models_init_parameters[model_id] = [init_means, init_stdevs, init_lambd, init_weights]
    fcEM.plot_init_comp(ax_peaks_ks, ax_peaks_logks, init_means, init_stdevs, init_lambd, init_weights)
    num_comp = len(init_means) + 1

    # Performing EM algorithm, computing the BIC and plotting the fitted components
    bic_peaks, means_peaks, stdevs_peaks, lambd_peaks, weights_peaks  = fcEM.em(num_comp, max_EM_iterations, deconvoluted_data, init_lambd, init_means, 
                                                                                init_stdevs, init_weights, model_id, 
                                                                                num_EM_initializations, max_num_comp, parameter_table, outfile,
                                                                                reduced_gaussians_flag=reduced_gaussians, EM_data=True)
    all_models_fitted_parameters[model_id] = [means_peaks, stdevs_peaks, lambd_peaks, weights_peaks]
    bic_dict[model_id] = bic_peaks
    fcEM.plot_fitted_comp(ax_peaks_ks, ax_peaks_logks, means_peaks, stdevs_peaks, lambd_peaks, weights_peaks, x_max_lim, peak_stats, correction_table_available, plot_correction_arrows)

    # -----------------------------------------------------------------

    logging.info("Performing EM algorithm with initialization from Ks paranome data plus a random lognormal component")
    model_id = 2
    ax_peaks2_ks.set_title(f"Model {model_id}")
    # Plotting background Ks paranome histograms for the figures (original and log-transformed)
    fcEM.plot_histograms_mixture(ax_peaks2_ks, ax_peaks2_logks, ks_data, ks_weights, ks_data_log, ks_weights_log, bin_list, bin_width, y_max_lim)

    # Adding a random lognormal to the components initialized from Ks data
    # It is done multiple times and only the best result is retained
    bic_from_same_num_comp = []
    start_parameters, final_parameters = [], []
    for i in range(num_EM_initializations):
      if len(init_means) > 4:
        # Limiting the number of total lognormal to 5:
        # 1 buffer lognormal, plus a maximum of 3 "peak" lognormals from data, plus the 1 random lognormal that is going to be added;
        # if there are already more than 4 lognormal, some of the peak lognormal are removed, while the buffer lognormal is always left
        updated_means = list(list(random.choice(init_means[:-1], size=4, replace=False)) + [init_means[-1]])
        updated_stdevs = list(list(random.choice(init_stdevs[:-1], size=4, replace=False)) + [init_stdevs[-1]])
        reduced_gaussians = True
      else:
        updated_means, updated_stdevs = init_means.copy(), init_stdevs.copy()
        reduced_gaussians = False

      updated_means.append(round(random.choice(arange(-0.5, 1, 0.1)), 1))
      updated_stdevs.append(round(random.choice(arange(0.3, 0.9, 0.1)), 1))
      num_comp = len(updated_means) + 1
      updated_weights = [1/num_comp] * num_comp
      start_parameters.append([updated_means, updated_stdevs, init_lambd, updated_weights])
      
      # Performing EM algorithm and computing the BIC
      bic_peaks, means_peaks, stdevs_peaks, lambd_peaks, weights_peaks  = fcEM.em(num_comp, max_EM_iterations, deconvoluted_data, init_lambd,  
                                                                                  updated_means, updated_stdevs, updated_weights, model_id, 
                                                                                  num_EM_initializations, max_num_comp,
                                                                                  parameter_table, outfile, reduced_gaussians_flag=reduced_gaussians,
                                                                                  model_iteration=i+1, EM_data_random=True)
      bic_from_same_num_comp.append(bic_peaks)
      final_parameters.append([means_peaks, stdevs_peaks, lambd_peaks, weights_peaks])

    # Get the initial and fitted parameters of the best model and plot them; get its BIC score
    updated_means, updated_stdevs, init_lambd, updated_weights = start_parameters[argmin(bic_from_same_num_comp)]
    all_models_init_parameters[model_id] = [updated_means, updated_stdevs, init_lambd, updated_weights]
    fcEM.plot_init_comp(ax_peaks2_ks, ax_peaks2_logks, updated_means, updated_stdevs, init_lambd, updated_weights)

    final_means, final_stdevs, final_lambd, final_weights = final_parameters[argmin(bic_from_same_num_comp)]
    all_models_fitted_parameters[model_id] = [final_means, final_stdevs, final_lambd, final_weights]
    fcEM.plot_fitted_comp(ax_peaks2_ks, ax_peaks2_logks, final_means, final_stdevs, final_lambd, final_weights, x_max_lim, peak_stats, correction_table_available, plot_correction_arrows)

    bic_dict[model_id] = min(bic_from_same_num_comp)
    parameters_list[model_id] = final_parameters[argmin(bic_from_same_num_comp)]

    plt.close()
    fig_peaks.savefig(os.path.join("rate_adjustment", f"{species}", output, f"elmm_{species}_models_data_driven.pdf"),
                      bbox_inches="tight", bbox_extra_artists=(sup_peaks,), format="pdf")
    logging.info(f"Saving PDF figure of mixture models [{species}/{output}/elmm_{species}_models_data_driven.pdf]")
    logging.info("")

  # -----------------------------------------------------------------------------
      
    logging.info("Performing EM algorithm with (almost) random initialization")

    num_comp_list = arange(min_num_comp, max_num_comp + 1) # range(3, 6) is 2, 3, 4, 5
    axes_ids = num_comp_list - min_num_comp # 0, 1, 2, 3
    model_ids = num_comp_list - (min_num_comp - 3) # 3, 4, 5, 6 # model 1 and 2 are from data; from model 3 on is from random method
    for num_comp, ax_id, model_id in zip(num_comp_list, axes_ids, model_ids):
      logging.info(f" - using {num_comp} components")
      ax_rand_ks, ax_rand_logks = axes_random[ax_id][0], axes_random[ax_id][1]
      # Plotting background Ks paranome histograms for the figures (original and log-transformed)
      fcEM.plot_histograms_mixture(ax_rand_ks, ax_rand_logks, ks_data, ks_weights, ks_data_log, ks_weights_log, bin_list, bin_width, y_max_lim)

      bic_from_same_num_comp = []
      start_parameters, final_parameters = [], []
      for i in range(num_EM_initializations):
        # Initializing parameters for the given number of components
        init_means, init_stdevs, init_lambd, init_weights = fcEM.init_parameters_randomly(model_id, num_comp, ax_rand_ks, max_ks_EM)
        start_parameters.append([init_means, init_stdevs, init_lambd, init_weights])

        # Performing EM algorithm and computing the BIC
        bic_random, means_random, stdevs_random, lambd_random, weights_random = fcEM.em(num_comp, max_EM_iterations, deconvoluted_data,
                                                                                        init_lambd, init_means, init_stdevs, init_weights,
                                                                                        model_id, num_EM_initializations,
                                                                                        max_num_comp, parameter_table, outfile, model_iteration=i+1, EM_random=True)
        bic_from_same_num_comp.append(bic_random)
        final_parameters.append([means_random, stdevs_random, lambd_random, weights_random])

      # Get the initial and fitted parameters of the best result for the given number of components and plot them; get its BIC score
      init_means, init_stdevs, init_lambd, init_weights = start_parameters[argmin(bic_from_same_num_comp)]
      all_models_init_parameters[model_id] = [init_means, init_stdevs, init_lambd, init_weights]
      fcEM.plot_init_comp(ax_rand_ks, ax_rand_logks, init_means, init_stdevs, init_lambd, init_weights)

      final_means, final_stdevs, final_lambd, final_weights = final_parameters[argmin(bic_from_same_num_comp)]
      all_models_fitted_parameters[model_id] = [final_means, final_stdevs, final_lambd, final_weights]
      fcEM.plot_fitted_comp(ax_rand_ks, ax_rand_logks, final_means, final_stdevs, final_lambd, final_weights, x_max_lim, peak_stats, correction_table_available, plot_correction_arrows)

      bic_dict[model_id] = min(bic_from_same_num_comp)
      parameters_list[model_id] = final_parameters[argmin(bic_from_same_num_comp)]
        
    plt.close()
    fig_random.savefig(os.path.join("rate_adjustment", f"{species}", output, f"elmm_{species}_models_random.pdf"), bbox_inches="tight",
                                    bbox_extra_artists=(sup_random,), format="pdf")
    logging.info(f"Saving PDF figure of mixture models [{species}/{output}/elmm_{species}_models_random.pdf]")
    logging.info("")

    # Generating tabular text file with all model parameters 
    fcEM.make_parameter_table_file(parameter_table, species)

    # BIC evaluation of all models (from data, from data with random lognormal and from random components)
    logging.info("Models are evaluated according to their BIC score.")
    # Get best model by lowest BIC score and plot it; print comparison with the other models
    best_model_id = fcEM.eval_best_model(bic_dict, outfile)
    fcEM.plot_best_model(fig_best_model, ax_best_ks, species, ks_data, ks_weights, bin_list, bin_width, x_max_lim,
                        y_max_lim, best_model_id, all_models_init_parameters, all_models_fitted_parameters, 
                        correction_table, correction_table_available, consensus_peak_for_multiple_outgroups,
                        peak_stats, color_list, plot_correction_arrows, deconvoluted_data, max_ks_EM)
  logging.info(f"Saving PDF figure of best mixture model [mixed_{species}_elmm.pdf]")
  logging.info("")
  logging.info("All done")