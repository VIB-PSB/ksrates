import os
import pandas
from pandas import DataFrame
import random
import logging
import numpy as np
from scipy import stats
import ksrates.fc_check_input as fcCheck
import ksrates.fc_extract_ks_list as fc_extract_ks_list


def compute_kde(ks_list, max_ks, bin_width, bandwidth_scaling=None):
    """
    Returns the KDE object and the KDE data point from the x axis and from the y axis.
    The KDE is rescaled because the default Scott's factor is always making the curve too smooth:
    - If bootstrap_random is given as bandwidth_scaling, the bandwidth factor is multiplied by a random number
    between 0.5 and 0.9, so that the resulting KDE becomes tighter to the underlying distribution with some
    variability.
    - If bootstrap_average is given as bandwidth_scaling, the bandwidth factor is multiplied by the average
    fraction from the previous option, which is 0.7, to have a single line representing 
    the average modification of the KDE.

    :param ks_list: list of ortholog Ks values to be plotted
    :param max_ks: maximum x value allowed in orthologs distribution plots
    :param bin_width: width of histogram bins in ortholog Ks plots (default value in config file: 0.1;\\
    there are 10 bins per tick in the x axis)
    :param bandwidth_scaling: if set, specifies the method used to scale the default (Scott's)
        bandwidth, either "bootstrap_random" or "bootstrap_average"
    :return: kde, gaussian_kde Scipy object
    :return: kde_x, array of x data points of KDE line
    :return: kde_y, array of y data points of KDE line
    """
    kde = stats.gaussian_kde(ks_list, bw_method='scott')
    if bandwidth_scaling:
        if bandwidth_scaling == "bootstrap_random":
            variable_bw = kde.factor * np.random.uniform(0.5, 0.9)
            kde.set_bandwidth(variable_bw)
        elif bandwidth_scaling == "bootstrap_average":
            average_bw = kde.factor * 0.7
            kde.set_bandwidth(average_bw)
        else:
            logging.warning(f"fc_kde_bootstrap.compute_kde(): unsupported value [\"{bandwidth_scaling}\"]\n\
        for parameter bandwidth_scaling, KDE bandwidth will not be scaled.")

    kde_x = np.linspace(-0.1, max_ks + 0.1, num=512)
    kde_y = kde(kde_x)
    kde_y = kde_y * len(ks_list) * bin_width  # rescaling KDE to match histogram distribution
    return kde, kde_x, kde_y


def get_mode_from_kde(kde_x, kde_y):
    """
    Takes a KDE object from SciPy and returns the x coordinate associated
    to the highest y value (i.e. the KDE mode).

    :param kde_x: x axis data point of KDE object
    :param kde_y: y axis data point of KDE object
    :return: the x coordinate associated to the KDE mode
    :return: the y coordinate associated to the KDE mode
    """
    index_max_y = np.argmax(kde_y)
    mode_x_value = kde_x[index_max_y]
    return mode_x_value, max(kde_y)


def bootstrap_peak(ks_list, n_replicates, max_ks, bin_width):
    """
    Provides a statistical estimate of the peak (indended as mode or median) of an ortholog Ks distribution.
    It obtains N (n_replicates) replicates with replacements from the ortholog Ks list.
    Then it computes the KDE and finds its mode and median.
    Returns the mean mode and the mean median across the replicates, together with their standard deviations.
    
    :param ks_list: list of the original ortholog Ks value used for boostrap
    :param n_replicates: number of bootstrap iterations
    :param max_ks: maximum x value allowed in ortholog Ks distribution plots (default: 10)
    :param bin_width: width of histogram bins in ortholog Ks plots (default value in config file: 0.1;
        there are 10 bins per tick in the x axis)
    :return: mean_mode: mean mode obtained from all bootstrap replicates
    :return: std_mode: standard deviation associated to mean_mode
    :return: mean_median: mean median obtained from all bootstrap replicates
    :return: std_median: standard deviation associated to mean_median
    """
    mode_list = []
    median_list = []
    for r in range(n_replicates):
        sample_list = random.choices(ks_list, k=len(ks_list))
        __, kde_x, kde_y = compute_kde(sample_list, max_ks, bin_width, bandwidth_scaling="bootstrap_random")

        mode, __ = get_mode_from_kde(kde_x, kde_y)
        mode_list.append(mode)

        median_list.append(np.median(sample_list))

    mean_mode = np.mean(mode_list, dtype=np.float64)
    std_mode = np.std(mode_list, dtype=np.float64)

    mean_median = np.mean(median_list, dtype=np.float64)
    std_median = np.std(median_list, dtype=np.float64)
    return mean_mode, std_mode, mean_median, std_median


def bootstrap_KDE(ks_list, n_replicates, max_ks, bin_width):
    """
    Returns a list containing data points of some KDE lines computed via bootstrap.
    They are plotted to have an idea of the bootstrap iterations.

    :param ks_list: list of ortholog Ks values for bootstrap
    :param n_replicates: number of iterations for boostrap (default: 20) 
    :param max_ks: maximum x value allowed for ortholog Ks distribution plots
    :param bin_width: width of histogram bins in ortholog Ks plots (default value in config file: 0.1;
        there are 10 bins per tick in the x axis)
    :return boostrap_kde: list of x and y data points of 20 KDEs 
    """
    bootstrap_kde = []
    for r in range(n_replicates):
        sample_list = random.choices(ks_list, k=len(ks_list))
        __, kde_x, kde_y = compute_kde(sample_list, max_ks, bin_width, bandwidth_scaling="bootstrap_random")
        bootstrap_kde.append([kde_x, kde_y])
    return bootstrap_kde



def estimate_peak(species1, species2, latinSp1, latinSp2, max_ks_ortho, n_iter, x_lim_ortho, bin_width_ortho, ks_list_db_path, db_path, flag_not_in_peak_db=True, flag_not_in_ks_db=True):
    """
    Extracts the Ks list of a species pair, estimates its peak and updates databases. 
    The species pair was previously checked and either is missing the peak or the ks_list or both in the databases.

    :param species1: first species of the pair (alphabetically ordered, case insensitive)
    :param species2: second species of the pair (alphabetically ordered, case insensitive)
    :param latinSp1: latin name of the first species of the pair
    :param latinSp2: latin name of the second species of the pair
    :param max_ks_ortho: maximum ortholog Ks value accepted when extracting Ks values from the TSV file (default: 10)
    :param n_iter: number of bootstrap iterations when estimating the ortholog distribution peak
    :param x_lim_ortho: upper limit of the x axis range for the ortholog distribution plot (default: 5)
    :param bin_width_ortho: width of the ortholog Ks histogram bins (default: 0.1, 10 per unit)
    :param ks_list_db_path: filename/path to the database of ortholog Ks lists
    :param db_path: filename/path to the database of ortholog peaks
    :param flag_not_in_peak_db: flag for the presence/absence of the species pair in the ortholog peak database (True/False)
    :param flag_not_in_ks_db: flag for the presence/absence of the species pair in the ortholog Ks list database (True/False)
    :return: a flag to state if the current peak computation failed due to missing ortholog Ks TSV file (True/False)
    """
    compute_peak_failed = False

    logging.info(f"{latinSp1} and {latinSp2}:")
    logging.info("- Extracting ortholog Ks list")
    default_path_tsv_file = os.path.join("ortholog_distributions", f"wgd_{species1}_{species2}", f"{species1}_{species2}.ks.tsv")
    tsv_path = fcCheck.check_file_existence_and_content_in_default_paths(default_path_tsv_file,  "Ortholog Ks TSV file")
    if tsv_path != "": # if the TSV is present
        ks_list = fc_extract_ks_list.ks_list_from_tsv(tsv_path, max_ks_ortho, "orthologs")

        if flag_not_in_peak_db:    # if the pair is missing in the ortholog peak database
            logging.info(f"- Computing distribution peak through bootstrap ({n_iter} iterations)")
            mean_peak, std_peak, mean_median, std_median = bootstrap_peak(ks_list, n_iter, x_lim_ortho,
                                                                                bin_width_ortho)
            db_new_row = DataFrame([[latinSp1, latinSp2, mean_peak, std_peak]],
                                    index=[f"{latinSp1}_{latinSp2}"])
            with open(db_path, "a+") as outfile_db:
                logging.info("- Adding peak to ortholog peak database")
                outfile_db.write(db_new_row.to_csv(sep="\t", header=None))

        if flag_not_in_ks_db:      # if the pair is missing in the ks list database
            logging.info("- Adding Ks list to ortholog Ks list database")
            ks_list_new_row = DataFrame([[latinSp1, latinSp2, ks_list]], index=[f"{latinSp1}_{latinSp2}"])
            with open(ks_list_db_path, "a+") as outfile_ksList:
                outfile_ksList.write(ks_list_new_row.to_csv(sep="\t", header=None))
    else:
        logging.warning(f"  Ortholog Ks TSV file not found [{species1}_{species2}.ks.tsv]. Skipping peak estimate.") 
        compute_peak_failed = True

    return compute_peak_failed