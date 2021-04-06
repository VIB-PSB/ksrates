import sys
import os
import pandas
from ast import literal_eval
import logging
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import ksrates.fc_check_input as fcCheck
import ksrates.fc_kde_bootstrap as fcPeak
import ksrates.fc_plotting as fcPlot
import ksrates.fc_configfile as fcConf
from ksrates.utils import init_logging
import matplotlib
# Use the Agg (Anti-Grain Geometry) backend to avoid needing a graphical user interface (X11 backend)
matplotlib.use('Agg')


def plot_orthologs_distr(config_file, trios_file):
    # INPUT
    config = fcConf.Configuration(config_file)
    init_logging("Plotting ortholog distributions for all ortholog trios", config.get_logging_level())
    logging.info("Loading parameters and input files")

    # Get parameters from configuration file
    species_of_interest = config.get_species()
    latin_names = config.get_latin_names()
    max_ks_ortho = config.get_max_ks_ortho()
    bin_width_ortho = config.get_bin_width_ortho()
    bin_list_ortho = fcPlot.get_bins(max_ks_ortho, bin_width_ortho)
    x_lim = config.get_x_lim_ortho()

    # Get input file listing the trios
    default_path_trios_file = os.path.join("rate_adjustment", f"{species_of_interest}", f"ortholog_trios_{species_of_interest}.tsv")
    trios_file = fcCheck.get_argument_path(trios_file, default_path_trios_file, "Trios TSV file")
    if trios_file == "":
        logging.error(f"Trios TSV file not found at default position [{default_path_trios_file}]")
        logging.error("Exiting")
        sys.exit(1)
    with open(trios_file, 'r') as f1:
        trios = pandas.read_csv(f1, sep="\t")

    # Get the ortholog Ks list database (to plot histograms; mandatory input file)
    ks_list_db_path = config.get_ks_db()
    fcCheck.check_inputfile(ks_list_db_path, "Ortholog Ks list database")
    with open(ks_list_db_path, 'r') as f2:
        ks_list_db = pandas.read_csv(f2, sep="\t", index_col=0)

    # Get the ortholog peak database (to plot distribution mode and median; not mandatory)
    db_path = config.get_ortho_db()
    no_peak_db = False
    try:
        with open(db_path, 'r') as f3:
            db = pandas.read_csv(f3, sep="\t", index_col=0)
    except Exception:
        no_peak_db = True
        logging.warning(f"Ortholog Ks peak database empty or not found at the path provided in the config file: distribution peaks will not be shown")

    # -----------------------------------------------------------------------------

    # GENERATING PDF FIGURE with ortholog distributions FOR EACH TRIO
    outgroups_per_divergent_pair_dict = {}
    missing_pairs_ks_list, missing_pairs_peaks = [], []

    for __, row in trios.iterrows():
        species, sister, out = row['Species'], row['Sister_Species'], row['Out_Species']
        # Generate dictionary of divergent pairs linked with their outgroups
        divergent_pair_key = f"{species}_{sister}"
        if divergent_pair_key not in outgroups_per_divergent_pair_dict.keys():
            outgroups_per_divergent_pair_dict[divergent_pair_key] = [out]
        else:
            outgroups_per_divergent_pair_dict[divergent_pair_key].append(out)


    # PLOTTING THE DISTRIBUTIONS
    for divergent_pair in outgroups_per_divergent_pair_dict.keys():
        with PdfPages(os.path.join("rate_adjustment", f"{species_of_interest}", f"orthologs_{divergent_pair}.pdf")) as pdf:
            species = divergent_pair.split("_")[0]
            latinSpecies = latin_names[species]
            sister = divergent_pair.split("_")[1]
            latinSister = latin_names[sister]

            logging.info("")
            logging.info(f"Plotting ortholog Ks distributions for species pair [{latinSpecies} - {latinSister}]")

            # tags, e.g. A.filiculoides_S.cucullata
            species_sister = "_".join(sorted([latinSpecies, latinSister], key=str.casefold))

            if species_sister in list(ks_list_db.index):
                ks_list_species_sister = literal_eval(ks_list_db.at[species_sister, 'Ks_Values'])
                # Getting 20 KDE curves through bootstrap
                logging.info(f"  Plotting data for the two sister species [{latinSpecies} - {latinSister}]")
                bootstrap_kde_species_sister = fcPeak.bootstrap_KDE(ks_list_species_sister, 20, x_lim,
                                                                            bin_width_ortho)
            else:  # If pair not in Ks list database
                if species_sister not in missing_pairs_ks_list:
                    missing_pairs_ks_list.append(species_sister)

            out_list = outgroups_per_divergent_pair_dict[divergent_pair]
            for out in out_list:
                latinOut = latin_names[out]
                logging.info(f"- Using outspecies [{latinOut}]:")

                fig, axes = fcPlot.generate_orthologs_figure(latinSpecies, latinSister, latinOut, x_lim)

                # tags, e.g. A.filiculoides_S.cucullata
                species_out = "_".join(sorted([latinSpecies, latinOut], key=str.casefold))
                sister_out = "_".join(sorted([latinSister, latinOut], key=str.casefold))

                # Plotting Ks lists and their KDE lines
                if species_sister in list(ks_list_db.index):
                    fcPlot.plot_orthologs_histogram_kdes(ks_list_species_sister, bin_list_ortho, axes[0],
                                                                bootstrap_kde_species_sister)

                if species_out in list(ks_list_db.index):
                    ks_list = literal_eval(ks_list_db.at[species_out, 'Ks_Values'])
                    # Getting 20 KDE curves through bootstrap
                    logging.info(f"  Plotting data for focal species and outspecies [{latinSpecies} - {latinOut}]")
                    bootstrap_kde = fcPeak.bootstrap_KDE(ks_list, 20, x_lim, bin_width_ortho)
                    # Plotting Ks lists and their KDE lines
                    fcPlot.plot_orthologs_histogram_kdes(ks_list, bin_list_ortho, axes[1], bootstrap_kde)
                else:  # If pair not in Ks list database
                    if species_out not in missing_pairs_ks_list:
                        missing_pairs_ks_list.append(species_out)                    

                if sister_out in list(ks_list_db.index):
                    ks_list = literal_eval(ks_list_db.at[sister_out, 'Ks_Values'])
                    # Getting 20 KDE curves through bootstrap
                    logging.info(f"  Plotting data for sister species and outspecies [{latinSister} - {latinOut}]")
                    bootstrap_kde = fcPeak.bootstrap_KDE(ks_list, 20, x_lim, bin_width_ortho)
                    # Plotting Ks lists and their KDE lines
                    fcPlot.plot_orthologs_histogram_kdes(ks_list, bin_list_ortho, axes[2], bootstrap_kde)
                else:  # If pair not in Ks list database
                    if sister_out not in missing_pairs_ks_list:
                        missing_pairs_ks_list.append(sister_out)   

                # Plotting estimated mode of the orthologs distributions as vertical lines
                y_upper_lim = axes[0].get_ylim()[1] 
                if not no_peak_db: # if the peak database is available
                    
                    if species_sister in list(db.index):
                        fcPlot.plot_orthologs_peak_lines(db, species_sister, axes[0], y_upper_lim)
                    else:
                        if species_sister not in missing_pairs_peaks:
                            missing_pairs_peaks.append(species_sister)

                    if species_out in list(db.index):
                        fcPlot.plot_orthologs_peak_lines(db, species_out, axes[1], y_upper_lim)
                    else:
                        if species_out not in missing_pairs_peaks:
                            missing_pairs_peaks.append(species_out)

                    if sister_out in list(db.index):
                        fcPlot.plot_orthologs_peak_lines(db, sister_out, axes[2], y_upper_lim)
                    else:
                        if sister_out not in missing_pairs_peaks:
                            missing_pairs_peaks.append(sister_out)

                object_list = []
                for ax in axes:
                    lgd = ax.get_legend()
                    object_list.append(lgd)
                sup = fig._suptitle
                object_list.append(sup)

                pdf.savefig(fig, transparent=True, bbox_extra_artists=(sup,), bbox_inches='tight')
                plt.close()
        logging.info(f"- Saving PDF figure [orthologs_{divergent_pair}.pdf]")

    # Report if species are missing from any of the two ortholog databases
    if len(missing_pairs_ks_list) != 0 or len(missing_pairs_peaks) != 0:
        logging.warning("")
        logging.warning("One or more species pairs are missing from the ortholog databases,")
        logging.warning("their plot will be incomplete.")
        logging.warning("")

        missing_in_both_dbs = list((set(missing_pairs_peaks) & set(missing_pairs_ks_list)))
        if len(missing_in_both_dbs) != 0:
            logging.warning("Species missing from both Ks peak and Ks list ortholog databases:")
            for pair in sorted(missing_in_both_dbs):
                logging.warning(f"- {pair}")
            logging.warning("")

        missing_pairs_peaks = list(set(missing_pairs_peaks) - set(missing_in_both_dbs))
        if len(missing_pairs_peaks) != 0:
            logging.warning("Species missing from the ortholog Ks peak database:")
            for pair in sorted(missing_pairs_peaks):
                logging.warning(f"- {pair}")
            logging.warning("")

        missing_pairs_ks_list = list(set(missing_pairs_ks_list) - set(missing_in_both_dbs))
        if len(missing_pairs_ks_list) != 0:
            logging.warning("Species missing from the ortholog Ks list database:")
            for pair in sorted(missing_pairs_ks_list):
                logging.warning(f"- {pair}")
            logging.warning("")

        logging.warning("Please compute their ortholog Ks data and/or add the ortholog data to the databases,")
        logging.warning("then rerun this step.")

    logging.info("")
    logging.info("All done")