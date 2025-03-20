from math import sqrt
import logging

# Filenames
_ADJUSTMENT_TABLE_ALL = "adjustment_table_{}_all.tsv"
_ADJUSTMENT_TABLE = "adjustment_table_{}.tsv"


def decompose_ortholog_ks(ortholog_db, idx_species_sister, idx_species_outgroup, idx_sister_outgroup, peak_stats):
    """
    Starting from a trio of species that includes the focal species, a sister species and their outgroup,
    the function computes the branch-specific Ks contributions of the focal species and of the sister species
    by applying principles from the relative rate test. It also computes the standard deviation associated to the
    branch-specific Ks contributions through error propagation rules.
    The branch-specific Ks contribution of the focal species accounts for the synonymous substitutions accumulated in the focal species
        since the divergence with the sister species, it is therefore a relative measure of the absolute substitution rate of the focal species.
    The branch-specific Ks contribution of the sister species accounts for the synonymous substitutions accumulated in the sister species 
        since the divergence with the focal species, it is therefore a relative measure of the absolute substitution rate of the sister species.

    :param ortholog_db: database of ortholog Ks distribution peaks
    :param idx_species_sister: string composed by species and sister latin names separated by underscore;
        used to search in the database (format example: "Azolla filiculoides_Salvinia cucullata")
    :param idx_species_outgroup: string composed by species and outspecies latin names separated by underscore;
        used to search in the database (same format)
    :param idx_sister_outgroup: string composed by sister and outspecies latin names separated by underscore;
        used to search in the database (same format)
    :param peak_stats: the statistics used to get the ortholog distribution peak (either "mode" or "median")
    :return: rel_rate_species, the branch-specific Ks contribution of the focal species
    :return: rel_rate_species_sd, the standard deviation associated to rel_rate_species
    :return: rel_rate_sister, the branch-specific Ks contribution of the sister species
    :return: rel_rate_sister_sd, the standard deviation associated to rel_rate_sister
    """
    peak_sp_sis = ortholog_db.at[idx_species_sister, 'Mode']
    sd_sp_sis = ortholog_db.at[idx_species_sister, 'Mode_SD']

    peak_sp_out = ortholog_db.at[idx_species_outgroup, 'Mode']
    sd_sp_out = ortholog_db.at[idx_species_outgroup, 'Mode_SD']

    peak_sis_out = ortholog_db.at[idx_sister_outgroup, 'Mode']
    sd_sis_out = ortholog_db.at[idx_sister_outgroup, 'Mode_SD']

    # RRT formulas
    rel_rate_species = (peak_sp_out + peak_sp_sis - peak_sis_out) / 2.0 # also called k_AO
    rel_rate_sister = (peak_sp_sis + peak_sis_out - peak_sp_out) / 2.0 # also called k_BO
    # Error propagation rules
    rel_rate_species_sd = sqrt(pow(sd_sp_out, 2) + pow(sd_sp_sis, 2) + pow(sd_sis_out, 2)) / 2.0 # also called k_AO_sd
    rel_rate_sister_sd = sqrt(pow(sd_sp_sis, 2) + pow(sd_sis_out, 2) + pow(sd_sp_out, 2)) / 2.0 # also called k_BO_sd

    focal_species = idx_species_sister.split("_")[0]
    sister_species = idx_species_sister.split("_")[1]
    
    # Checkpoint to spot negative Ks distances for the focal and/or sister species.
    # Prints a warning, doesn't stop the pipeline
    if rel_rate_species < 0:
        logging.warning("")
        logging.warning(f"ksrates has returned a negative number ({round(rel_rate_species, 5)}) for the genetic distance accumulated by")
        logging.warning(f"focal species [{focal_species}] since its divergence with sister species {sister_species}")
        logging.warning(f"(refer to the 'K_OA' segment in Supplementary materials).")
    if rel_rate_sister < 0:
        logging.warning("")
        logging.warning(f"ksrates has returned a negative number ({round(rel_rate_sister, 5)}) for the genetic distance accumulated by")
        logging.warning(f"sister species {sister_species} since its divergence with focal species [{focal_species}]")
        logging.warning(f"(refer to the 'K_OB' segment in Supplementary materials).")
    if rel_rate_species < 0 or rel_rate_sister < 0:
        logging.warning("")
        logging.warning("Negative distances are an artefact due to poor estimate of the ortholog distribution peaks")
        logging.warning(f"used as input values for the rate-adjustment formulas.")
        logging.warning("")
        logging.warning("Poor peak estimation is often due to the very old divergence age between")
        logging.warning(f"the two chosen species and/or outgroup, or by bad genome/sequence quality.")
        logging.warning("It is thus recommended to visually check the quality of the ortholog distributions and their peak estimates")
        logging.warning(f"in the 'orthologs_species1_species2.pdf' output files.")
        logging.warning("")
        logging.warning("TIPS! Try rerunning your analysis with one or more of the following changes:")
        logging.warning(" - limit the number of distant outgroups by decreasing 'max_number_outgroups' (default is 4)")
        logging.warning(" - set 'consensus_mode_for_multiple_outgroups' to 'best outgroup' (see Documentation)")
        logging.warning(" - remove species responsible for very old divergences in the input phylogeny")
        logging.warning("")

    return rel_rate_species, rel_rate_species_sd, rel_rate_sister, rel_rate_sister_sd


def compute_corrected_ks_species_sister(rel_rate_species, rel_rate_species_sd):
    """
    Performs substitution rate-adjustment on the original peak of the ortholog distribution between
    the focal species (species A) and sister species (species B).

    :param rel_rate_species: branch-specific Ks contribution of the focal species, namely the synonymous substitutions occurred 
        in the focal species since the speciation event with the sister species; also called K_OA
    :param rel_rate_species_sd: standard deviation associated to rel_rate_species; also called K_OA_sd
    :return: adjusted_peak, rate-adjusted ortholog Ks distribution peak of focal species and sister species
    :return: adjusted_peak_sd, standard deviation associated to adjusted_peak
    """
    corrected_peak = rel_rate_species * 2
    corrected_peak_sd = sqrt(2) * rel_rate_species_sd

    return corrected_peak, corrected_peak_sd


def interpretation_adjusted_plot(focal_species_latin, consensus_peak_for_multiple_outgroups, 
                                wgd_peaks, correction_table):
    """
    Prints an automatic interpretation of the rate-adjusted mixed Ks plots in which for each
    inferred (putative) WGD peak it is said which species share it with the focal species.
    This is done by comparing the x-coordinate of WGD peaks coming from mixture
    modeling analyses and the rate-adjusted ortholog modes.
    The WGD x-coordinate (its Ks age) is either coming from a component when the method is 
    exponential-lognormal mixture modeling and lognormal mixture modeling, or from an 
    anchor cluster when the method is anchor clustering.
    
    :param focal_species_latin: latin name of focal species
    :param consensus_peak_for_multiple_outgroups: choice of how to deal with multiple adjustments
       for the same divergence (expected values: either "mean among outgroups" or "best outgroup")
    :param wgd_peaks: dictionary associating to a wgd letter the coordinate of its peak (either component or anchor cluster)
    :param correction_table: adjustment results in DataFrame format (contains both possible types of consensus strategy for how to deal with multiple outgroups)
    """
    if consensus_peak_for_multiple_outgroups == "mean among outgroups":
        # Then consider the columns in the correction_table that are generated using the average method
        column_header_peak = "Adjusted_Mode_Mean"
    elif consensus_peak_for_multiple_outgroups == "best outgroup":
        # Then consider the columns in the correction_table that are generated using the best outgroup method
        column_header_peak = "Adjusted_Mode_Best"

    species_interpretation = {}
    shared_wgd_interpretation = {}

    # Generate dictionary associating each divergent species to its mode
    ortholog_modes = {}
    for diverging_species, mode in zip(correction_table["Sister_Species"], correction_table[column_header_peak]):
        ortholog_modes[diverging_species] = mode
        # Also initialize a dictionary for the final interpretation logging output
        species_interpretation[diverging_species] = []

    logging.info("")
    logging.info(f"Automatic interpretation of rate-adjusted mixed Ks plot:")
    logging.info("Please revise it manually due to known overfitting of mixture modeling methods!")
    # Figuring out which species share which (putative) WGDs
    if len(wgd_peaks) == 0:
        logging.info(" - There are no inferred (putative) WGD peaks.")
        return

    for peak in wgd_peaks:  # For each wgd peak
        shared_wgd_interpretation[peak] = []
        for diverging_species in ortholog_modes:  # For each divergent species
            # What about when it's "=="?
            if ortholog_modes[diverging_species] < wgd_peaks[peak]:
                species_interpretation[diverging_species].append(peak)
                shared_wgd_interpretation[peak].append(diverging_species)
    
    # For every peak, say which species share it with the focal species
    for peak_letter in wgd_peaks:
        # If wgd peak is not shared with any other species
        if len(shared_wgd_interpretation[peak_letter]) == 0:
            logging.info(f' - Inferred (putative) WGD signature "{peak_letter}" is ' +
                         f'specific to focal species {focal_species_latin}')
            continue
        # If wgd peak is shared with all of the other species in the plots
        elif len(shared_wgd_interpretation[peak_letter]) == len(ortholog_modes):
            logging.info(f' - Inferred (putative) WGD signature "{peak_letter}" is shared ' + 
                        f"with all involved species in the provided phylogeny")
            continue
        # If wgd peak is shared with only one other species
        elif len(shared_wgd_interpretation[peak_letter]) == 1:
            logging.info(f' - Inferred (putative) WGD signature "{peak_letter}" is shared with {shared_wgd_interpretation[peak_letter][0]}')
        else: # If wgd peak is shared with some other species
            logging.info(f' - Inferred (putative) WGD signature "{peak_letter}" is shared ' + 
                         f"with the following species:")
            for species in shared_wgd_interpretation[peak_letter]:
                logging.info(f"   {species}")
    logging.info("")
    return