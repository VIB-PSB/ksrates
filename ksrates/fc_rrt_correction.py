from math import sqrt

# Filenames
_ADJUSTMENT_TABLE_ALL = "adjustment_table_{}_all.tsv"
_ADJUSTMENT_TABLE = "adjustment_table_{}.tsv"


def compute_ks_distances(ortholog_db, idx_species_sister, idx_species_outgroup, idx_sister_outgroup, peak_stats):
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
    if peak_stats == "mode": # choosing the MODE as peak Ks for the WGD event
        peak, sd = 'Ortholog_Mode', 'Ortholog_Mode_SD'
    elif peak_stats == "median": # choosing the MEDIAN as peak Ks for the WGD event
        peak, sd = 'Ortholog_Median', 'Ortholog_Median_SD'        

    peak_sp_sis = ortholog_db.at[idx_species_sister, peak]
    sd_sp_sis = ortholog_db.at[idx_species_sister, sd]

    peak_sp_out = ortholog_db.at[idx_species_outgroup, peak]
    sd_sp_out = ortholog_db.at[idx_species_outgroup, sd]

    peak_sis_out = ortholog_db.at[idx_sister_outgroup, peak]
    sd_sis_out = ortholog_db.at[idx_sister_outgroup, sd]

    # RRT formulas
    rel_rate_species = (peak_sp_out + peak_sp_sis - peak_sis_out) / 2.0 # also called k_AO
    rel_rate_sister = (peak_sp_sis + peak_sis_out - peak_sp_out) / 2.0 # also called k_BO
    # Error propagation rules
    rel_rate_species_sd = sqrt(pow(sd_sp_out, 2) + pow(sd_sp_sis, 2) + pow(sd_sis_out, 2)) / 2.0 # also called k_AO_sd
    rel_rate_sister_sd = sqrt(pow(sd_sp_sis, 2) + pow(sd_sis_out, 2) + pow(sd_sp_out, 2)) / 2.0 # also called k_BO_sd

    return rel_rate_species, rel_rate_species_sd, rel_rate_sister, rel_rate_sister_sd


def compute_corrected_ks_species_sister(rel_rate_species, rel_rate_species_sd):
    """
    Performs substitution rate-adjustment on the original peak of the ortholog distribution between
    the focal species (species A) and sister species (species B).

    :param rel_rate_species: branch-specific Ks contribution of the focal species, namely the synonymous substitutions occurred 
        in the focal species since the speciation event with the sister species; also called K_OA
    :param rel_rate_species_sd: standard deviation associated to rel_rate_species; also called k_AO_sd
    :return: corrected_peak, rate-corrected ortholog Ks distribution peak of focal species and sister species
    :return: corrected_peak_sd, standard deviation associated to corrected_peak
    """
    corrected_peak = rel_rate_species * 2
    corrected_peak_sd = sqrt(2) * rel_rate_species_sd

    return corrected_peak, corrected_peak_sd
