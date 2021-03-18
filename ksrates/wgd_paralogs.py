import os
import sys
import logging
import datetime
import ksrates.fc_configfile as fcConf
import ksrates.fc_check_input as fcCheck
import ksrates.fc_wgd as fc_wgd
from ksrates.utils import init_logging


def wgd_paralogs(config_file, n_threads):
    # INPUT
    # Get parameters and FASTA files from configuration file
    config = fcConf.Configuration(config_file)
    species = config.get_species()
    init_logging(f"Paralog wgd analysis for species [{species}]", config.get_logging_level())

    latin_names = config.get_latin_names()
    max_gene_family_size = config.get_max_gene_family_size()
    paranome = config.get_paranome()
    colinearity = config.get_colinearity()

    logging.info(f"Checking if sequence data files exist and if sequence IDs are compatible with wgd pipeline...")
    # Will exit if FASTA or GFF files are missing or empty or if GFF feature/attribute are missing
    trigger_exit = False

    if colinearity:  # if colinearity analysis is required, load related parameters
        gff = config.get_gff(species)
        if fcCheck.check_file_nonexistent_or_empty(gff, "GFF file"):
            trigger_exit = True

        gff_feature = config.get_feature()
        gff_gene_attribute = config.get_attribute()
        if gff_feature == "":
            logging.error("No GFF attribute provided in configuration file. Will exit.")
            trigger_exit = True
        if gff_gene_attribute == "":
            logging.error("No GFF feature provided in configuration file. Will exit.")
            trigger_exit = True

    # Checking if FASTA file exists and if sequence IDs are compatible with wgd pipeline (paml)
    fasta_names_dict = config.get_fasta_dict()
    species_fasta_file = config.get_fasta_name(fasta_names_dict, species)
    if fcCheck.check_file_nonexistent_or_empty(species_fasta_file, "FASTA file"):  # if missing/empty
        trigger_exit = True
    else: # If FASTA file exists, check for ID compatibility
        if colinearity:
            fcCheck.check_IDs(species_fasta_file, latin_names[species], gff)
        else:
            fcCheck.check_IDs(species_fasta_file, latin_names[species])

    if trigger_exit:
        logging.error("Please add the missing information to the configuration file and rerun the analysis. Exiting.")
        sys.exit(1)
    logging.info("Completed")

    # Creating folder for output files of wgd paralog pipeline
    paralog_dists_dir = os.path.join("paralog_distributions", "")
    if not os.path.isdir(paralog_dists_dir):
        logging.info(f"Creating directory [{paralog_dists_dir}]")
        os.makedirs(paralog_dists_dir)

    # -----------------------------------------------------------------------------

    # ESTIMATING PARANOME Ks VALUES
    logging.info("Running wgd paralog Ks pipeline...")
    fc_wgd.ks_paralogs(species, species_fasta_file, max_gene_family_size=max_gene_family_size, base_dir=paralog_dists_dir, n_threads=n_threads)

    # EXTRACTING COLINEARITY/SYNTENY ANCHOR PAIRS Ks VALUES
    if colinearity:
        logging.info('---')
        logging.info("Running wgd colinearity Ks pipeline...")
        fc_wgd.ks_colinearity(species, gff, base_dir=paralog_dists_dir, gff_feature=gff_feature,
                            gff_gene_attribute=gff_gene_attribute, n_threads=n_threads)

    logging.info(datetime.datetime.today().ctime())
    logging.info("Done")
