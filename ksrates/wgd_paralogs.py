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
    fasta_names_dict = config.get_fasta_dict()
    species_fasta_file = config.get_fasta_name(fasta_names_dict, species)
    fcCheck.check_inputfile(species_fasta_file, "FASTA file")
    max_gene_family_size = config.get_max_gene_family_size()

    paranome = config.get_paranome()
    colinearity = config.get_colinearity()

    if colinearity:  # if colinearity analysis is required, load extra parameters
        gff_names_dict = config.get_gff_dict()
        species_gff_file = config.get_gff_name(gff_names_dict, species)
        fcCheck.check_inputfile(species_gff_file, "GFF file")

        gff_feature = config.get_feature()
        gff_gene_attribute = config.get_attribute()
        if gff_feature == "":
            logging.error("No GFF attribute provided in configuration file.")
        if gff_gene_attribute == "":
            logging.error("No GFF feature provided in configuration file.")
        if gff_feature == "" or gff_gene_attribute == "":
            logging.error("Will exit the colinearity analysis.")
            sys.exit(1)
    # Checking if IDs in FASTA (and in GFF if applicable) are compatible with wgd pipeline (paml)
    logging.info("")
    if colinearity:
        fcCheck.check_IDs(species_fasta_file, latin_names[species], species_gff_file)
    else:
        fcCheck.check_IDs(species_fasta_file, latin_names[species])


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
        fc_wgd.ks_colinearity(species, species_gff_file, base_dir=paralog_dists_dir, gff_feature=gff_feature,
                            gff_gene_attribute=gff_gene_attribute, n_threads=n_threads)

    logging.info(datetime.datetime.today().ctime())
    logging.info("Done")
