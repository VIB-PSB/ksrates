import os
import logging
import datetime
import ksrates.fc_configfile as fcConf
import ksrates.fc_check_input as fcCheck
import ksrates.fc_wgd as fc_wgd
from ksrates.utils import init_logging


def wgd_orthologs(config_file, species_one, species_two, n_threads):
    # INPUT
    species_pair = sorted([species_one, species_two], key=str.casefold)
    species1, species2 = species_pair[0], species_pair[1] # sorted!

    config = fcConf.Configuration(config_file)
    init_logging(f"Ortholog wgd analysis for species pair [{species1} - {species2}]", config.get_logging_level())

    # Get parameters and FASTA files from configuration file
    newick_tree = config.get_newick_tree()
    latin_names = config.get_latin_names(newick_tree)

    fasta_names_dict = config.get_fasta_dict()
    species1_fasta_file = config.get_fasta_name(fasta_names_dict, species1)
    fcCheck.check_inputfile(species1_fasta_file, "FASTA file")
    fcCheck.check_IDs(species1_fasta_file, latin_names[species1])

    species2_fasta_file = config.get_fasta_name(fasta_names_dict, species2)
    fcCheck.check_inputfile(species2_fasta_file, "FASTA file")
    fcCheck.check_IDs(species2_fasta_file, latin_names[species2])

    # Creating folder for output files of wgd ortholog pipeline.
    # Note: since in Nextflow mode there are multiple wgdOrtholog processes running in parallel,
    # this "if-try-except" prevents that almost-simultaneous checks rise an error: a slower process 
    # would rise an error if in the meanwhile a faster process had already created the folder.
    ortholog_dists_dir = os.path.join("ortholog_distributions", "")
    if not os.path.exists(ortholog_dists_dir):
        logging.info(f"Creating directory {ortholog_dists_dir}")
        os.makedirs(ortholog_dists_dir, exist_ok=True)

    # -----------------------------------------------------------------------------

    # ESTIMATING ORTHOLOG Ks VALUES
    logging.info("Running wgd ortholog Ks pipeline...")
    fc_wgd.ks_orthologs(species1, species2, species1_fasta_file, species2_fasta_file, base_dir=ortholog_dists_dir,
                        n_threads=n_threads)

    logging.info(datetime.datetime.today().ctime())
    logging.info("Done")
