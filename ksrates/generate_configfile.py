import sys
import configparser
import logging
from ksrates.utils import init_logging


def generate_configfile(configfile_name):

    cfgfile = open(configfile_name, "w+")

    Config = configparser.ConfigParser(allow_no_value=True)

    Config.add_section("SPECIES")
    Config.set("SPECIES", "focal_species", "")
    Config.set("SPECIES", "# informal name of the species that will be used to perform the rate-adjustment\n")

    Config.set("SPECIES", "newick_tree", "();")
    Config.set("SPECIES", "# input phylogenetic tree in Newick format; use the informal names of the species\n")

    Config.set("SPECIES", "latin_names", "")
    Config.set("SPECIES", "# informal names associated to their scientific names through a colon and separated by comma\n")

    Config.set("SPECIES", "fasta_filenames", "")
    Config.set("SPECIES", "gff_filename", "")
    Config.set("SPECIES", "# informal names associated to their FASTA or GFF filenames/paths through a colon and separated by commas\n")

    Config.set("SPECIES", "peak_database_path", "ortholog_peak_db.tsv")
    Config.set("SPECIES", "ks_list_database_path", "ortholog_ks_list_db.tsv")
    Config.set("SPECIES", "# filenames/paths of the ortholog data databases\n")

    Config.add_section("ANALYSIS SETTING")

    Config.set("ANALYSIS SETTING", "paranome", "yes")
    Config.set("ANALYSIS SETTING", "collinearity", "no")
    Config.set("ANALYSIS SETTING", "# analysis type for paralog data; allowed values: 'yes' and 'no'\n")

    Config.set("ANALYSIS SETTING", "gff_feature", "")
    Config.set("ANALYSIS SETTING", "# keyword to parse the sequence type from the GFF file (column 3); can be 'gene', 'mRNA'...\n")

    Config.set("ANALYSIS SETTING", "gff_attribute", "")
    Config.set("ANALYSIS SETTING", "# keyword to parse gene ID from the GFF file (column 9); can be 'ID', 'Name'...\n")

    Config.set("ANALYSIS SETTING", "max_number_outgroups", "4")
    Config.set("ANALYSIS SETTING", "# maximum number of outspecies/trios selected to correct each divergent species pair (default: 4)\n")

    Config.set("ANALYSIS SETTING", "consensus_mode_for_multiple_outgroups", "mean among outgroups")
    Config.set("ANALYSIS SETTING", "# allowed values: 'mean among outgroups' or 'best outgroup' (default: 'mean among outgroups')\n")

    Config.add_section("PARAMETERS")

    Config.set("PARAMETERS", "x_axis_max_limit_paralogs_plot", "5")
    Config.set("PARAMETERS", "# highest value of the x axis in the mixed distribution plot (default: 5)\n")

    Config.set("PARAMETERS", "bin_width_paralogs", "0.1")
    Config.set("PARAMETERS", "# bin width in paralog Ks histograms (default: 0.1, ten bins per unit)\n")

    Config.set("PARAMETERS", "y_axis_max_limit_paralogs_plot", "None")
    Config.set("PARAMETERS", "# highest value of the y axis in the mixed distribution plot  (default: None)\n")

    Config.set("PARAMETERS", "num_bootstrap_iterations", "200")
    Config.set("PARAMETERS", "# number of bootstrap iterations for ortholog peak estimate\n")

    Config.set("PARAMETERS", "divergence_colors", "Red, MediumBlue, DarkGoldenrod, ForestGreen, HotPink, DarkCyan, SaddleBrown, Black")
    Config.set("PARAMETERS", "# color of the divergence lines drawn in correspondence of the ortholog peaks")
    Config.set("PARAMETERS", "# use color names/codes separated by comma and use at least as many colors as the number of divergence nodes\n")

    Config.set("PARAMETERS", "x_axis_max_limit_orthologs_plots", "5")
    Config.set("PARAMETERS", "# highest value of the x axis in the ortholog distribution plots (default: 5)\n")

    Config.set("PARAMETERS", "bin_width_orthologs", "0.1")
    Config.set("PARAMETERS", "# bin width in ortholog Ks histograms (default: 0.1, ten bins per unit)\n")

    Config.set("PARAMETERS", "max_ks_paralogs", "5")
    Config.set("PARAMETERS", "# maximum paralog Ks value accepted from Ks data table (default: 5)\n")

    Config.set("PARAMETERS", "max_ks_orthologs", "10")
    Config.set("PARAMETERS", "# maximum ortholog Ks value accepted from Ks data table (default: 10)")

    logging.basicConfig(format='%(levelname)s\t%(message)s', level="INFO", stream=sys.stdout)
    logging.info(f"Configuration file [{configfile_name}] not found: it will be now generated")
    logging.info(f"Please fill in the required input parameters:")
    logging.info(f"species, newick_tree, latin_names, fasta_filenames and if applicable gff_filename, gff_feature and gff_attribute")
    Config.write(cfgfile)
    cfgfile.close()
