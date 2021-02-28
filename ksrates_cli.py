import click
import logging
from sys import argv

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.version_option("0.1", prog_name="ksrates", help="Print version number.")
def cli():
    """
    Welcome to ksrates!
    """


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Generates configuration file.")
@click.argument('filename')
def generate_config(filename):
    """
    Generates the configuration file for the rate-adjustment.
    The configuration file name is given by argument FILENAME.

    FILENAME: configuration file name
    """
    from ksrates.generate_configfile import generate_configfile
    generate_configfile(filename)
    

@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Initializes rate-adjustment.")
@click.argument('config_file', type=click.Path(exists=True))
@click.option("-n", "--nextflow", is_flag=True, help="Flag for Nextflow pipeline (Default: False)")
def init(config_file, nextflow):
    """
    Initializes rate-adjustment from CONFIG_FILE.

    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    """
    from ksrates.setup_correction import setup_correction
    click.format_filename(config_file)
    setup_correction(config_file, nextflow)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Performs paralog Ks estimation.")
@click.argument('config_file', type=click.Path(exists=True))
@click.option("--n-threads", type=int, default=4, help="Number of threads (default: 4)")
def paralogs_ks(config_file, n_threads):
    """
    Performs paralog Ks estimation for the focal species through wgd.

    Takes parameters from CONFIG_FILE. 

    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    """
    from ksrates.wgd_paralogs import wgd_paralogs
    click.format_filename(config_file)
    wgd_paralogs(config_file, n_threads)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Performs ortholog Ks estimation.")
@click.argument('config_file', type=click.Path(exists=True))
@click.argument("species1")
@click.argument("species2")
@click.option("--n-threads", type=int, default=4, help="Number of threads (default: 4)")
def orthologs_ks(config_file, species1, species2, n_threads):
    """
    Performs ortholog Ks estimation for SPECIES1 and SPECIES2 through wgd.

    Takes parameters from CONFIG_FILE. 

    \b
    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    SPECIES1: first of the two species involved in the ortholog Ks estimation
    SPECIES2: second of the two species involved in the ortholog Ks estimation
    """
    from ksrates.wgd_orthologs import wgd_orthologs
    click.format_filename(config_file)
    wgd_orthologs(config_file, species1, species2, n_threads)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, 
             short_help="Computes ortholog divergence times Ks estimates.")
@click.argument('config_file', type=click.Path(exists=True))
@click.option('--ortholog-pairs', type=click.Path(exists=True), help="User-defined path to file containing the ortholog pairs with missing ortholog Ks peak in database (default: correction_analysis/species/ortholog_pairs_species.tsv)")
def orthologs_analysis(config_file, ortholog_pairs):
    """
    Computes ortholog Ks distribution mode (or median) and updates the ortholog databases.

    Takes parameters from CONFIG_FILE. 

    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    """
    from ksrates.compute_peaks import compute_peaks
    click.format_filename(config_file)
    if ortholog_pairs:
        click.format_filename(ortholog_pairs)
    compute_peaks(config_file, ortholog_pairs)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, 
             short_help="Performs ortholog substitution rate-adjustment.")
@click.argument('config_file', type=click.Path(exists=True))
@click.option("--trios", type=click.Path(exists=True), help="User-defined path to file containing the ortholog trios (default: correction_analysis/species/orthologs_trios_species.tsv)")
def orthologs_adjustment(config_file, trios):
    """
    Performs substitution rate-adjustment relative to the focal species.

    Takes parameters from CONFIG_FILE. 

    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    """
    from ksrates.correct import correct
    click.format_filename(config_file)
    if trios:
        click.format_filename(trios)
    correct(config_file, trios)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Generates rate-adjusted mixed Ks plot.")
@click.argument('config_file', type=click.Path(exists=True))
@click.option("--correction-table", type=click.Path(exists=True), help="User-defined path to file containing adjustment results (default: correction_analysis/species/correction_table_species.tsv)")
@click.option("--paranome-table", type=click.Path(exists=True), help="User-defined path to file containing paranome Ks (default: paralog_distributions/wgd_species/species.ks.tsv)")
@click.option("--anchors-table", type=click.Path(exists=True), help="User-defined path to file containing anchor pair Ks (default: paralog_distribution/wgd_species/species.ks_anchors.tsv)")
def plot_paralogs(config_file, correction_table, paranome_table, anchors_table):
    """
    Plots rate-adjusted mixed paralog-ortholog Ks distribution.
        
    Takes parameters from CONFIG_FILE. 

    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    """
    from ksrates.plot_paralogs import plot_paralogs_distr
    click.format_filename(config_file)
    if correction_table:
        click.format_filename(correction_table)
    if paranome_table:
        click.format_filename(paranome_table)
    if anchors_table:
        click.format_filename(anchors_table)
    plot_paralogs_distr(config_file, correction_table, paranome_table, anchors_table)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Generates phylogram with Ks-unit branch lengths.")
@click.argument('config_file', type=click.Path(exists=True))
@click.option("--correction-table", type=click.Path(exists=True), help="User-defined path to file containing adjustment results (default: correction_analysis/species/correction_table_species.tsv)")
@click.option("-n", "--nextflow", is_flag=True, help="Flag for Nextflow pipeline (Default: False)")
def plot_tree(config_file, correction_table, nextflow):
    """
    Generates a phylogram of the input dataset with branch lengths set to\
    Ks distances estimated from ortholog KS distributions.

    Takes parameters from CONFIG_FILE. 

    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    """
    from ksrates.plot_tree import plot_tree_rates
    click.format_filename(config_file)
    if correction_table:
        click.format_filename(correction_table)
    plot_tree_rates(config_file, correction_table, nextflow)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Generates ortholog Ks distributions plot.")
@click.argument('config_file', type=click.Path(exists=True))
@click.option("--trios", type=click.Path(exists=True), help="User-defined path to file containing the ortholog trios (default: correction_analysis/species/orthologs_trios_species.tsv)")
def plot_orthologs(config_file, trios):
    """
    Plots ortholog Ks distributions used for rate-adjustment.

    Takes parameters from CONFIG_FILE. 

    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    """
    from ksrates.plot_orthologs import plot_orthologs_distr
    click.format_filename(config_file)
    if trios:
        click.format_filename(trios)
    plot_orthologs_distr(config_file, trios)


@cli.command(context_settings={'help_option_names': ['-h', '--help']}, short_help="Detects WGD signatures in paralog Ks distribution.")
@click.argument('config_file', type=click.Path(exists=True))
@click.option("--paranome-table", type=click.Path(exists=True), help="User-defined path to file containing paranome Ks (default: paralog_distributions/wgd_species/species.ks.tsv)")
@click.option("--anchors-table", type=click.Path(exists=True), help="User-defined path to file containing anchor pair Ks (default: paralog_distribution/wgd_species/species.ks_anchors.tsv)")
@click.option("--correction-table", type=click.Path(exists=True), help="User-defined path to file containing adjustment results (default: correction_analysis/species/correction_table_species.tsv)")
@click.option("--anchorpoints", type=click.Path(exists=True), help="User-defined path to i-ADHoRe file anchorpoints.txt (default: paralog_distributions/wgd_species/species_i-adhore/anchorpoints.txt)")
@click.option("--multiplicons", type=click.Path(exists=True), help="User-defined path to i-ADHoRe file multiplicons.txt (default: paralog_distributions/wgd_species/species_i-adhore/multiplicons.txt)")
@click.option("--segments", type=click.Path(exists=True), help="User-defined path to i-ADHoRe file segments.txt (default: paralog_distributions/wgd_species/species_i-adhore/segments.txt)")
@click.option("--list-elements", type=click.Path(exists=True), help="User-defined path to i-ADHoRe file list_elements.txt (default: paralog_distributions/wgd_species/species_i-adhore/list_elements.txt)")
@click.option("--multiplicon-pairs", type=click.Path(exists=True), help="User-defined path to i-ADHoRe file multiplicons_pairs.txt (default: paralog_distributions/wgd_species/species_i-adhore/multiplicons_pairs.txt)")
def paralogs_analyses(config_file, paranome_table, anchors_table, correction_table, anchorpoints, multiplicons, segments, list_elements, multiplicon_pairs):
    """
    Reconstructs potential WGD peaks in the paralog Ks distributions.

    Performs anchor Ks clustering if "collinearity" analysis is switched on, otherwise performs exponential-lognormal mixture model on paranome.
    
    If extra methods are asked through the expert configuration file, performs all methods available for the analysis type(s) selected.

    Takes parameters from CONFIG_FILE. 

    CONFIG_FILE: configuration file to set up the rate-adjustment relative to the focal species
    """
    from ksrates.paralogs_analyses import paralogs_analyses_methods
    click.format_filename(config_file)
    if paranome_table:
        click.format_filename(paranome_table)
    if anchors_table:
        click.format_filename(anchors_table)
    if correction_table:
        click.format_filename(correction_table)
    if anchorpoints:
        click.format_filename(anchorpoints)
    if multiplicons:
        click.format_filename(multiplicons)
    if segments:
        click.format_filename(segments)
    if list_elements:
        click.format_filename(list_elements)
    if multiplicon_pairs:
        click.format_filename(multiplicon_pairs)
    paralogs_analyses_methods(config_file, paranome_table, anchors_table, correction_table, 
                  anchorpoints, multiplicons, segments, list_elements, multiplicon_pairs)



# For debugging
# Syntax: python3 ksrates_cli.py [command] [args]
if __name__ == "__main__":
    cli(argv[1:])