import configparser
import os
from ete3 import Tree
import ksrates.fc_check_input as fcCheck
from matplotlib.colors import is_color_like
import logging
import sys
from ast import literal_eval

class Configuration:
    
    def __init__(self, config_path, expert_config_file):
        """
        Initializes the configuration file and the expert configuration file.
        This latter is always named "config_expert.txt", the code looks for it
        and if it is found then takes the customized expert parameters from there.
        """
        fcCheck.check_inputfile(config_path, "Configuration file")
        # Configuration file
        self.config = configparser.ConfigParser()
        self.config.read(config_path)

        # Expert configuration file
        self.expert_config = configparser.ConfigParser()
        # If there is no user-defined expert config file given through "--expert" in the command line,
        # variable "expert_config_file" was set to an empty string in the CLI code block
        if expert_config_file == "":
            # If there is in the launching folder an expert config file called with the default name "config_expert.txt", fallback to this latter
            if os.path.exists("config_expert.txt"):
                self.expert_config.read("config_expert.txt")
            else:
                # Else set the variable to None (scripts will use default expert parameters)
                self.expert_config = None
        # Else if there is a user-defined expert config file given through "--expert",
        # variable expert_config_file has already been set to such given filename (so, it's not an empty string)
        else:
            # If the user-defined file exists, read it
            if os.path.exists(expert_config_file):
                self.expert_config.read(expert_config_file)
            else:
                # Else if the user-defined file doesn't exists (e.g. misspelled), print an error and exit
                # (This case is actually already caught by the CLI option definition of "--expert")
                logging.error(f"User-defined expert configuration file {expert_config_file} not found:")
                logging.error("please check the input path after the '--expert' parameter in the command line and rerun the analysis")
                sys.exit(1)
    

    def _get_clean_dict(self, dict_like_string, parameter):
        """This method reads a dictionary-like field from the configuration file \\
        (FASTA filenames, GFF filenames) as a python dictionary \\
        despite the absence of braces and quotes and the presence of whitespaces\\
        (spaces, line breaks) around the pairs and next to the colon symbol.

        Some malformed types are not tolerated (no key; one word with no colon)
        i.e. { :filename.fasta, species, filename.fasta}.
        One malformed type is instead tolerated (key with colon and no value) because it will be later\\
        filled in by a fallback value.
        i.e. {species: } will become {species: filename.fasta}

        SINGLE-LINE OPTION to build the dictionary:
        (limitation: spaces between the colon symbol and the species names raise error!)
        fasta_names_dict = dict(f2.split(':') for f2 in (f1.strip() for f1 in fasta_names_list))

        :param dict_like_string: dictionary-like string from config file; it's not a Python dictionary
        :param parameter: reference to the dictionary that is malformed (either "FASTA" or "GFF")
        :return: a Python dictionary
        """
        clean_list = []
        for key_value_string in dict_like_string.split(","):
            key_value_pair = key_value_string.strip().split(":")
            # The following check catches if the user has left the informal species name without both colon and latin name.
            # However this check doesn't work with the FASTA filename system because it would assign to Elaeis fasta filename "elaeis", 
            # which is never correct.
            # I suggest to remove this check and to force the user to insert all latin names, which are needed to do a good check in the
            # ortholog database and for further pipeline steps anyways
            # if len(key_value_pair) < 2:
            #     logging.warning(f"No mapping found for key [{key_value_pair[0]}] in config parameter [{parameter_name}], "
            #           "will use key as value.")
            #     key_value_pair.append(key_value_pair[0])

            # Malformed type { :filename.fasta, species, filename.fasta} --> exit
            if key_value_pair[0] == "" or len(key_value_pair) == 1:
                logging.error(f"Malformed {parameter} dictionary: {key_value_pair}.")
                logging.error("Exiting")
                sys.exit(1)
            # It covers the malformed type {species: } --> convert into {species: ""} and also the good (expected) case {species:filename.fasta}
            else:
                clean_list.append((key_value_pair[0].strip(), key_value_pair[1].strip()))
        clean_dict = dict(clean_list)
        return clean_dict

    def _get_clean_dict_stringent(self, dict_like_string, parameter):
        """
        Converts the dictionary-like structure of latin_names and divergence_colors
        into a Python dictionary. If the dictionary-like structure is malformed or
        has a missing value, it exits.

        This function is more stringent than "_get_clean_dict" because ALL malformed types are not tolerated
        e.g. {species: , :filename.fasta, species, filename.fasta}

        :param dict_like_string: dictionary-like string from config file; it's not a Python dictionary
        :param parameter: config file filed (either "latin_names" or "divergence_colors")
        :return: a Python dictionary
        """
        clean_list = []
        for key_value_string in dict_like_string.split(","):
            key_value_pair = key_value_string.strip().split(":")
            # if one latin name not found, exit the script
            if len(key_value_pair) == 1 or "" in key_value_pair:
                logging.error(f"Malformed {parameter} dictionary in configuration file: {key_value_pair}.")
                logging.error(f"Exiting")
                sys.exit(1)
            else:
                clean_list.append((key_value_pair[0].strip(), key_value_pair[1].strip()))
        clean_dict = dict(clean_list)
        return clean_dict

    def get_species(self):
        """
        Gets the config file field of the name of the focal species.
        The focal species is one of the leaves of the input tree and is the one\\
        whose paralog distribution is plotted and which the rate-adjustment will be relative to. 

        :return species: informal species name 
        """
        species = self.config.get("SPECIES", "focal_species")
        if species == "":
            logging.error('Field "focal_species" in configuration file is empty, please fill in')
            sys.exit(1)
        elif len(species.split()) != 1 or "_" in species:
            logging.error(f'Field "focal_species" [{species}] should be a short name and must not contain any spaces or underscores, please change accordingly')
            sys.exit(1)
        return species

    def get_newick_tree(self):
        """
        Gets the config file field of the Newick tree.
        Checks and exits if the species' names in the Newick tree contain illegal characters (underscore or spaces).

        :return tree_string: the tree object by ete3
        """
        tree_string = self.config.get("SPECIES", "newick_tree")
        if not (tree_string.endswith(';')):
            tree_string += ";"
        if tree_string == "();" or tree_string == ";":
            logging.error('Field "newick_tree" in configuration file is empty, please fill in')
            sys.exit(1)
        try:
            tree = Tree(tree_string)
        except Exception:
            logging.error('Unrecognized format for field "newick_tree" in configuration file (for example, parentheses do not match)')
            sys.exit(1)

        # Check if species' informal names contain illegal characters (underscore or spaces)
        species_illegal_char=[]
        for informal_name in tree.get_leaf_names():
            if "_" in informal_name or " " in informal_name:
                species_illegal_char.append(informal_name)
        if len(species_illegal_char) != 0:
            logging.error(f"Informal species' names must not contain any spaces or underscores. Please change the following names in the configuration file:")
            for informal_name in species_illegal_char:
                logging.error(f"- {informal_name}")
            sys.exit(1)
        return tree


    def check_complete_latin_names_dict(self, dictionary):
        """
        Checks if a dictionary field from latin_names contains all the species present in the Newick tree. 
        If one or more species are missing, it exits.

        :param dictionary: the dictionary coming from latin_names
        """
        all_leaves = []
        for leaf in self.get_newick_tree().get_leaves():
            all_leaves.append(leaf.name)
        missing_species = list(set.difference(set(all_leaves), set(dictionary.keys())))
        if len(missing_species) != 0:
            if len(missing_species) == 1:
                logging.error(f'The following species is missing from the "latin_names" configuration file field:')
            else:
                logging.error(f'The following species are missing from the "latin_names" configuration file field:')
            for missing_name in missing_species:
                logging.error(f" - {missing_name}")
            
            logging.error(f"Please add the missing information and restart the analysis.")
            logging.error("Exiting.")
            sys.exit(1)

    def get_latin_names(self):
        """
        Gets the config file field of the dictionary that associates the informal species name to its latin (scientific) name.                                      

        :return latin_names_dict: python dictionary
        """
        latin_names = self.config.get("SPECIES", "latin_names")
        if latin_names != "":
            latin_names_dict = self._get_clean_dict_stringent(latin_names, "latin_names")
        else:
            logging.error('Configuration file field "latin_names" is empty, please fill in and restart the analysis.')
            logging.error("Exiting.")
            sys.exit(1)
        # Check if latin_names contains all the species present in the Newick tree; if not, exits
        self.check_complete_latin_names_dict(latin_names_dict)
        return latin_names_dict

    def get_ortho_db(self):
        """
        Gets the config file field of the peak ortholog database path.

        :return db_path: path to the ortholog peak database
        """
        db_path = self.config.get("SPECIES", "peak_database_path", fallback="ortholog_peak_db.tsv")
        if not db_path:
            db_path = "ortholog_peak_db.tsv"
        return db_path
    
    def get_ks_db(self):
        """
        Gets the config file field of the ortholog Ks list database path.  

        :return ks_list_db_path: path to the ortholog Ks list database
        """
        ks_list_db_path = self.config.get("SPECIES", "ks_list_database_path", fallback="ortholog_ks_list_db.tsv")
        if not ks_list_db_path:
            ks_list_db_path = "ortholog_ks_db.tsv"
        return ks_list_db_path

    def get_fasta_dict(self):
        """
        Gets the config file field of the dictionary that associates the informal species names to the their FASTA files.

        :return fasta_names_dict: python dictionary
        """
        fasta_names_string = self.config.get("SPECIES", "fasta_filenames")
        if fasta_names_string != "":
            fasta_names_dict = self._get_clean_dict(fasta_names_string, "FASTA file")
        else:
            logging.warning("Configuration file field [fasta_filenames] is empty.")
            fasta_names_dict = {}
        return fasta_names_dict

    def get_fasta_name(self, fasta_dict, species):
        """
        Gets the path to the FASTA file of a species from the dictionary.
        If the value is an empty string, it applies the fallback filename "species + .fasta".

        :param fasta_dict: Python dictionary that associates each informal species name to the path of its FASTA file 
        :param species: the species informal name
        :return: the FASTA file path
        """
        if species in fasta_dict:
            if fasta_dict[species] != "": # if the fasta filename is an acceptable string (not empty)
                fasta = fasta_dict[species]
            else:
                logging.warning(f"FASTA filename for {species} not found in configuration file; assuming default one ({species}.fasta)")
                fasta = f"{species}.fasta"  # fallback name
        else:  # if species is missing from fasta_dict
            logging.warning(f"FASTA filename for {species} not found in configuration file; assuming default one ({species}.fasta)")
            fasta = f"{species}.fasta"  # fallback name 
        return fasta

    def get_gff(self, species):
        """
        Gets the config file field of the focal species's GFF file.

        :species: focal species
        :return gff: GFF filename
        """
        gff = self.config.get("SPECIES", "gff_filename")
        if gff == "":
            logging.warning(f"GFF filename for focal species [{species}] not found in configuration file; assuming default one ({species}.gff)")
            gff = f"{species}.gff"   # fallback name
        return gff

    def get_feature(self):
        """
        Assures that the entered feature word for parsing the GFF file matches
        the right lowercase/uppercase standards for GFF format.\\
        Rises a warning if the given feature is not among the commonly used for protein-coding genes.\\
        Rises an error and exits if the GFF feature field is left empty.

        :return f: GFF file feature corrected for letter case sensitivity
        :return feature: GFF file feature as it was provided by the user
        """
        feature = self.config.get("ANALYSIS SETTING", "gff_feature")
        if feature == "":
            return "" # the empty string will be used to exit the colinearity analysis through the sys command in the main script
        # Setting the correct lower/upper cases in the most common terms:
        f = feature.lower()
        if f == "gene": return f # standard feature "gene" is lowercase, just return it 
        elif f == "mrna":
            f = "mRNA"
            return f
        elif f == "cds":
            f = "CDS"
            return f
        elif f == "cdna_match":
            f = "cDNA_match"
            return f
        elif f == "transcript": return f
        elif f == "region": return f
        elif f == "exon": return f
        else:
            logging.warning(f"The provided GFF feature field [{feature}] is not a common choice (mRNA, gene) for protein-coding genes.")
            return feature

    def get_attribute(self):
        """
        Assures that the entered attribute word for parsing the GFF file matches 
        the right lowercase/uppercase standards for GFF format.\\
        Rises a warning if the given attribute is not among the commonly used for protein-coding genes.\\
        Rises an error and exits if the GFF attribute field is left empty.

        :return a: GFF file attribute corrected for letter case sensitivity
        :return attribute: GFF file attribute as it was provided by the user
        """
        attribute = self.config.get("ANALYSIS SETTING", "gff_attribute")
        if attribute == "":
            return "" # the empty string will be used to exit the colinearity analysis through the sys command in the main script
        # Setting the correct lower/upper cases in the most common terms:
        a = attribute.lower()
        if a == "id":
            a = "ID"
            return a
        elif a == "name":
            a = "Name"
            return a
        elif a == "parent":
            a = "Parent"
            return a
        else:
            logging.warning(f"The provided GFF attribute field [{attribute}] is not a commonly used one (ID, Name, Parent).")
            return attribute

    def get_max_num_outspecies(self):
        """
        Gets the config file field specifying the maximum number of outgroup species to be used when adjusting a divergence.

        :return max_outspecies: integer
        """
        max_outspecies = self.config.get("ANALYSIS SETTING", "max_number_outgroups")
        return max_outspecies

    def get_paranome(self):
        """
        Gets the config file field specifying if the mixed distribution will show the whole-paranome of the focal species or not.

        :return boolean: flags whether the paranome analysis is required
        """
        paranome = self.config.get("ANALYSIS SETTING", "paranome").lower()
        if paranome == "yes":
            return True
        elif paranome == "no":
            return False
        else:
            logging.error('Unrecognized "paranome" parameter in configuration file; please choose between "yes" and "no"')
            sys.exit(1)

    def get_colinearity(self):
        """
        Gets the config file field specifying if the mixed distribution will show the anchor pairs of the focal species or not.

        :return boolean: flags whether the colinearity analysis is required
        """
        colinearity = self.config.get("ANALYSIS SETTING", "collinearity").lower()

        if colinearity == "yes":
            return True
        elif colinearity == "no":
            return False
        else:
            logging.error('Unrecognized "collinearity" parameter in configuration file; please choose between "yes" and "no"')
            sys.exit(1)

    def get_reciprocal_retention(self):
        """
        Gets the config file field specifying if the mixed distribution will show the reciprocally retained gene families of the focal species or not.

        :return boolean: boolean
        """
        reciprocal_retention = self.config.get("ANALYSIS SETTING", "reciprocal_retention").lower()

        if reciprocal_retention == "yes":
            return True
        elif reciprocal_retention == "no":
            return False
        else:
            logging.error('Unrecognized "reciprocal_retention" parameter in configuration file; please choose between "yes" and "no"')
            sys.exit(1)

    def get_consensus_peak_for_multiple_outgroups(self):
        """
        Gets the config file field of the consensus method for how to deal with multiple adjustments for the same divergence.
        Checks if the user has misspelled or left empty the field in the configuration file for
        the choice on how to deal with multiple outgroups when adjusting a divergent pair.

        :return consensus_peak_for_multiple_outgroups: a valid string for the consensus method
        """
        consensus_peak_for_multiple_outgroups = self.config.get("ANALYSIS SETTING", "consensus_mode_for_multiple_outgroups")

        if consensus_peak_for_multiple_outgroups == "mean among outgroups" or consensus_peak_for_multiple_outgroups == "best outgroup":
            return consensus_peak_for_multiple_outgroups
        else:
            logging.warning("Unrecognized choice in 'consensus_mode_for_multiple_outgroups' filed in configuration file.")
            logging.warning("  Please choose between 'mean among outgroups' and 'best outgroup'")
            logging.warning("  The default option will be executed ('mean among outgroups').")
            consensus_peak_for_multiple_outgroups = "mean among outgroups"
            return consensus_peak_for_multiple_outgroups

    
    def get_max_ks_ortho(self):
        """
        Gets the config file field specifying the maximum ortholog Ks value accepted for the analysis.

        :return max_ks_ortho: integer
        """
        max_ks_ortho = float(self.config.get("PARAMETERS", "max_ks_orthologs"))
        return max_ks_ortho

    def get_max_ks_para(self):
        """
        Gets the config file field specifying the maximum paralog Ks value accepted for the analysis.

        :return max_ks_para: integer
        """
        max_ks_para = float(self.config.get("PARAMETERS", "max_ks_paralogs"))
        return max_ks_para

    def get_num_iteration(self):
        """
        Gets the config file field specifying the number of bootstrap iterations for the ortholog peak estimate.

        :return n_inter: integer
        """
        n_iter = int(self.config.get("PARAMETERS", "num_bootstrap_iterations"))
        return n_iter

    def get_bin_width_para(self):
        """
        Gets the config file field specifying the width of the paralog histogram bins.

        :return bin_width_para: number (float or integer) for bin width in paralog Ks histogram
        """
        bin_width_para = float(self.config.get("PARAMETERS", "bin_width_paralogs"))
        return bin_width_para

    def get_bin_width_ortho(self):
        """
        Gets the config file field specifying the width of the ortholog histogram bins.

        :return bin_width_ortho: number (float or integer) for bin width in ortholog Ks histogram
        """
        bin_width_ortho = float(self.config.get("PARAMETERS", "bin_width_orthologs"))
        return bin_width_ortho

    def get_x_lim_ortho(self):
        """
        Gets the config file field specifying the upper limit of the x axis range in the ortholog Ks distribution plots.

        :return x_lim_ortho: integer or float
        """
        x_lim_ortho = float(self.config.get("PARAMETERS", "x_axis_max_limit_orthologs_plots", fallback="5"))
        return x_lim_ortho

    def get_x_max_lim(self):
        """
        Gets the config file field specifying the upper limit of the x axis range in the paralog/mixed distribution plot.

        :return x_max_lim: integer or float
        """
        x_max_lim = float(self.config.get("PARAMETERS", "x_axis_max_limit_paralogs_plot", fallback="5"))
        return x_max_lim

    def get_y_lim(self):
        """
        Gets the config file field specifying the upper limit of the y axis range in the paralog/mixed distribution plot.

        :return: the upper limit as a floating number or as None string
        """
        y_lim = self.config.get("PARAMETERS", "y_axis_max_limit_paralogs_plot")  # by default it's "None"
        try:
            y_lim = float(y_lim) # returning the y_lim as a float
        except Exception:
            y_lim = None # returning the y_lim as a string (either the default "None" or empty string or else)
        return y_lim

    def get_color_list(self):
        """
        Gets the config file field of the color list for the divergence lines.
        All the divergence lines coming from the same divergence node in the tree share the
        same color. The first color of the list is assigned to the first internal node 
        encountered when going from the focal species leaf up to the root. The second color 
        is assigned to the second internal node encountered along this path, and so on.
        There must be at least as many colors as the number of divergence nodes.

        Checks if there are colors whose name is not recognized by matplotlib, e.g. misspelled.

        :return colors: list of colors
        """
        color_list_string = self.config.get("PARAMETERS", "divergence_colors")
        colors = [c.strip() for c in color_list_string.split(',')]
        if len(colors) == 1 and colors[0] == "":
            logging.error('Field "divergence_colors" in configuration file is empty, please fill in')
            logging.error("Exiting.")
            sys.exit(1)

        # Check if color names are recognized by matplotlib
        faulty_color_names = []
        for color in colors:
            if not is_color_like(color):
                faulty_color_names.append(color)
        if len(faulty_color_names) != 0:
            logging.error('Field "divergence_colors" in configuration file contains color names not recognized by Matplotlib, please adjust the following:')
            for color in faulty_color_names:
                logging.error(f"- {color}")
            sys.exit(1)
        return colors


    def get_logging_level(self):
        """
        Checks the logging message level provided in the expert configuration file.
        If the level is not among the available ones, prints an message and sets INFO as default level.
        Available level options: CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET.

        :return level: logging message level that will be used in the pipeline
        """
        if self.expert_config is not None:
            try:
                level = self.expert_config.get("EXPERT PARAMETERS", "logging_level").upper()
                if level not in ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG", "NOTSET"]:
                    logging.basicConfig(format='%(levelname)s\t%(message)s', level="INFO", stream=sys.stdout)
                    logging.warning(f'Unrecognized logging level in configuration file ["{level}"]; please choose among "CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET. Default logging level "INFO" will be applied.')
                    level = "INFO"
            except Exception:
                logging.basicConfig(format='%(levelname)s\t%(message)s', level="INFO", stream=sys.stdout)
                logging.warning(f'Missing field in expert configuration file [logging level]; please choose among "CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET. Default logging level "INFO" will be applied')
                level = "INFO"
        else:
                level = "INFO"
        return level


    def get_preserve_ks_tmp_files(self):
        """
        Checks whether "preserve" is required or not in the expert configuration file.
        If set to "yes", temporary files produced during Ks estimate won't be removed.
        Default is "no", i.e. removing.

        :return preserve: boolean for preserving (True) or not (False) Ks temporary files
        """
        if self.expert_config is not None:
            try:
                preserve = self.expert_config.get("EXPERT PARAMETERS", "preserve_ks_tmp_files").lower()
                if preserve not in ["yes", "no"]:
                    logging.warning(f'Unrecognized field in expert configuration file [preserve = {preserve}]. Please choose between "yes" and "no". Default choice will be applied [no]')
                    preserve = False
                else:
                    if preserve == "yes":
                        preserve = True
                    elif preserve == "no":
                        preserve = False
            except Exception:
                logging.warning(f'Missing field [preserve] in expert configuration file. Please choose between "yes" and "no". Default choice will be applied [no]')
                preserve = False
        else:
            preserve = False
        return preserve

    
    def get_peak_stats(self):
        """
        Checks the statistical measure used to get a representative Ks value for an ortholog Ks distribution that
        is decided by the user in the expert configuration file.
        The "mode" option gets the mode of the distribution as distribution peak, while the "median" option will take its median. 
        If the field choice is unrecognized, a message is returned and mode is used [Default: mode].

        :return peak_stats: defines the statistics measure that will be used to get the peaks (mode or median)
        """
        if self.expert_config is not None:
            try:
                peak_stats = self.expert_config.get("EXPERT PARAMETERS", "distribution_peak_estimate").lower()
                if peak_stats != "mode" and peak_stats != "median":
                    logging.warning(f'Unrecognized field in expert configuration file [distribution_peak_estimate = {peak_stats}]. Please choose between "mode" and "median". Default choice will be applied [mode]')
                    peak_stats = "mode"
            except Exception:
                peak_stats = "mode"     
        else:
            peak_stats = "mode"
        return peak_stats


    def plot_correction_arrows(self):
        """
        Checks whether the user wants to show the shifts of the adjusted lines as an arrow
        starting from the original position and ending on the adjusted position.

        :return correction_arrows: flag that states whether to plot the arrows or not (True/False)
        """
        if self.expert_config is not None:
            try:
                correction_arrows = self.expert_config.get("EXPERT PARAMETERS", "plot_adjustment_arrows").lower()
                if correction_arrows not in ["yes", "no"]:
                    logging.warning(f'Unrecognized field in expert configuration file [plot_adjustment_arrows = {correction_arrows}]. Please choose between "yes" and "no". Default choice will be applied [no]')
                    correction_arrows = False
                else:
                    if correction_arrows == "yes":
                        correction_arrows = True
                    elif correction_arrows == "no":
                        correction_arrows = False
            except Exception:
                logging.warning(f'Missing field in expert configuration file [plot_adjustment_arrows]. Please choose between "yes" and "no". Default choice will be applied [no]')
                correction_arrows = False
        else:
            correction_arrows = False
        return correction_arrows


    def get_kde_bandwidth_modifier(self):
        """
        The default KDE bandwidth is computed by applying Scott's rule, but the
        resulting KDE is usually too smooth on the actual distribution. That's why
        the code by default reduces the kde.factor multiplying it by 0.4 (a "modifier").

        However, for very steep WGD peaks the user might need to reduce the factor even more and
        can edit here the modifier (e.g. set it to 0.2). Note that a tighter distribution will 
        likely catch distribution noise. 
        To make instead the KDE smoother multiply by a modifier greater than 0.4;
        to set the original factor set the modifier to 1.
        Allowed floats and integer values  [Default: 0.4].

        :return modifier: number between 0 excluded and 1 included.
        """
        if self.expert_config is not None:
            try:
                modifier = self.expert_config.get("EXPERT PARAMETERS", "kde_bandwidth_modifier")
                try:
                    modifier = literal_eval(modifier)
                except Exception:
                    pass
                if (not isinstance(modifier, int) and not isinstance(modifier, float)) or modifier == 0:
                    logging.warning(f"Unrecognized field in expert configuration file [kde_bandwidth_modifier = {modifier}]. Please enter a non-zero number (0 < n <= 1). Default choice will be applied [0.4]")
                    modifier = 0.4
            except Exception:
                logging.warning(f"Missing field in expert configuration file [kde_bandwidth_modifier]. Please enter a non-zero number (0 < n <= 1). Default choice will be applied [0.4]")
                modifier = 0.4
        else:
            modifier = 0.4
        return modifier


    def get_max_EM_iterations(self):
        """
        Gets the maximum number of EM iterations to apply during mixture modeling [Default: 600].

        :return max_mixture_model_iterations: max number of iterations for the exponential-maximization
        algorithm 
        """
        if self.expert_config is not None:
            try:
                max_iter = self.expert_config.get("EXPERT PARAMETERS", "max_mixture_model_iterations")
                try:
                    max_iter = literal_eval(max_iter)
                except Exception:
                    pass
                if not isinstance(max_iter, int):
                    logging.warning(f'Unrecognized field in expert configuration file [max_mixture_model_iterations = {max_iter}]. Please choose a positive integer. Default choice will be applied [600]')
                    max_iter = 600
                elif max_iter <= 0:
                    logging.warning(f'Unrecognized field in expert configuration file [max_mixture_model_iterations = {max_iter}]. Please choose a positive integer. Default choice will be applied [600]')
                    max_iter = 600
                elif max_iter <= 300:
                    logging.warning(f"A small maximum number of mixture model iterations [max_mixture_model_iterations = {max_iter}] can produce poor fitting.")
                elif max_iter > 600:
                    logging.warning(f"A large maximum number of mixture model iterations [max_mixture_model_iterations = {max_iter}] can increase the runtime.")
            except Exception:
                logging.warning(f'Missing field in expert configuration file [max_mixture_model_iterations]. Please choose a positive integer. Default choice will be applied [600]')
                max_iter = 600
        else:
            max_iter = 600
        return max_iter


    def get_num_EM_initializations(self):
        """
        Gets the number of EM initialization iterations to be applied to
        the k-means during lognormal mixture modeling and to the 
        random initialization during the exponential-lognormal mixture modeling.
        [Default: 10].

        :return num_mixture_model_initializations: number of times that the expectation-maximization
        algorithm is initialized
        """
        if self.expert_config is not None:
            try:
                n_init = self.expert_config.get("EXPERT PARAMETERS", "num_mixture_model_initializations")
                try:
                    n_init = literal_eval(n_init)
                except Exception:
                    pass
                if not isinstance(n_init, int):
                    logging.warning(f'Unrecognized field in expert configuration file [num_mixture_model_initializations = {n_init}]. Please choose a positive integer. Default choice will be applied [10]')
                    n_init = 10
                elif  n_init <= 0:
                    logging.warning(f'Unrecognized field in expert configuration file [num_mixture_model_initializations = {n_init}]. Please choose a positive integer. Default choice will be applied [10]')  
                    n_init = 10    
                elif n_init < 5:
                    logging.warning(f"A small number of mixture model initializations [num_mixture_model_initializations = {n_init}] can produce poor fitting.")
                elif n_init > 10:
                    logging.warning(f"A large number of mixture model initializations [num_mixture_model_initializations = {n_init}] can increase the runtime.")
            except Exception:
                logging.warning(f'Missing field in expert configuration file [num_mixture_model_initializations]. Please choose a positive integer. Default choice will be applied [10]')
                n_init = 10
        else:
            n_init = 10
        return n_init


    def get_extra_paralogs_analyses_methods(self):
        """
        A flag to state whether to perform the additional peak calling methods beside the
        default ones.

        :return extra_paralogs_analyses_methods: boolean
        """
        if self.expert_config is not None:
            try:
                extra_paralogs_analyses_methods = self.expert_config.get("EXPERT PARAMETERS", "extra_paralogs_analyses_methods").lower()
                if extra_paralogs_analyses_methods not in ["yes", "no"]:
                    logging.warning(f'Unrecognized field in expert configuration file [extra_paralogs_analyses_methods = {extra_paralogs_analyses_methods}]. Please choose between "yes" and "no". Default choice will be applied [no]')
                    extra_paralogs_analyses_methods = False
                else:
                    if extra_paralogs_analyses_methods == "yes":
                        extra_paralogs_analyses_methods = True
                    elif extra_paralogs_analyses_methods == "no":
                        extra_paralogs_analyses_methods = False
            except Exception:
                logging.warning(f'Missing field in expert configuration file [extra_paralogs_analyses_methods]. Please choose between "yes" and "no". Default choice will be applied [no]')
                extra_paralogs_analyses_methods = False
        else:
            extra_paralogs_analyses_methods = False
        return extra_paralogs_analyses_methods


    def get_max_mixture_model_components(self):
        """
        Gets the maximum number of components used in the mixture models (i.e. the 
        exp-lognormal mixture model with random components and the lognormal mixture model).
        Minimum required is 3, which includes the exponential component, one buffer lognormal
        component. 
        [Default: 5]. A higher number of components may be useful when the paralog
        distribution is believed to have undergone many WGDs (suggested: more than 3), but 
        comes along with increased chance of overfitting and thus over-interpretation of WGM signals.

        :return max_mixture_model_components: number of times that the expectation-maximization
        algorithm is initialized
        """
        if self.expert_config is not None:
            try:
                max_comp = self.expert_config.get("EXPERT PARAMETERS", "max_mixture_model_components")
                try:
                    max_comp = literal_eval(max_comp)
                except Exception:
                    pass
                if not isinstance(max_comp, int):
                    logging.warning(f'Unrecognized field in expert configuration file [max_mixture_model_components = {max_comp}]. Please choose a positive integer >= 2. Default choice will be applied [5]')
                    max_comp = 5
                elif  max_comp <= 0:
                    logging.warning(f'Unrecognized field in expert configuration file [max_mixture_model_components = {max_comp}]. Please choose a positive integer >= 2. Default choice will be applied [5]')  
                    max_comp = 5
                elif max_comp == 1:
                    logging.warning(f'Field "max_mixture_model_components" has been changed from {max_comp} to the minimum required, 2.')
                    max_comp = 2 # exponential + buffer
                elif max_comp <= 3:
                    logging.warning(f"A low number of mixture model components [max_mixture_model_components = {max_comp}] can produce poor fitting.")
                elif max_comp >= 7:
                    logging.warning(f"A high number of mixture model components [max_mixture_model_components = {max_comp}] increases overfitting risk.")
            except Exception:
                logging.warning(f'Missing field in expert configuration file [max_mixture_model_components]. Please choose a positive integer. Default choice will be applied [5]')
                max_comp = 5
        else:
            max_comp = 5
        return max_comp


    def get_max_ks_for_mixture_model(self, max_ks_para):
        """
        Gets the upper limit for the Ks range in which the mixture model algorithm will
        perform the fitting. This upper Ks value should be placed around the visible
        Ks coordinate where the paralog distribution starts showing only a flat right tail
        without any WGM trace. Species with low substitution rates are likely to have 
        visible signals only up around 3 Ks, thus this parameter makes sure that the mixture model
        is performed only on the relative range 0-to-3 Ks. Species with high rates are instead
        likely to have some signal up to a higher Ks, such as 5 Ks. It is not advised to
        set this parameter at more than 5 Ks since high Ks are not reliable.
        [Default: 5].

        :param max_ks_para: maximum paralog Ks accepted in the analysis
        :return max_ks_EM: number of times that the expectation-maximization
        algorithm is initialized
        """
        if self.expert_config is not None:
            try:
                max_ks_EM = self.expert_config.get("EXPERT PARAMETERS", "max_mixture_model_ks")
                try:
                    max_ks_EM = literal_eval(max_ks_EM)
                except Exception:
                    pass
                if (not isinstance(max_ks_EM, int) and not isinstance(max_ks_EM, float)) or max_ks_EM <= 0:
                    logging.warning(f'Unrecognized field in expert configuration file [max_mixture_model_ks = {max_ks_EM}]. Please enter a positive integer or float. Default choice will be applied [5]')
                    max_ks_EM = 5
            except Exception:
                logging.warning(f'Missing field in expert configuration file [max_mixture_model_ks]. Please enter a positive integer or float. Default choice will be applied [5]')
                max_ks_EM = 5
        else:
            max_ks_EM = 5
        return max_ks_EM


    def get_min_ks_anchors(self):
        """
        Gets in the expert configuration file the lower limit for the Ks range used to build the anchor pair Ks distribution.
        Removing very small Ks values (<0.05) in the anchor distribition can help removing noise at recent Ks age that
        could be mistakenly interpreted as a very recent WGM. However, some WGMs with actual very young age might thereby be
        erroneously cut off. This parameter allows user customization. [Default: 0.05].

        :param min_ks_anchors: minimum anchor pair Ks value used to build the anchor pair distribution
        """
        if self.expert_config is not None:
            try:
                min_ks_anchors = self.expert_config.get("EXPERT PARAMETERS", "min_ks_anchor_pairs")
                try:
                    min_ks_anchors = literal_eval(min_ks_anchors)
                except Exception:
                    pass
                if (not isinstance(min_ks_anchors, int) and not isinstance(min_ks_anchors, float)) or min_ks_anchors < 0:
                    logging.warning(f'Unrecognized field in expert configuration file [min_ks_anchor_pairs = {min_ks_anchors}]. Please enter a non-negative integer or float.')
                    exit()
            except Exception:
                logging.warning(f'Missing field in expert configuration file [min_ks_anchor_pairs]. Please enter a positive integer or float. Default choice will be applied [0.05]')
                min_ks_anchors = 0.05
        else:
            min_ks_anchors = 0.05
        return min_ks_anchors
    

    def get_max_gene_family_size(self):
        """
        Gets the maximum size for a gene family in order to be considered and analyzed by wgd. 

        :return max_gene_family_size: maximum number of genes in a gene family accepted for wgd analysis
        """
        if self.expert_config is not None:
            try:
                max_size = self.expert_config.get("EXPERT PARAMETERS", "max_gene_family_size")
                try:
                    max_size = literal_eval(max_size)
                except Exception:
                    pass
                if not isinstance(max_size, int):
                    logging.warning(f'Unrecognized field in expert configuration file [max_gene_family_size = {max_size}]. Please choose a positive integer. Default choice will be applied [200]')
                    max_size = 200
                elif max_size <= 0:
                    logging.warning(f'Unrecognized field in expert configuration file [max_gene_family_size = {max_size}]. Please choose a positive integer. Default choice will be applied [200]')
                    max_size = 200
            except Exception:
                logging.warning(f'Missing field in expert configuration file [max_gene_family_size]. Please choose a positive integer. Default choice will be applied [200]')
                max_size = 200
        else:
            max_size = 200
        return max_size


    def get_reciprocal_retention_top(self, reciprocal_retention):
        """
        Gets the number of top reciprocally retained gene families to be considered out of the total ranked 9178 ones.

        :return top: integer or float
        """
        if reciprocal_retention:
            if self.expert_config is not None:
                # Get user-defined top value in field "top_reciprocally_retained_gfs"
                # NOTE: If FIELD "top_reciprocally_retained_gfs" is missing in config file, "top" variable falls back to 2000
                # NOTE: If instead only the related VALUE is not present, "top" variable is an empty string
                top = self.expert_config.get("EXPERT PARAMETERS", "top_reciprocally_retained_gfs", fallback="2000")
                # Convert to integer the user-defined value or the fallback value
                try:
                    top = int(top)
                # If the VALUE was left empty in by user ("top" variable is empty string), assume again 2000
                except ValueError:
                    logging.warning("Field [top_reciprocally_retained_gfs] in expert config file was left empty: assuming top 2000 GFs")
                    top = 2000
            else:
                top = 2000
        else:
            top = None
        return top


    def get_reciprocal_retention_rank_type(self, reciprocal_retention):
        """
        Returns lambda ranking or None if reciprocal retention analysis not requested.
        
        Possible integration in the future:
        Gets the rank type ("lambda" or "combined") for the reciprocally retained gene family list. The lambda rank is based 
        on the lambda parameter (the lower this is, the more the gene family is duplicated mainly through WGM with no gene loss),
        while the combined rank takes also into account the block duplicate percentage as 
        an independent source of evidence for reciprocal retention.

        :return rank_type: string "lambda" or None.
        """
        if reciprocal_retention:
            return "lambda"
        else:
            return None
        # The next block is ignored, as only the lambda ranking is currently supported
        if self.expert_config is not None:
            rank_type = self.expert_config.get("EXPERT PARAMETERS", "reciprocal_retention_rank_type", fallback="lambda").lower()
            if rank_type not in ["lambda", "combined"]:
                logging.warning(f'Unrecognized field in expert configuration file [reciprocal_retention_rank_type = {rank_type}]. \
                                  Please choose between "lambda" and "combined". Default choice will be applied ["lambda"]')
                rank_type = "lambda"
        else:
            rank_type = "lambda"
        return rank_type


    def get_orthomcl_inflation(self, reciprocal_retention):
        """
        Gets the inflation parameter for OrthoMCL-based orthogroup clustering.
        It affects cluster granularity (tightness):
        - low values  --> few large clusters
        - high values --> many small clusters
        Default and recommended is 1.5.
        
        :return top: inflation, as float
        """
        if reciprocal_retention:
            if self.expert_config is not None:
                inflation = float(self.expert_config.get("EXPERT PARAMETERS", "orthomcl_inflation", fallback="1.5"))
            else:
                inflation = 1.5
        else:
            inflation = None
        return inflation
    
    def get_orthomcl_version(self):
        """
        A flag to state whether to run the original OrthoMCL v1.4 or the edited version
        called OrthoMCLight that significantly speeds up the runtime and reduces memory usage.
        By default, the edited version is used.

        :return use_original_orthomcl_version: boolean
        """
        if self.expert_config is not None:
            try:
                use_original_orthomcl_version = self.expert_config.get("EXPERT PARAMETERS", "use_original_orthomcl_version").lower()
                if use_original_orthomcl_version not in ["yes", "no"]:
                    logging.warning(f'Unrecognized field in expert configuration file [use_original_orthomcl_version = {use_original_orthomcl_version}]. Please choose between "yes" and "no". Default choice will be applied [no]')
                    use_original_orthomcl_version = False
                else:
                    if use_original_orthomcl_version == "yes":
                        use_original_orthomcl_version = True
                    elif use_original_orthomcl_version == "no":
                        use_original_orthomcl_version = False
            except Exception:
                logging.warning(f'Missing field in expert configuration file [use_original_orthomcl_version]. Please choose between "yes" and "no". Default choice will be applied [no]')
                use_original_orthomcl_version = False
        else:
            use_original_orthomcl_version = False
        return use_original_orthomcl_version


    def get_max_extra_original_genes_in_new_gfs(self, reciprocal_retention):
        """
        In the context of reconstructing the original rec.ret. GFs for the focal species,
        the genes of an original top GF ("original genes") are compared to the 
        genes of a newly obtained OrthoMCL GF ("new GF").
        Some of the original genes are shared with the new GF, while other genes are not ("extra"). 
        This function gets the maximum allowed number of extra original genes
        ("max_extra_original_genes_in_new_gfs") [Default: 15].
        
        E.g. There are 50 genes in the original GF, and 45 of them are in common with 
        a certain new GF. There are therefore only 5 original genes not shared (extra)
        with the new GF ("extra"). This is an acceptable number, according to 
        the default threshold (15).
        
        :return top: inflation, as float
        """
        if reciprocal_retention:
            if self.expert_config is not None:
                max_extra_original_genes_in_new_gfs = int(self.expert_config.get("EXPERT PARAMETERS", "max_extra_original_genes_in_new_gfs", fallback="15"))
            else:
                max_extra_original_genes_in_new_gfs = 15
        else:
            max_extra_original_genes_in_new_gfs = None
        return max_extra_original_genes_in_new_gfs


    def min_common_old_genes_in_new_gfs(self):
        """
        This is current not effectively used.        
        :return top: inflation, as float
        """
        if self.expert_config is not None:
            min_common_old_genes_in_new_gfs = int(self.expert_config.get("EXPERT PARAMETERS", "min_common_old_genes_in_new_gfs", fallback="3"))
        else:
            min_common_old_genes_in_new_gfs = 3
        return min_common_old_genes_in_new_gfs


    # NOTE: The following getters are unused (diamond-based reciprocally retained pipeline not implemented)
    # TODO: To be removed
    
    # def get_reciprocal_retention_pident_threshold(self):
    #     """
    #     Gets the minimum identity percentage required in the local alignment between a main species sequence and 
    #     the top reciprocally retained sequences from the top orthogroups. This number is given to diamond "blastx"
    #     command as the "-id" argument. Diamond hits with less than e.g. 40 won't be part of the diamond output table.

    #     :return pident_threshold: integer or float. Default is 40
    #     """
    #     if self.expert_config is not None:
    #         pident_threshold = int(self.expert_config.get("EXPERT PARAMETERS", "pident_threshold", fallback="40"))
    #     else:
    #         pident_threshold = 40
    #     return pident_threshold

    # def get_reciprocal_retention_min_qcov_threshold(self):
    #     """
    #     Gets the minimum coverage percentage in the main species sequence that is required in the alignment between 
    #     the main species sequence and the top reciprocally retained sequences from the top orthogroups. Diamond hits
    #     with less than e.g. 50 will be removed from the diamond output table, so to remove the hits with bad coverage.

    #     :return top: integer or float
    #     """
    #     if self.expert_config is not None:
    #         min_qcov_threshold = int(self.expert_config.get("EXPERT PARAMETERS", "min_qcov_threshold", fallback="50"))
    #     else:
    #         min_qcov_threshold = 50
    #     return min_qcov_threshold

    # def get_filter_by_false_predictions(self):
    #     """
    #     Gets whether to use all the top gene families (set to "no") or a subset of gene families (set to "yes") which for the three test species
    #     (A. thaliana, P. dactylifera and V. vinifera) had few false predictions in the confusion matrix (the total amount of false positives 
    #     and negatives was less than half of the total amount of true positives).

    #     :return boolean: boolean
    #     """
    #     if self.expert_config is not None:
    #         filter_false_predictions = self.expert_config.get("EXPERT PARAMETERS", "filter_by_false_predictions", fallback="no").lower()

    #         if filter_false_predictions == "yes":
    #             return True
    #         elif filter_false_predictions == "no":
    #             return False
    #         else:
    #             logging.warning(f'Unrecognized field in expert configuration file [filter_by_false_predictions = {filter_false_predictions}]. \
    #                                 Please choose between "yes" and "no". Default choice will be applied ["no"]')
    #             return False
    #     else:
    #         return False
        

    # def get_filter_by_zero_false_predictions(self):
    #     """
    #     Gets whether to use all the top gene families (set to "no") or a subset of gene families (set to "yes") which for the three test species
    #     (A. thaliana, P. dactylifera and V. vinifera) had no false positives or negatives.

    #     :return boolean: boolean
    #     """
    #     if self.expert_config is not None:
    #         filter_zero_false_predictions = self.expert_config.get("EXPERT PARAMETERS", "filter_by_zero_false_predictions", fallback="no").lower()

    #         if filter_zero_false_predictions == "yes":
    #             return True
    #         elif filter_zero_false_predictions == "no":
    #             return False
    #         else:
    #             logging.warning(f'Unrecognized field in expert configuration file [filter_by_zero_false_predictions = {filter_zero_false_predictions}]. \
    #                                 Please choose between "yes" and "no". Default choice will be applied ["no"]')
    #             return False
    #     else:
    #         return False
