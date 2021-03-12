import sys
import os
import pandas
import logging
from ast import literal_eval
from pandas import DataFrame
import ksrates.fc_configfile as fcConf
import ksrates.fc_check_input as fcCheck
import ksrates.fc_manipulate_trees as fcTree
from ksrates.utils import init_logging


def setup_correction(config_file, nextflow_flag):
    config = fcConf.Configuration(config_file)

    init_logging("Setting up the analysis from configuration file", config.get_logging_level())
    logging.info("Loading parameters and input files")

    # Check configfile
    species_of_interest = config.get_species()
    original_tree = config.get_newick_tree()
    tree = fcTree.reorder_tree_leaves(original_tree, species_of_interest)  # focal species is the top leaf
    latin_names = config.get_latin_names()

    trigger_exit = False # will trigger exit if FASTA or GFF files are missing or empty
    fasta_dict = config.get_fasta_dict()
    gff_dict = config.get_gff_dict(warn_empty_dict=False)
    # If a GFF is provided, check existence and content
    if species_of_interest in gff_dict:
        gff = config.get_gff_name(gff_dict, species_of_interest)
        if fcCheck.check_file_nonexistent_or_empty(gff, "GFF file") == True:
            trigger_exit = True

    db_path = config.get_ortho_db()
    ks_list_db_path = config.get_ks_db()
    max_num_outspecies = config.get_max_num_outspecies()
    try:
        max_num_outspecies = int(max_num_outspecies)
    except Exception:
        pass

    # Checking if the IDs in FASTA (and GFF) files are likely to raise an error when using the PAML package.
    # For example, paml 4.4 gives "Node" Keyerror if the len(IDs) > 50, and/or there are special characters,
    # and/or there are two spaces in a row.
    logging.info(f"Checking if sequence data files exist and if sequence IDs are compatible with wgd pipeline...")
    all_species = []
    for leaf in tree.get_leaves():
        all_species.append(leaf.name)

    for species in all_species:
        # Check existence and content of FASTA file
        fasta = config.get_fasta_name(fasta_dict, species)
        if fcCheck.check_file_nonexistent_or_empty(fasta, "FASTA file") == True:
            trigger_exit = True
            continue
        # Check the IDs in FASTA file
        if species == species_of_interest and species in gff_dict:  # Warn about both FASTA and GFF files
            fcCheck.check_IDs(fasta, latin_names[species], gff)
        else:  # Warn only about FASTA file
            fcCheck.check_IDs(fasta, latin_names[species])
            
    if trigger_exit:
        logging.error("Please add the missing information in [fasta_filenames] configuration file field and rerun the analysis. Exiting.")
        exit(1)
    logging.info("Completed")

    # Creating folders for correction output files
    if not os.path.exists('rate_adjustment'):
        os.mkdir('rate_adjustment')
    if not os.path.exists(os.path.join("rate_adjustment", f"{species_of_interest}")):
        os.mkdir(os.path.join("rate_adjustment", f"{species_of_interest}"))
        logging.info(f"Creating output folder [rate_adjustment/{species_of_interest}]")

    # -----------------------------------------------------------------------------

    # 1) FINDING SISTER AND OUTGROUP SPECIES FOR EACH NODE OF INTEREST

    logging.info("")
    logging.info(f"Extracting ortholog trios for rate-adjustment [ortholog_trios_{species_of_interest}.tsv]")

    if isinstance(max_num_outspecies, int):
        logging.info(f"- Each divergent species pair will be rate-adjusted with at maximum {max_num_outspecies} "
                     f"trios by using the closest outspecies")
        logging.info(f"  (as required in configuration file field 'maximum_number_outspecies')")
    else:
        logging.info(f"- Each divergent species pair will be rate-adjusted by using all the possible outspecies "
                     f"found in the tree.")

    # get tree node object of the focal species
    species_of_interest_node = fcTree.get_species_node(species_of_interest, tree)
    # get the list of ancestors (as tree node objects) in the lineage that lead to the focal species
    sp_history = fcTree.get_species_history(species_of_interest_node)
    # Checking if the focal species has at least one outgroup in the provided tree
    if len(sp_history)-2 == 0:
        logging.error("")
        logging.error(f"Species [{species_of_interest}] has no outgroup in the provided Newick tree "
                      f"and the rate-adjustment can't be performed.")
        logging.error(f"Please add at least one outgroup species or change the focal species.")
        sys.exit(1)

    trios_array = []  # list of trios
    outfile_drawing_path = os.path.join("rate_adjustment", f"{species_of_interest}",
                                        f"tree_{species_of_interest}.txt")
    with open(outfile_drawing_path, "w+") as outfile_drawing:
        outfile_drawing.write(f"Focal species: {species_of_interest}\n\n")

        # Obtaining the numeric labels for internal nodes relevant in the species analysis
        fcTree.labeling_internal_nodes(species_of_interest_node)
        
        node = 0
        while node < len(sp_history)-2:
            # the name label to be shown in the ASCII tree will start from 1 and not from 0
            outfile_drawing.write(f"Node {node+1}:\n")
            currentnode = sp_history[node]

            # GETTING SISTER SPECIES
            sisters = fcTree.get_sister_species_of_a_node(currentnode)
            outfile_drawing.write(f"Sister species:      {', '.join(sisters)}\n")

            # GETTING OUTSPECIES
            outspecies = fcTree.get_outspecies_of_a_node(currentnode, max_num_outspecies)
            outfile_drawing.write(f"Outgroup species:    {', '.join(outspecies)}\n\n")

            # APPENDING TRIOS (a trio is composed of focal species, sister species and outgroup species)
            for s in sisters:
                for o in outspecies:
                    trios_array.append([node+1, species_of_interest, s, o])
            node += 1

        print_tree = tree.get_ascii(attributes=["name"], show_internal=True)
        outfile_drawing.write(f"{print_tree}\n\n")

    logging.info(f"- Total number of trios: {len(trios_array)}")
    logging.info("")

    # Generate trios DataFrame from trios array
    trios_df = DataFrame.from_records(trios_array, columns=["Node", "Species", "Sister_Species", "Out_Species"])
    outfile_trios_path = os.path.join("rate_adjustment", f"{species_of_interest}",
                                         f"ortholog_trios_{species_of_interest}.tsv")
    with open(outfile_trios_path, "w+") as outfile:
        outfile.write(trios_df.to_csv(sep="\t", index=False))

    # -----------------------------------------------------------------------------

    # 2) FINDING UNKNOWN PAIRS FOR wgd ORTHOLOG RUNS
    logging.info(f"Extracting species pairs for wgd ortholog pipeline [ortholog_pairs_{species_of_interest}.tsv]")
    if isinstance(max_num_outspecies, int):
        logging.info(f"- Only species pairs required by the selected trios will be considered.")

    # THIS CODE JUST TAKES THE PAIRS THAT ARE NECESSARY FOR THE TRIOS
    species_pairs = []
    for trio in trios_array:
        combos = [[trio[1], trio[2]], [trio[1], trio[3]], [trio[2], trio[3]]]
        for pair in combos:
            if pair not in species_pairs:
                species_pairs.append(pair)

    species_pairs = [sorted(x, key=str.casefold) for x in species_pairs]
    species_pairs_unknown = []

    # Getting possible extra species pairs whose ortholog data needs to be computed in order to
    # be able to fill in all branch lengths in the tree figure
    missing_pairs_with_latin_names, missing_pairs = fcTree.find_missing_pairs_for_tree_rates(tree,
                                                                                             species_of_interest_node,
                                                                                             sp_history, latin_names)
    # Adding these species pairs to the list obtained from the trios
    species_pairs.extend(missing_pairs)

    tags_list = [] 
    missing_latin_names = []
    flag_exit = False
    for pair in species_pairs:
        sp1, sp2 = pair[0], pair[1]
        try:
            latinSp1 = latin_names[sp1]
        except KeyError:
            if sp1 not in missing_latin_names:
                missing_latin_names.append(sp1)
                logging.error(f"Latin name for [{sp1}] not found in the configuration file")
                flag_exit = True
        try:
            latinSp2 = latin_names[sp2]
        except KeyError:
            if sp2 not in missing_latin_names:
                missing_latin_names.append(sp2)
                logging.error(f"Latin name for [{sp2}] not found in the configuration file")
                flag_exit = True
        try:
            latin_tag = "_".join(sorted([latinSp1, latinSp2], key=str.casefold))
            tags_list.append([latin_tag, [sp1, sp2]])
        except Exception:
            pass
    if flag_exit:  # exits if one or more scientific names are missing
        logging.error("Exiting")
        sys.exit(1)

    try:
        with open(db_path, "r") as f:
            db = pandas.read_csv(f, sep="\t", index_col=0)
            for tags in tags_list:
                if tags[0] not in db.index:
                    species_pairs_unknown.append([tags[1][0], tags[1][1]])
            no_db_file = False
    except Exception:
        no_db_file = True

    try:
        with open(ks_list_db_path, "r") as f:
            ks_list_db = pandas.read_csv(f, sep="\t", index_col=0)
            # When imported from csv format, all the Ks lists in the df are read as
            # plain text (strings) and must be converted back to value lists
            ks_list_db.loc[:, 'Ks_Values'] = ks_list_db.loc[:, 'Ks_Values'].apply(literal_eval)
            for tags in tags_list:
                if tags[0] not in ks_list_db.index:
                    if [tags[1][0], tags[1][1]] not in species_pairs_unknown:
                        species_pairs_unknown.append([tags[1][0], tags[1][1]])
        no_ks_list_db_file = False
    except Exception:
        no_ks_list_db_file = True

    if no_db_file and not no_ks_list_db_file:
        logging.info("Ortholog peak database empty or not found at the path provided in the configuration file.")
        logging.info("All of the species pairs will be used to compute their ortholog distribution and its peak.")
        species_pairs_unknown = species_pairs.copy()
    elif not no_db_file and no_ks_list_db_file:
        logging.info("Ortholog Ks list database empty or not found at the path provided in the configuration file.")
        logging.info("All of the species pairs will be used to compute their ortholog distribution and its peak.")
        species_pairs_unknown = species_pairs.copy()
    elif no_db_file or no_ks_list_db_file:
        logging.info("Ortholog databases empty or not found at the paths provided in the configuration file.")
        logging.info("All of the species pairs will be used to compute their ortholog distribution and its peak.")
        species_pairs_unknown = species_pairs.copy()
    else:  # False and False
        pass

    logging.info(f"- Total number of species pairs not in database(s): {len(species_pairs_unknown)}")
    logging.info("")

    species_pairs_unknown_df = DataFrame.from_records(species_pairs_unknown, columns=["Species1", "Species2"])
    outfile_pairs_path = os.path.join("rate_adjustment", f"{species_of_interest}",
                                      f"ortholog_pairs_{species_of_interest}.tsv")
    with open(outfile_pairs_path, "w+") as outfile_pairs:
        outfile_pairs.write(species_pairs_unknown_df.to_csv(sep="\t", index=False))

    # -----------------------------------------------------------------------------

    # 3) PLOTTING THE ORIGINAL UN-CORRECTED TREE in PDF FORMAT
    logging.info(f"Plotting input phylogenetic tree [{fcTree._TREE.format(species_of_interest)}]")
    logging.info("")
    fcTree.plot_uncorrected_phylogeny(tree, species_of_interest, latin_names, sp_history)

    # -----------------------------------------------------------------------------

    # 4) IF OUTSIDE NEXTFLOW PIPELINE:
    # Generating an output file containing the list of wgd runs that needs to be done by wgd_orthologs.py
    # By default it is assumed no Nextflow and the list of commands will be created
    if not nextflow_flag:
        logging.info(f"Generating the list of wgd paralog and ortholog runs to be manually executed "
                     f"[wgd_runs_{species_of_interest}.txt]")
        logging.info("")
        with open(os.path.join(f"wgd_runs_{species_of_interest}.txt"), "w+") as wgd_runs:
            default_threads = 1
            wgd_runs.write(f"# Note: the number of threads can be increased depending on the number of cores [default: {default_threads}]\n\n")
            wgd_runs.write(f"ksrates paralogs-ks config_{species_of_interest}.txt --n-threads={default_threads}\n")
            for pair in species_pairs_unknown:
                sp1, sp2 = pair[0], pair[1]
                wgd_runs.write(f"ksrates orthologs-ks config_{species_of_interest}.txt {sp1} {sp2} "
                               f"--n-threads={default_threads}\n")

    logging.info("All done")
