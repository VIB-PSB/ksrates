import sys
import os
import numpy as np
import logging

def check_IDs(fasta, species, gff=None):
    """
    Checks the IDs of the FASTA file of a single species. IDs longer than 50 characters, or with special characters 
    or with two space in a row may cause error in codeml in paml v.4.4. This function warns the user of that
    and suggests how to change the IDs to make them compatible.

    :param fasta: path to the FASTA file of the species
    :param species: the name of the species to which the FASTA file belongs
    :param gff: path to the GFF file of the species
    """
    too_long_ID, two_spaces = False, False
    list_of_triggered_warnings = []
    with open(fasta, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                line = line.rstrip()
                if len(line) > 50:
                    too_long_ID = True
                if "  " in line:
                    two_spaces = True
                not_allowed_characters = ""
                for character in [",", ":", "#", "(", ")", "$", "="]:
                    if character in line:
                        not_allowed_characters += f"{character} "
                if too_long_ID or two_spaces or len(not_allowed_characters) != 0:
                    if len(list_of_triggered_warnings) == 0:
                        if gff is None:
                            logging.warning(f"{species}: sequence IDs in FASTA file [{fasta}] could raise an error due to:")
                        else:
                            logging.warning(f"{species}: sequence IDs in FASTA and possibly GFF files [{fasta}, {gff}] could raise an error due to:")
                if too_long_ID:
                    if "too_long_ID" not in list_of_triggered_warnings:
                        list_of_triggered_warnings.append("too_long_ID")
                        logging.warning(" - ID length longer than 50 characters, it is advised to shorten them")
                if two_spaces:
                    if "two_spaces" not in list_of_triggered_warnings:
                        list_of_triggered_warnings.append("two_spaces")
                        logging.warning(" - ID name contains two spaces in a row: please remove at least one of them.")
                if len(not_allowed_characters) != 0:
                    if "not_allowed_characters" not in list_of_triggered_warnings:
                        list_of_triggered_warnings.append("not_allowed_characters")
                        logging.warning(f" - ID name contains one or more characters that are not allowed: {not_allowed_characters}")


def check_existence(filename, msg_prefix=""):
    """
    Gives an error message and exits when the provided file doesn't exist.
    
    :param filename: name of the file to be checked
    :param msg_prefix: (default: empty string) string to describe the file content (e.g. "Ortholog peak database")
    """
    if not os.path.exists(filename):
        # print error message separate to ensure print out
        logging.error(f"{msg_prefix} [{filename}] not found. Exiting.")
        sys.exit() # not with error code 1 (successive parts will call exit code 1 in case this file was really necessary) 


def check_inputfile(filename, msg_prefix=""):
    """
    Returns an error message and exits when the provided file doesn't exist or is empty.
    
    :param filename: name of the file to be checked
    :param msg_prefix: (default: empty string) string to describe the file content (e.g. "Ortholog peak database")
    """
    check_existence(filename, msg_prefix)
    if os.path.getsize(filename) == 0:
        # print error message separate to ensure print out
        logging.error(f"{msg_prefix} [{filename}] is empty. Exiting.")
        sys.exit() # not with error code 1 (successive parts will call exit code 1 in case this file was really necessary) 
        

def get_possible_subpaths_for_file(path):
    """
    Given a file path, it returns a list of all the possible subpaths leading to the file name.
    (Example: "path/to/file.txt" produces ["path/to/file.txt", "to/file.txt", "file.txt"]).

    :param path: the default path of an input file, where it is by default produced by the pipeline
    :return: a list with all the possible subpaths
    """
    path_list = path.split(os.sep)
    possible_subpaths_list = []
    for i in range(len(path_list)):
        possible_subpath = (os.sep).join(path_list[i:])
        possible_subpaths_list.append(possible_subpath)
    return possible_subpaths_list


def check_file_existence_and_content_in_default_paths(default_path, message_prefix):
    """
    Looks for a file in its default path by searching in all possible subpaths in order to take into account
    from which folder the user might be running the script.
    The validation step checks both the existence and that the file is not empty.

    :param default_path: the path at which the file is by default created by the pipeline
    :param message_prefix: a short definition of the file to be printed in the error message
    :return: the default (sub)path at which the file has been found or an empty string if file not found 
    """
    list_of_possible_subpaths = get_possible_subpaths_for_file(default_path)
    filename = list_of_possible_subpaths[-1]
    argument_path = "" # will store the default (sub)path in which the file will be found
    for path in list_of_possible_subpaths: 
        if os.path.exists(path):
            if message_prefix == "Ortholog pairs file": # only for this file, the existence check is enough
                argument_path = path
                return argument_path
            else: # for all the other input files let's check the file content too
                if os.path.getsize(path) == 0:
                    logging.error(f"{message_prefix} [{filename}] is empty. Please check the file at default path [{path}].")
                    sys.exit() # not with error code 1 (successive parts will call exit code 1 in case this file was really necessary) 
                if os.path.getsize(path) > 0:
                    argument_path = path
                    return argument_path
    if argument_path == "":
        return argument_path


def get_argument_path(argument_path, default_path, message_prefix):
    """
    Checks if the user has provided a path for the input file, otherwise the function will look for the file
    in the default path where the pipeline creates it.

    :param argument_path: the optional argument given to the script for the input file position (can be user-defined or is "None" by default)
    :param default_path: the path at which the file is by default created by the pipeline
    :param message_prefix: a short definition of the file to be printed in the error message
    :return: the path at which the file has been found (either user-defined or default path)
    """
    if argument_path: # search in user-defined path
        check_inputfile(argument_path, message_prefix)
        return argument_path
    else: # search in default path
        argument_path = check_file_existence_and_content_in_default_paths(default_path, message_prefix)
        return argument_path