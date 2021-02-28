import os
import sys
import logging
import datetime
import subprocess
import shutil
import math
import tempfile

def init_logging(log_heading, logging_level):
    """
    Initialize the logging environment and print a header to the log.

    :param log_heading: heading for the log
    :return: nothing
    """
    logging.basicConfig(format='%(levelname)s\t%(message)s', level=logging_level, stream=sys.stdout)
    length = math.ceil((len(log_heading)/2))
    logging.info('- ' * length)
    logging.info(log_heading)
    logging.info(datetime.datetime.today().ctime())
    logging.info('- ' * length)


def merge_dicts(dict1, dict2):
    """
    Merge two dictionaries into a new dictionary. Keys in the second dictionary
    overwrite existing keys in the first dictionary.

    :param dict1: first dictionary
    :param dict2: second dictionary
    :return: a new dictionary
    """
    new_dict = {}
    new_dict.update(dict1)
    new_dict.update(dict2)
    return new_dict


def concat_files(file1, file2, new_file):
    """
    Concatenates two files into a new file.

    :param file1: path to first file
    :param file2: path to second file
    :param new_file: path to new concatenated file
    :return: nothing
    """
    with open(new_file, 'wb') as nf:
        with open(file1, 'rb') as f1:
            shutil.copyfileobj(f1, nf)
        with open(file2, 'rb') as f2:
            shutil.copyfileobj(f2, nf)


def can_i_run_software(software):
    """
    Copied and modified from wgd.utils

    Test if external software is executable
    :param software: list or string of executable(s)
    :return: 1 (failure) or 0 (success)
    """
    if type(software) == str:
        software = [software]
    ex = 0
    for s in software:
        # codeml needs input otherwise it prompts the user for input, so a dummy
        # file is created
        if s == 'codeml':
            tmp_file = "codeml.ctl"
            command = ['codeml', tmp_file]
        elif s == 'prank':
            command = [s, '--help']
        elif s == 'FastTree':
            command = s
        elif s in ['blastp', 'makeblastdb', 'blast', 'muscle', 'i-adhore']:
            command = [s, '-version']
        else:
            command = [s, '--version']
        try:
            # logging.info(command)
            if s == "codeml":
                # Since the Nextflow pipeline processes multiple wgd runs at the same time,
                # let's generate for each run the dummy "codeml.ctl" files within a different
                # temporary directory named with a unique ID. This way, the several parallel
                # wgd runs will not interfere with each other when creating and removing 
                # the dummy codeml.ctl file.
                with tempfile.TemporaryDirectory(dir = ".") as tmp_dir:
                    with open(os.path.join(tmp_dir, tmp_file), 'w') as o:  # write the codeml.ctl file in it
                        o.write('seqfile = test')
                    sp = subprocess.run(command, cwd=tmp_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    # tmp_dir is removed both if the subprocess succeeds or fails, thanks to the "with" statement
            else:
                sp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            out = sp.stdout.decode("utf-8").strip()
            err = sp.stderr.decode("utf-8").strip()
            if out:
                logging.info(out.splitlines()[0])
            elif err:
                logging.info(err.splitlines()[0])

        except FileNotFoundError:
            logging.error('{} executable not found!'.format(s))
            ex = 1
    return ex


def translate_cds(sequence_dict, skip_invalid=False):
    """
    Copied and modified from wgd.utils

    Just another CDS to protein translater. Will give warnings when in-frame
    stop codons are found, invalid codons are found, or when the sequence length
    is not a multiple of three. Will translate the full sequence or until an
    unspecified or in-frame codon is encountered.
    :param sequence_dict: dictionary with gene IDs and CDS sequences
    :param skip_invalid: bool, skip invalid CDS? (default translates to first
        stop codon or end)
    :return: dictionary with gene IDs and proteins sequences
    """
    # TODO I should just use the Biopython translator
    aa_dict = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '',  'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '',  'TGG': 'W',
    }
    protein_dict = {}

    j = 0
    total = 0
    for key, val in sequence_dict.items():
        j += 1
        aa_seq = ''
        if len(val) % 3 != 0:
            logging.warning('Sequence length != multiple of 3 for {}!'.format(key))
            total += 1
        invalid = False
        for i in range(0, len(val), 3):
            if val[i:i + 3] not in aa_dict.keys():
                logging.warning('Invalid codon {0:>3} in {1}'.format(val[i:i+3], key))
                invalid = True
                total += 1
                break
            else:
                if aa_dict[val[i:i + 3]] == '' and i+3 != len(val):
                    logging.warning('In-frame STOP codon in {0} at position {1}:{2}'.format(key, i, i+3))
                    invalid = True
                    total += 1
                    break
                aa_seq += aa_dict[val[i:i + 3]]
        if invalid and skip_invalid:
            continue
        protein_dict[key] = aa_seq

    if total:
        logging.warning("There were {} warnings during translation".format(total))
    return protein_dict


def write_fasta(seq_dict, output_file, id_prefix=None, append=False):
    """
    Copied and modified from wgd.utils

    Write/append a sequence dictionary to a fasta file.
    :param seq_dict: sequence dictionary, see :py:func:`read_fasta`
    :param output_file: output file name
    :param id_prefix: prefix to add to gene IDs
    :param append: append to file
    :return: nothing
    """
    mode = 'w'
    if append:
        mode = 'a'
    with open(output_file, mode) as o:
        for key, val in seq_dict.items():
            if id_prefix and id_prefix != '':
                o.write('>' + id_prefix + '|' + key + '\n')
            else:
                o.write('>' + key + '\n')
            o.write(val + '\n')
