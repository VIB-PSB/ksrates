import os
import sys
import logging
import shutil
import subprocess
import pandas as pd
from numpy import zeros
from wgd.utils import read_fasta
from wgd.blast_mcl import run_mcl_ava, ava_blast_to_abc, get_one_v_one_orthologs_rbh
from wgd.ks_distribution import ks_analysis_paranome, ks_analysis_one_vs_one
from wgd.colinearity import gff_parser
from ksrates.utils import merge_dicts, concat_files, can_i_run_software, translate_cds, write_fasta

_OUTPUT_BLAST_FILE_PATTERN = '{}.blast.tsv'
_OUTPUT_MCL_FILE_PATTERN = '{}.mcl.tsv'
_OUTPUT_KS_FILE_PATTERN = '{}.ks.tsv'
_PARALOGS_OUTPUT_DIR_PATTERN = 'wgd_{}'
_ORTHOLOGS_OUTPUT_DIR_PATTERN = 'wgd_{}_{}'

_TMP_BLAST = '{}.blast_tmp'
_TMP_KS = '{}.ks_tmp'

def ks_paralogs(species_name, cds_fasta, base_dir='.', eval_cutoff=1e-10, inflation_factor=2.0, aligner='muscle',
                min_msa_length=100, codeml='codeml', codeml_times=1, pairwise=False,
                max_gene_family_size=200, weighting_method='fasttree', n_threads=4, overwrite=False):
    """
    Modified from wgd_cli.py

    All vs. all Blast and one-to-one ortholog delineation (reciprocal best hits)
    :param species_name: name or ID of the species, used in output file names
    :param cds_fasta: CDS fasta file
    :param base_dir: directory where the output directory will be created in
    :param eval_cutoff: e-value cut off for blastp analysis
    :param inflation_factor: inflation factor for MCL clustering
    :param aligner: aligner to use
    :param min_msa_length: minimum multiple sequence alignment length
    :param codeml: path to codeml executable for ML estimation of Ks
    :param codeml_times: number of times to iteratively perform codeml ML estimation of Ks
    :param pairwise: run in pairwise mode
    :param max_gene_family_size: maximum size of a gene family to be included in analysis
    :param weighting_method: weighting method (fasttree, phyml or alc)
    :param n_threads: number of threads to use
    :param overwrite: overwrite existing results
    :return: nothing
    """

    # input checks
    if not species_name:
        logging.error('No species name provided. Exiting.')
        sys.exit(1)

    output_blast_file = _OUTPUT_BLAST_FILE_PATTERN.format(species_name)
    output_mcl_file = _OUTPUT_MCL_FILE_PATTERN.format(species_name)
    output_ks_file = _OUTPUT_KS_FILE_PATTERN.format(species_name)

    output_dir = os.path.abspath(os.path.join(base_dir, _PARALOGS_OUTPUT_DIR_PATTERN.format(species_name)))

    tmp_blast_paralogs = os.path.join(output_dir, _TMP_BLAST.format(species_name))
    tmp_ks_paralogs = os.path.join(output_dir, _TMP_KS.format(species_name))

    # determine which parts to run
    do_blast = True
    do_mcl = True
    do_ks = True
    if os.path.exists(output_dir):
        if os.path.exists(os.path.join(output_dir, output_blast_file)):
                
            # Check if blast tmp folder of previous failed run exists (in this case the blast table is incomplete)
            if os.path.exists(tmp_blast_paralogs):
                logging.error(f'tmp directory [{tmp_blast_paralogs}] was found: leftover from a failed earlier run?')
                logging.error(f"Paralog blast data may be incomplete, please delete [{output_blast_file}] and the "
                              f"tmp directory and rerun the analysis. Exiting.")
                sys.exit(1)

            if overwrite:
                logging.warning(f'Paralog blast file {output_blast_file} exists, will overwrite')
            else:
                logging.info(f'Paralog blast data {output_blast_file} already exists, '
                             f'will skip wgd all versus all Blastp')
                do_blast = False
        else:
            logging.info('No paralog blast data, will run wgd all versus all Blastp')
        if os.path.exists(os.path.join(output_dir, output_mcl_file)):
            if overwrite or do_blast:
                logging.warning(f'Paralog gene families file {output_mcl_file} exists, will overwrite')
            else:
                logging.info(f'Paralog gene family data {output_mcl_file} already exists, '
                             f'will skip wgd mcl')
                do_mcl = False
        else:
            logging.info('No paralog gene family data, will run wgd mcl')

        # Check if Ks tmp folder of previous failed run exists (in this case the Ks data will be incomplete)
        if os.path.exists(tmp_ks_paralogs): # if there is a Ks tmp folder from previous failed run
            logging.error(f'tmp directory [{tmp_ks_paralogs}] was found: leftover from a failed earlier run?')
            logging.error('Paralog Ks data may be incomplete. Please delete the tmp directory and rerun '
                          'the analysis. Exiting.')
            sys.exit(1)

        if os.path.exists(os.path.join(output_dir, output_ks_file)):
            if overwrite or do_blast or do_mcl:
                logging.warning(f'Paralog Ks file {output_ks_file} exists, will overwrite')
            else:
                logging.info(f'Paralog Ks data {output_ks_file} already exists, will skip wgd Ks analysis')
                do_ks = False
        else:
            logging.info('No paralog Ks data, will run wgd Ks analysis')

    if not do_blast and not do_mcl and not do_ks:
        logging.info('All paralog data already exist, nothing to do')
        logging.info('Done')
        return

    # software checks
    logging.info('---')
    logging.info('Checking external software...')
    software = []
    if do_blast:
        software += ['makeblastdb', 'blastp']
    if do_mcl:
        software += ['mcl']
    if do_ks:
        software += [aligner, codeml]
        if weighting_method == 'fasttree':
            software += ['FastTree']
        elif weighting_method == 'phyml':
            software += ['phyml']
    if can_i_run_software(software) == 1:
        logging.error('Could not run all required external software. Exiting.')
        sys.exit(1)

    # input checks
    if not cds_fasta:
        logging.error('No CDS fasta file provided. Exiting.')
        sys.exit(1)

    if not os.path.exists(output_dir):
        logging.info(f'Creating output directory {output_dir}')
        os.makedirs(output_dir)

    # read and translate CDS file
    cds_sequences = None
    protein_sequences = None
    if do_blast or do_ks:
        logging.info(f'Translating CDS file {cds_fasta}...')
        cds_sequences = read_fasta(os.path.abspath(cds_fasta))
        protein_sequences = translate_cds(cds_sequences)

    # all vs. all blast
    if do_blast:
        logging.info('---')
        logging.info('Running all versus all Blastp')

        # tmp directory
        if not os.path.exists(tmp_blast_paralogs):
            os.mkdir(tmp_blast_paralogs)

        # all-vs.-all, i.e. species is both db and query
        logging.info(f'Writing protein Blastdb sequences to {tmp_blast_paralogs}/...')
        db = os.path.join(tmp_blast_paralogs, f'{species_name}.db.fasta')
        write_fasta(protein_sequences, db)
        logging.info(f'Writing protein query sequences to {tmp_blast_paralogs}/...')
        query = os.path.join(tmp_blast_paralogs, f'{species_name}.query.fasta')
        write_fasta(protein_sequences, query)

        # start the blast
        logging.info('Performing all versus all Blastp (this might take a while)...')
        all_v_all_blast(query, db, output_dir, output_file=output_blast_file, eval_cutoff=eval_cutoff,
                        n_threads=n_threads)

        # remove temporary files
        logging.info('Removing tmp directory')
        shutil.rmtree(tmp_blast_paralogs)

    # get paranome (MCL)
    if do_mcl:
        logging.info('---')
        logging.info(f'Running gene family construction (MCL clustering with inflation factor = {inflation_factor})')
        ava_graph = ava_blast_to_abc(os.path.join(output_dir, output_blast_file))
        run_mcl_ava(ava_graph, output_dir=output_dir, output_file=output_mcl_file, inflation=inflation_factor)

    # whole paranome Ks analysis
    if do_ks:
        logging.info('---')
        logging.info('Running whole paranome Ks analysis...')

        # tmp directory
        cw_dir = os.getcwd()
        if not os.path.exists(tmp_ks_paralogs):
            os.mkdir(tmp_ks_paralogs)
        os.chdir(tmp_ks_paralogs)  # change dir to tmp dir, as codeml writes non-unique file names to the working dir

        output_mcl_path = os.path.join(output_dir, output_mcl_file)
        results_df = ks_analysis_paranome(cds_sequences, protein_sequences, output_mcl_path, tmp_ks_paralogs,
                                          output_dir, codeml, times=codeml_times, aligner=aligner,
                                          ignore_prefixes=False, min_length=min_msa_length, pairwise=pairwise,
                                          max_gene_family_size=max_gene_family_size, method=weighting_method,
                                          preserve=False, n_threads=n_threads)
        if os.path.exists(tmp_ks_paralogs):
            logging.error(f'tmp directory {tmp_ks_paralogs} still exists after analysis, '
                          f'paralog Ks data may be incomplete and will not be written to paralog Ks file!')
            logging.error('Please delete the tmp directory and rerun the analysis. Exiting.')
            sys.exit(1)
        if isinstance(results_df, type(None)) or results_df.empty:
            logging.warning('No paralog Ks data computed, will write empty paralog Ks file!')
        with open(os.path.join(output_dir, output_ks_file), 'w+') as o:
            o.write(results_df.round(5).to_csv(sep='\t'))
        # change back to current directory as tmp dir got deleted and subsequent os.getcwd() may fail
        os.chdir(cw_dir)

    logging.info('---')
    logging.info('Done')


def ks_orthologs(species1, species2, cds_fasta1, cds_fasta2, base_dir='.', eval_cutoff=1e-10, aligner='muscle',
                 codeml='codeml', codeml_times=1, n_threads=4, overwrite=False):
    """
    Modified from wgd_cli.py

    All vs. all Blast and one-to-one ortholog delineation (reciprocal best hits)
    :param species1: name or ID of the first species, used in output file names
    :param species2: name or ID of the second species, used in output file names
    :param cds_fasta1: CDS fasta file of the first species
    :param cds_fasta2: CDS fasta file of the second species
    :param base_dir: directory where the output directory will be created in
    :param eval_cutoff: e-value cut off for blastp analysis
    :param aligner: aligner to use
    :param codeml: path to codeml executable for ML estimation of Ks
    :param codeml_times: number of times to iteratively perform codeml ML estimation of Ks
    :param n_threads: number of threads to use
    :param overwrite: overwrite existing results
    :return: nothing
    """

    # input checks
    if not species1 and not species2:
        logging.error('No ortholog species names provided. Exiting.')
        sys.exit(1)

    output_blast_file = f'{species1}_{species2}.blast.tsv'
    output_orthologs_file = f'{species1}_{species2}.orthologs.tsv'
    output_ks_file = f'{species1}_{species2}.ks.tsv'

    output_dir = os.path.abspath(os.path.join(base_dir, _ORTHOLOGS_OUTPUT_DIR_PATTERN.format(species1, species2)))

    tmp_blast_orthologs = os.path.join(output_dir, _TMP_BLAST.format(f"{species1}_{species2}"))
    tmp_ks_orthologs = os.path.join(output_dir, _TMP_KS.format(f"{species1}_{species2}"))

    # determine which parts to run
    do_blast = True
    do_rbh = True
    do_ks = True
    if os.path.exists(output_dir):
        if os.path.exists(os.path.join(output_dir, output_blast_file)):

            # Check if blast tmp folder of previous failed run exists (in this case the blast table is incomplete)
            if os.path.exists(tmp_blast_orthologs):
                logging.error(f'tmp directory [{tmp_blast_orthologs}] was found: leftover from a failed earlier run?')
                logging.error(f"Ortholog blast data may be incomplete, please delete [{output_blast_file}] and the "
                              f"tmp directory and rerun the analysis. Exiting.")
                sys.exit(1)

            if overwrite:
                logging.warning(f'Ortholog blast file {output_blast_file} exists, will overwrite')
            else:
                logging.info(f'Ortholog blast data {output_blast_file} already exists, '
                             f'will skip wgd all versus all Blastp(s)')
                do_blast = False
        else:
            logging.info('No ortholog blast data, will run wgd all versus all Blastp(s)')
        if os.path.exists(os.path.join(output_dir, output_orthologs_file)):
            if overwrite or do_blast:
                logging.warning(f'Ortholog pairs file {output_orthologs_file} exists, will overwrite')
            else:
                logging.info(f'Ortholog pairs data {output_orthologs_file} already exists, '
                             f'will skip wgd one-to-one ortholog detection')
                do_rbh = False
        else:
            logging.info('No ortholog pairs data, will run wgd one-to-one ortholog detection')

        # Check if Ks tmp folder of previous failed run exists (in this case the Ks data will be incomplete)
        if os.path.exists(tmp_ks_orthologs):  # if there is a Ks tmp folder from previous failed run
            logging.error(f'tmp directory [{tmp_ks_orthologs}] was found: leftover from a failed earlier run?')
            logging.error('Ortholog Ks data may be incomplete. Please delete the tmp directory and rerun '
                          'the analysis. Exiting.')
            sys.exit(1)

        if os.path.exists(os.path.join(output_dir, output_ks_file)):
            if overwrite or do_blast or do_rbh:
                logging.warning(f'Ortholog Ks file {output_ks_file} exists, will overwrite')
            else:
                logging.info(f'Ortholog Ks data {output_ks_file} already exists, will skip wgd Ks analysis')
                do_ks = False
        else:
            logging.info('No ortholog Ks data, will run wgd Ks analysis')

    if not do_blast and not do_rbh and not do_ks:
        logging.info('All ortholog data already exist, nothing to do')
        logging.info('Done')
        return

    # software checks
    logging.info('---')
    logging.info('Checking external software...')
    software = []
    if do_blast:
        software += ['makeblastdb', 'blastp']
    if do_ks:
        software += [aligner, codeml]
    if can_i_run_software(software) == 1:
        logging.error('Could not run all required external software. Exiting.')
        sys.exit(1)

    # input checks
    if not cds_fasta1 and not cds_fasta2:
        logging.error('No CDS fasta files provided for ortholog species. Exiting.')
        sys.exit(1)

    if not os.path.exists(output_dir):
        logging.info(f'Creating output directory {output_dir}')
        os.makedirs(output_dir)

    # read and translate CDS files
    cds_sequences1 = None
    cds_sequences2 = None
    protein_sequences1 = None
    protein_sequences2 = None
    if do_blast or do_ks:
        logging.info(f'Translating CDS file {cds_fasta1}...')
        cds_sequences1 = read_fasta(os.path.abspath(cds_fasta1))
        protein_sequences1 = translate_cds(cds_sequences1)
        logging.info(f'Translating CDS file {cds_fasta2}...')
        cds_sequences2 = read_fasta(os.path.abspath(cds_fasta2))
        protein_sequences2 = translate_cds(cds_sequences2)

    # reciprocal all vs. all blast
    if do_blast:
        logging.info('---')
        logging.info('Running all versus all Blastp(s)')

        # tmp directory
        if not os.path.exists(tmp_blast_orthologs):
            os.mkdir(tmp_blast_orthologs)

        all_v_all_blast_1way(species1, species2, protein_sequences1, protein_sequences2, tmp_blast_orthologs, output_dir,
                             output_blast_file, eval_cutoff, n_threads)
        # all_v_all_blast_2ways(species1, species2, protein_sequences1, protein_sequences2, tmp_blast_orthologs, 
        #                       output_dir, output_blast_file, eval_cutoff, n_threads)

        # remove temporary files
        logging.info('Removing tmp directory')
        shutil.rmtree(tmp_blast_orthologs)

    # get one-to-one orthologs (RBHs)
    if do_rbh:
        logging.info('---')
        logging.info('Running one-to-one ortholog detection (reciprocal best blast hits)')
        output_file = get_one_v_one_orthologs_rbh(os.path.join(output_dir, output_blast_file), output_dir)
        os.rename(output_file, os.path.join(output_dir, output_orthologs_file))

    # one-to-one ortholog Ks analysis
    if do_ks:
        logging.info('---')
        logging.info('Running one-to-one ortholog Ks analysis...')

        # tmp directory
        cw_dir = os.getcwd()
        if not os.path.exists(tmp_ks_orthologs):
            os.mkdir(tmp_ks_orthologs)
        os.chdir(tmp_ks_orthologs)  # change dir to tmp dir, as codeml writes non-unique file names to the working dir

        cds_sequences = merge_dicts(cds_sequences1, cds_sequences2)
        protein_sequences = merge_dicts(protein_sequences1, protein_sequences2)
        output_orthologs_path = os.path.join(output_dir, output_orthologs_file)
        results_df = ks_analysis_one_vs_one(cds_sequences, protein_sequences, output_orthologs_path, tmp_ks_orthologs,
                                            output_dir, codeml, times=codeml_times, aligner=aligner, preserve=False,
                                            n_threads=n_threads)
        if os.path.exists(tmp_ks_orthologs):
            logging.error(f'tmp directory {tmp_ks_orthologs} still exists after analysis, '
                          f'ortholog Ks data may be incomplete and will not be written to ortholog Ks file!')
            logging.error('Please delete the tmp directory and rerun the analysis. Exiting.')
            sys.exit(1)
        if isinstance(results_df, type(None)) or results_df.empty:
            logging.warning('No ortholog Ks data computed, will write empty ortholog Ks file!')
        with open(os.path.join(output_dir, output_ks_file), 'w+') as o:
            o.write(results_df.round(5).to_csv(sep='\t'))
         # change back to current directory as tmp dir got deleted and subsequent os.getcwd() may fail
        os.chdir(cw_dir)

    logging.info('---')
    logging.info('Done')


def all_v_all_blast_1way(species1, species2, protein_seq_dict1, protein_seq_dict2, tmp_dir, output_dir, output_file,
                         eval_cutoff, n_threads):
    """
    Perform an all-versus-all Blastp of protein_seq_dict1 + protein_seq_dict2 as
    query against protein_seq_dict1 + protein_seq_dict2 as database.

    :param species1: name or ID of the first species, used in output file names
    :param species2: name or ID of the second species, used in output file names
    :param protein_seq_dict1: protein sequence dictionary for first species
    :param protein_seq_dict2: protein sequence dictionary for second species
    :param tmp_dir: temporary directory
    :param output_dir: output directory
    :param output_file: output file name
    :param eval_cutoff: e-value cut off for blastp analysis
    :param n_threads: number of threads to use
    :return: nothing
    """
    logging.info(f'Writing protein Blastdb sequences to {tmp_dir}/...')
    db = os.path.join(tmp_dir, f'{species1}_{species2}.db.fasta')
    write_fasta(protein_seq_dict1, db, id_prefix=species1)
    write_fasta(protein_seq_dict2, db, id_prefix=species2, append=True)

    logging.info(f'Writing protein query sequences to {tmp_dir}/...')
    query = os.path.join(tmp_dir, f'{species1}_{species2}.query.fasta')
    write_fasta(protein_seq_dict1, query, id_prefix=species1)
    write_fasta(protein_seq_dict2, query, id_prefix=species2, append=True)

    # start the blast
    logging.info('Performing all versus all Blastp (this might take a while)...')
    all_v_all_blast(query, db, output_dir, output_file=output_file, eval_cutoff=eval_cutoff, n_threads=n_threads)


def all_v_all_blast_2ways(species1, species2, protein_seq_dict1, protein_seq_dict2, tmp_dir, output_dir, output_file,
                          eval_cutoff, n_threads):
    """
    Perform two all-versus-all Blastp of 1) protein_seq_dict1 against
    protein_seq_dict2 and 2) protein_seq_dict2 against protein_seq_dict1,
    i.e. a reciprocal Blastp. Concatenate the results into output_file.

    :param species1: name or ID of the first species, used in output file names
    :param species2: name or ID of the second species, used in output file names
    :param protein_seq_dict1: protein sequence dictionary for first species
    :param protein_seq_dict2: protein sequence dictionary for second species
    :param tmp_dir: temporary directory
    :param output_dir: output directory
    :param output_file: output file name
    :param eval_cutoff: e-value cut off for blastp analysis
    :param n_threads: number of threads to use
    :return: nothing
    """
    logging.info(f'Writing protein sequences to {tmp_dir}/...')
    protein_fasta1 = os.path.join(tmp_dir, f'{species1}.protein.fasta')
    write_fasta(protein_seq_dict1, protein_fasta1, id_prefix=species1)
    protein_fasta2 = os.path.join(tmp_dir, f'{species2}.protein.fasta')
    write_fasta(protein_seq_dict2, protein_fasta2, id_prefix=species2)
    # do the 1st blast
    logging.info(f'Performing 1st all versus all Blastp {species1} - {species2} (this might take a while)...')
    output_blast_file_tmp12 = os.path.join(tmp_dir, f'{species1}_{species2}.blast.tsv')
    all_v_all_blast(protein_fasta1, protein_fasta2, tmp_dir, output_file=output_blast_file_tmp12,
                    eval_cutoff=eval_cutoff, n_threads=n_threads)
    # do the 2nd blast
    logging.info('-')
    logging.info(f'Performing 2nd all versus all Blastp {species2} - {species1} (this might take a while)...')
    output_blast_file_tmp21 = os.path.join(tmp_dir, f'{species2}_{species1}.blast.tsv')
    all_v_all_blast(protein_fasta2, protein_fasta1, tmp_dir, output_file=output_blast_file_tmp21,
                    eval_cutoff=eval_cutoff, n_threads=n_threads)
    # concatenate two blast results file
    logging.info('Concatenating results')
    concat_files(output_blast_file_tmp12, output_blast_file_tmp21, os.path.join(output_dir, output_file))


def all_v_all_blast(query, db, output_directory='./', output_file='blast.tsv', eval_cutoff=1e-10, n_threads=4):
    """
    Modified from wgd/blast_mcl.py

    Perform all-versus-all Blastp.
    Runs a blast of ``query`` vs. ``db``.
    :param query: query sequences fasta file
    :param db: database sequences fasta file
    :param output_directory: output directory
    :param output_file: output file name
    :param eval_cutoff: e-value cut-off
    :param n_threads: number of threads to use for blastp
    :return: blast file path
    """
    logging.info("Making Blastdb")
    command = ['makeblastdb', '-in', db, '-dbtype', 'prot']
    logging.info(' '.join(command))
    try:
        sp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8', check=True)
        out = sp.stdout.strip()
        err = sp.stderr.strip()
        if out:
            logging.info("makeblastdb output:\n" + out)
        if err:
            logging.error("makeblastdb standard error output:\n" + err)
    except subprocess.CalledProcessError as e:
        logging.error(f"makeblastdb execution failed with return code: {e.returncode}")
        if e.stderr:
            logging.error("makeblastdb standard error output:\n" + e.stderr.strip())
        logging.error("Exiting.")
        sys.exit(e.returncode)

    logging.info("Running Blastp")
    command = ['blastp', '-db', db, '-query', query, '-evalue', str(eval_cutoff), '-outfmt', '6',
               '-num_threads', str(n_threads), '-out', os.path.join(output_directory, output_file)]
    logging.info(' '.join(command))
    try:
        sp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8', check=True)
        out = sp.stdout.strip()
        err = sp.stderr.strip()
        if out:
            logging.info("blastp output:\n" + out)
        if err:
            logging.error("blastp standard error output:\n" + err)
    except subprocess.CalledProcessError as e:
        logging.error(f"blastp execution failed with return code: {e.returncode}")
        if e.stderr:
            logging.error("blastp standard error output:\n" + e.stderr.strip())
        if e.returncode == -11:
            logging.error("This looks like a segmentation fault, you may want to look for and "
                          "remove any core dump files (core.*).")
            logging.error("Try to increase the memory and/or lower the number of threads.")
        # remove output file as it is likely incomplete or broken
        os.remove(os.path.join(output_directory, output_file))
        ### subprocess.run(['rm', os.path.join(output_directory, output_file)])
        # remove database and other stuff
        os.remove(db + '.phr')
        os.remove(db + '.pin')
        os.remove(db + '.psq')
        ### subprocess.run(['rm', db + '.phr', db + '.pin', db + '.psq'])
        logging.error("Exiting.")
        sys.exit(e.returncode)
    logging.info("All versus all Blastp done")

    # remove database and other stuff
    os.remove(db + '.phr')
    os.remove(db + '.pin')
    os.remove(db + '.psq')
    ### subprocess.run(['rm', db + '.phr', db + '.pin', db + '.psq'])

    return os.path.join(output_directory, output_file)


def ks_colinearity(species_name, gff_file, base_dir='.', gff_feature='mRNA', gff_gene_attribute='Parent', n_threads=1,
                   overwrite=False):
    """
    Modified from wgd_cli.py

    Co-linearity analysis with i-ADHoRe 3.0.
    :param species_name: name or ID of the species, used in output file names
    :param gff_file: GFF3 annotation file (see the annotation files on PLAZA as
        an example)
    :param base_dir: directory where the output directory will be created in
    :param gff_feature: keyword for entities of interest in the GFF file, e.g.
        'CDS' or 'mRNA'
    :param gff_gene_attribute: attribute key for the gene ID in the GFF (9th
        column), e.g. 'ID' or 'Parent'
    :param n_threads: number of threads to use
    :param overwrite: overwrite existing results
    :return: nothing
    """

    # input checks
    if not species_name:
        logging.error('No species name provided. Exiting.')
        sys.exit(1)

    output_dir = os.path.abspath(os.path.join(base_dir, _PARALOGS_OUTPUT_DIR_PATTERN.format(species_name)))
    output_file = f'{species_name}.ks_anchors.tsv'

    if os.path.exists(os.path.join(output_dir, output_file)):
        if overwrite:
            logging.warning(f'Colinearity anchor pair Ks file {output_file} exists, will overwrite')
        else:
            logging.info(f'Colinearity anchor pair Ks data {output_file} already exists, '
                         f'will skip wgd colinearity Ks analysis')
            logging.info('All colinearity data already exist, nothing to do')
            logging.info('Done')
            return
    else:
        logging.info('No colinearity anchor pair Ks data, will run wgd colinearity Ks analysis')

    # software check
    logging.info('Checking external software...')
    if can_i_run_software(['i-adhore']) == 1:
        logging.error('Could not run all required external software. Exiting.')
        sys.exit(1)

    # input checks
    if not gff_file:
        logging.error('No GFF file provided. Exiting.')
        sys.exit(1)

    mcl_file = os.path.join(output_dir, _OUTPUT_MCL_FILE_PATTERN.format(species_name))
    if not os.path.exists(mcl_file):
        logging.error(f'Paralog gene families file {mcl_file} does not exist. Exiting.')
        sys.exit(1)

    ks_file = os.path.join(output_dir, _OUTPUT_KS_FILE_PATTERN.format(species_name))
    if not os.path.exists(ks_file):
        logging.error(f'Paralog Ks file {ks_file} does not exist. Exiting.')
        sys.exit(1)

    # i-ADHoRe tmp directory
    tmp_dir = os.path.join(output_dir, f'{species_name}.ks_anchors_tmp')

    if not os.path.exists(output_dir):
        logging.info(f'Creating output directory {output_dir}')
        os.makedirs(output_dir)
    if not os.path.exists(tmp_dir):
        logging.info(f'Creating i-ADHoRe tmp directory {tmp_dir}')
        os.mkdir(tmp_dir)

    # parse the gff
    logging.info("Parsing GFF file")
    try:
        genome, all_genes = gff_parser(gff_file, feature=gff_feature, gene_attribute=gff_gene_attribute)
    except IndexError:
        logging.error(f'Invalid GFF file {gff_file}: number of columns != 9. Exiting.')
        sys.exit(1)

    if not genome:
        logging.error("Gene lists for i-ADHoRe are empty, check GFF feature and (gene) attribute configurations. "
                      "Exiting.")
        shutil.rmtree(tmp_dir)
        sys.exit(1)

    # generate necessary files for i-ADHoRe
    logging.info("Writing gene lists for i-ADHoRe")
    _write_gene_lists(genome, os.path.join(tmp_dir, 'i-adhore_gene_lists'))

    logging.info("Writing families file for i-ADHoRe")
    n_matched = _write_families_file(mcl_file, all_genes, os.path.join(tmp_dir, 'families.tsv'))
    if not n_matched:
        logging.error("None of the (gene) IDs in the MCL gene families file match (gene) IDs in the GFF file, "
                      "check GFF feature and (gene) attribute configurations and/or (gene) IDs between FASTA and GFF "
                      "files. Exiting.")
        shutil.rmtree(tmp_dir)
        sys.exit(1)

    logging.info("Writing i-ADHoRe configuration file")
    _write_config_iadhore('i-adhore_gene_lists', 'families.tsv', base_dir=tmp_dir, genome=species_name,
                          number_of_threads=n_threads)

    # run i-ADHoRe
    logging.info("Running i-ADHoRe 3.0...")
    _run_iadhore(os.path.join(tmp_dir, 'i-adhore.conf'))

    # Ks distribution for anchors
    logging.info("Constructing Ks distribution for anchors")
    anchor_points_file = os.path.join(tmp_dir, 'i-adhore_out', 'anchorpoints.txt')
    _write_anchor_pairs_ks(anchor_points_file, ks_file, os.path.join(output_dir, output_file))

    logging.info('Removing i-ADHoRe tmp directory')
    subfolder = os.path.join(output_dir, f"{species_name}_i-adhore")
    if not os.path.exists(subfolder):
        os.mkdir(subfolder)
    shutil.copy2(os.path.join(tmp_dir, 'i-adhore_out', 'anchorpoints.txt'), subfolder)
    shutil.copy2(os.path.join(tmp_dir, "i-adhore_out", "multiplicons.txt"), subfolder)
    shutil.copy2(os.path.join(tmp_dir, "i-adhore_out", "segments.txt"), subfolder)
    shutil.copy2(os.path.join(tmp_dir, "i-adhore_out", "list_elements.txt"), subfolder)
    shutil.copy2(os.path.join(tmp_dir, "i-adhore_out", "multiplicon_pairs.txt"), subfolder)
    shutil.rmtree(tmp_dir)

    logging.info("Done")


def _write_gene_lists(genome, output_dir='gene_lists'):
    """
    Modified from wgd/colinearity.py

    Write out the gene lists for i-ADHoRe. See i-ADHoRe manual for information.
    :param genome: genome dictionary from gff_parser
    :param output_dir: output directory
    """
    if not genome:
        logging.warning("Gene lists are empty!")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for k, v in genome.items():
        with open(os.path.join(output_dir, k + '.lst'), 'w') as o:
            for feature in v:
                o.write(feature[0] + feature[-1] + '\n')

def _write_families_file(families_mcl_file, all_genes, output_file='families.tsv', split_on_pipe=True):
    """
    Modified from wgd/colinearity.py (to fix issue with pipe in gene names)

    Write out gene families file as expected by i-ADHoRe (one gene per line)
    :param families_mcl_file: gene families file from MCL (one family per line)
    :param all_genes: set object with all genes for the focal species
    :param output_file: output file name
    :param split_on_pipe: boolean, split gene IDs in the MCL gene families file on '|' character and take first part
    :return: number of (gene) IDs from MCL gene families file matched to GFF file
             (zero will indicate some problem with (gene) IDs)
    """
    counter = 1
    genes_seen = set()

    with open(families_mcl_file, 'r') as f:
        with open(output_file, 'w') as o:
            for line in f:
                line = line.strip().split('\t')
                for gene in line:
                    if split_on_pipe:
                        gene = gene.split("|")[0].strip()
                    o.write(gene + '\t' + str(counter) + '\n')
                    genes_seen.add(gene)
                counter += 1

            # add genes not seen to the families file as singletons
            # I-ADHoRe throws an error when there are genes in the gene
            # lists that are not in the blast table/families file
            rest = all_genes - genes_seen
            for gene in list(rest):
                o.write(gene + '\t' + str(counter) + '\n')
                counter += 1
    return len(all_genes) - len(rest)

def _write_config_iadhore(gene_lists_dir, families_file, base_dir='.', genome='genome', config_file='i-adhore.conf',
                          output_path='i-adhore_out', gap_size=30, tandem_gap=None, cluster_gap=35,
                          max_gaps_in_alignment=None, q_value=0.75, prob_cutoff=0.01, anchor_points=3,
                          alignment_method='gg2', level_2_only='false', table_type='family',
                          multiple_hypothesis_correction='FDR', visualize_ghm='false', visualize_alignment='true',
                          number_of_threads=1):
    """
    Modified from wgd/colinearity.py

    Write out the config file for i-ADHoRe. See i-ADHoRe manual for information
    on parameter settings.
    :param gene_lists_dir: directory with gene lists per chromosome
    :param families_file: name of file with gene to family mapping
    :param base_dir: directory where files are located/saved
    :param config_file: name for the config file to generate
    :param genome: genome name
    :param output_path: i-ADHoRe output path name
    :param gap_size: see i-ADHoRe 3.0 documentation
    :param tandem_gap: see i-ADHoRe 3.0 documentation
    :param cluster_gap: see i-ADHoRe 3.0 documentation
    :param q_value: see i-ADHoRe 3.0 documentation
    :param prob_cutoff: see i-ADHoRe 3.0 documentation
    :param anchor_points: see i-ADHoRe 3.0 documentation
    :param alignment_method: see i-ADHoRe 3.0 documentation
    :param level_2_only: see i-ADHoRe 3.0 documentation
    :param table_type: see i-ADHoRe 3.0 documentation
    :param multiple_hypothesis_correction: see i-ADHoRe 3.0 documentation
    :param visualize_ghm: see i-ADHoRe 3.0 documentation
    :param visualize_alignment: see i-ADHoRe 3.0 documentation
    :return: nothing
    """
    with open(os.path.join(base_dir, config_file), 'w') as o:
        o.write(f'genome={genome}\n')

        for gene_list in sorted(os.listdir(os.path.join(base_dir, gene_lists_dir))):
            o.write(f'{gene_list[:-4]} {os.path.join(base_dir, gene_lists_dir)}/{gene_list}\n')

        o.write(f'blast_table={os.path.join(base_dir, families_file)}\n')
        o.write(f'output_path={os.path.join(base_dir, output_path)}\n')
        o.write(f'gap_size={gap_size}\n')
        if tandem_gap:
            o.write(f'tandem_gap={tandem_gap}\n')
        else:
            o.write(f'tandem_gap={gap_size/2}\n')
        o.write(f'q_value={q_value}\n')
        o.write(f'cluster_gap={cluster_gap}\n')
        if max_gaps_in_alignment:
            o.write(f'max_gaps_in_alignment={max_gaps_in_alignment}\n')
        else:
            o.write(f'max_gaps_in_alignment={cluster_gap}\n')
        o.write(f'prob_cutoff={prob_cutoff}\n')
        o.write(f'anchor_points={anchor_points}\n')
        o.write(f'alignment_method={alignment_method}\n')
        o.write(f'level_2_only={level_2_only}\n')
        o.write(f'table_type={table_type}\n')
        o.write(f'multiple_hypothesis_correction={multiple_hypothesis_correction}\n')
        o.write(f'visualizeGHM={visualize_ghm}\n')
        o.write(f'visualizeAlignment={visualize_alignment}\n')
        o.write(f'number_of_threads={number_of_threads}\n')

def _run_iadhore(config_file):
    """
    Modified from wgd/colinearity.py

    Run i-ADHoRe for a given config file.
    :param config_file: path to i-ADHoRe configuration file
    """
    command = ['i-adhore', config_file]
    logging.info(' '.join(command))
    try:
        sp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8', check=True)
        out = sp.stdout.strip()
        err = sp.stderr.strip()
        if out:
            logging.info("i-ADHoRe output:\n" + out)
        if err:
            # i-ADHoRe can fail without a proper exit/return code, but just with an error on stderr
            logging.error("i-ADHoRe execution failed with standard error output:\n" + err)
            logging.error("Exiting.")
            sys.exit(1)
    except subprocess.CalledProcessError as e:
        logging.error(f"i-ADHoRe execution failed with return code: {e.returncode}")
        if e.stderr:
            logging.error("i-ADHoRe standard error output:\n" + e.stderr.strip())
        logging.error("Exiting.")
        sys.exit(e.returncode)

    return

def compute_weights_anchor_pairs(df, min_ks=0.05, max_ks=20, aln_id=0, aln_len=300,
        aln_cov=0):
    """
    Modified from wgd.
    Computes the weights of anchor pair Ks estimates.
    
    :param min_ks: minimum Ks value considered (hard coded to 0.05 Ks)
    :param max_ks: maximum Ks value considered
    :param aln_id: minimum alignment identity considered
    :param aln_len: minimum alignment length (with gaps) considered
    :param aln_cov: minimum alignment coverage considered
    :return: dataframe with updated weights ("outliers excluded", i.e.
             weights are updated only for pairs which have alignment lenght greater
             than 300 and with a Ks value between 0.05 and max_ks)
    """
    df = df[~df.index.duplicated()]  # for safety
    df_ = df[df["Ks"] <= max_ks]
    df_ = df_[df_["Ks"] >= min_ks]
    df_ = df_[df_["AlignmentCoverage"] >= aln_cov]
    df_ = df_[df_["AlignmentIdentity"] >= aln_id]
    df_ = df_[df_["AlignmentLength"] >= aln_len]
    df["WeightOutliersExcluded"] = zeros(len(df.index))
    df.loc[df_.index, "WeightOutliersExcluded"] = 1 / df_.groupby(
            ['Family', 'Node'])['Ks'].transform('count')
    
    df = df.drop(columns=["WeightOutliersIncluded"])  # it's only paranome related and not used anyways

    return df

def _write_anchor_pairs_ks(anchor_points_file, ks_file, out_file='ks_anchors.tsv'):
    """
    Modified from wgd/colinearity.py

    Get anchor pairs and their corresponding Ks values and write to file
    :param anchor_points_file: anchorpoints.txt file from i-ADHoRe 3.0
    :param ks_file: Ks data file
    :param out_file: output file name
    :return: nothing
    """
    # anchor_points = pd.read_csv(anchor_points_file, sep='\t')
    # anchor_points = anchor_points[['gene_x', 'gene_y']]
    # anchor_points.drop_duplicates()
    # ids = anchor_points.apply(lambda x: '__'.join(sorted(x)), axis=1)
    # print(ids)
    # ids.drop_duplicates()
    # print(ids)
    #
    # ks_distribution = pd.read_csv(ks_file, index_col=0, sep='\t')
    # # if (gene) IDs in the Ks file contain pipes remove them (i-ADHoRe anchor (gene) IDs should not contain pipes)
    # index_list = list(ks_distribution.index.values)
    # if index_list[0].count('|') == 2:
    #     logging.info("Removing pipes from gene IDs")
    #     nonpiped_index_list = list(index_list)
    #     for i, e in enumerate(nonpiped_index_list):
    #         genes = e.split("__")
    #         gene1 = genes[0].split("|")[0].strip()
    #         gene2 = genes[1].split("|")[0].strip()
    #         nonpiped_index_list[i] = f"{gene1}__{gene2}"
    #     ks_distribution = ks_distribution.rename(index=dict(zip(index_list, nonpiped_index_list)))
    #
    # ks_anchors = ks_distribution.loc[ks_distribution.index.intersection(ids)]
    #
    # if out_file:
    #     ks_anchors.to_csv(out_file, sep='\t')

    with open(anchor_points_file, "r") as apf, open(ks_file, "r") as ksf:
        anchor_points = pd.read_csv(apf, sep="\t")
        anchor_points = anchor_points[['gene_x', 'gene_y']]
        anchor_points.drop_duplicates()
        anchor_points_id_list = []
        for __, row in anchor_points.iterrows():
            anchor_point_pair = [row["gene_x"], row["gene_y"]]
            anchor_point_pair.sort()
            anchor_point_id = "__".join(anchor_point_pair)
            if not anchor_point_id in anchor_points_id_list:
                anchor_points_id_list.append(anchor_point_id)

        ks = pd.read_csv(ksf, index_col=0, sep='\t')
        # if (gene) IDs in the Ks file contain pipes remove them (i-ADHoRe anchor (gene) IDs should not contain pipes)
        index_list = list(ks.index.values)
        if isinstance(index_list[0], str) and index_list[0].count('|') == 2:
            logging.info("Removing pipes from gene IDs")
            nonpiped_index_list = list(index_list)
            for i, e in enumerate(nonpiped_index_list):
                genes = e.split("__")
                gene1 = genes[0].split("|")[0].strip()
                gene2 = genes[1].split("|")[0].strip()
                nonpiped_index_list[i] = f"{gene1}__{gene2}"
            ks = ks.rename(index=dict(zip(index_list, nonpiped_index_list)))

        ks_anchors = ks.loc[ks.index.intersection(pd.Series(anchor_points_id_list))]

        # (Re)calculate weights of anchor pairs
        ks_anchors_weighted = compute_weights_anchor_pairs(ks_anchors)

        if out_file:
            with open(out_file, "w+") as of:
                of.write(ks_anchors_weighted.to_csv(sep='\t'))
