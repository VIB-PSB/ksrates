"""
--------------------------------------------------------------------------------

Copyright (C) 2018 Arthur Zwaenepoel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contact: arzwa@psb.vib-ugent.be

--------------------------------------------------------------------------------
"""
# Note that the reweighting method is still performed here for compatibility
# reasons. However the recomputed weights are not used in other parts of the
# wgd suite!

# Ideas:
#  - add functionality to use custom alignments?
#  - add more advanced alignment stripping?

# IMPORTS
from .codeml import Codeml
from .alignment import prepare_aln, align, get_pairwise_alns
from .utils import process_gene_families, get_sequences, write_fasta, read_fasta
from .phy import run_phyml, phylogenetic_tree_to_cluster_format, run_fasttree
from .phy import average_linkage_clustering
from operator import itemgetter
from joblib import Parallel, delayed
import numpy as np
import pandas as pd
import os
import shutil
import sys
import logging
import subprocess as sp
import glob
import traceback


# maximum 1% of the paralog gene families can fail, e.g. 30 out of 3000
_THRESHOLD_FAILED_PARALOG_GF = 0.01
# maximum 5% of the ortholog gene families can fail, e.g. 400 out of 8000
_THRESHOLD_FAILED_ORTHOLOG_GF = 0.05


class GeneFamiliesError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


# HELPER FUNCTIONS -------------------------------------------------------------
def _get_nucleotide_sequences(family, nucleotide):
    """
    Get nucleotide sequences for codon alignment with PRANK

    :param family: sequence dictionary for gene family
    :param nucleotide: nucleotide sequences
    :return: nucleotide sequences for gene family
    """
    out = {}
    for k, v in family.items():
        if len(nucleotide[k]) % 3 != 0:
            logging.warning("Sequence length not a multiple of 3!")
        out[k] = nucleotide[k]
    return out


def _weighting(pairwise_estimates, msa=None, method='alc'):
    """
    Wrapper for different weighting methods. The fastest method is average
    linkage clustering based on Ks values, other methods use phylogenetic trees
    (phyml or fasttree)

    :param pairwise_estimates: results dictionary from
        :py:meth:`codeml.Codeml.run()` output
    :param msa: protein multiple sequence alignment file path (optional,
        not necessary for method=``alc``)
    :param method: method, one of ``alc|phyml|fasttree``
    :return: clustering data structure and pairwise leaf distances
        (not if method=``alc``)
    """
    if pairwise_estimates is None:
        return None, None, None

    if pairwise_estimates['Ks'].shape[0] < 2:
        return None, None, None

    pairwise_distances = None
    tree_path = None

    if method == 'phyml':
        # PhyML tree construction
        logging.debug('Constructing phylogenetic tree with PhyML')
        tree_path = run_phyml(msa)
        clustering, pairwise_distances = phylogenetic_tree_to_cluster_format(
                tree_path, pairwise_estimates['Ks'])

    elif method == 'fasttree':
        # FastTree tree construction
        logging.debug('Constructing phylogenetic tree with FastTree')
        tree_path = run_fasttree(msa)
        clustering, pairwise_distances = phylogenetic_tree_to_cluster_format(
                tree_path, pairwise_estimates['Ks'])

    else:
        # Average linkage clustering based on Ks
        logging.debug('Performing average linkage clustering on Ks values.')
        clustering = average_linkage_clustering(pairwise_estimates['Ks'])

    logging.debug('Clustering used for weighting: \n{}'.format(str(clustering)))
    return clustering, pairwise_distances, tree_path


def _calculate_weights(clustering, pairwise_estimates, pairwise_distances=None):
    """
    This is a patch for weight calculation in the pairwise approach

    :param clustering: clustering array
    :param pairwise_estimates: pairwise Ks estimates array
    :param pairwise_distances: pairwise distances array
    :return: data frame
    """
    # None -> None
    if pairwise_estimates is None or clustering is None:
        return None

    # process the clustering structure to get weights
    leaves = pairwise_estimates.shape[0]
    nodes = {i: [i] for i in range(leaves)}
    weights = {}
    out = set()

    for x in range(clustering.shape[0]):
        node_1, node_2, distance = \
            clustering[x, 0], clustering[x, 1], clustering[x, 2] * 2
        grouping_node = leaves + x
        nodes[grouping_node] = nodes[node_1] + nodes[node_2]

        # get for every pair the weight with outliers included and detect
        # outliers
        for i in nodes[node_1]:
            for j in nodes[node_2]:
                if pairwise_distances:
                    distance = pairwise_distances[i][j]
                pair = '__'.join(sorted([
                    pairwise_estimates.index[j], pairwise_estimates.index[i]
                ]))
                weights[pair] = {
                    'WeightOutliersIncluded': 0,
                    'WeightOutliersExcluded': 0,
                    'Distance': distance,
                    'Node': grouping_node
                }

    return pd.DataFrame.from_dict(weights, orient='index')


def _calculate_weighted_ks(clustering, pairwise_estimates,
                           pairwise_distances=None, family_id=None):
    """
    Calculate weighted Ks, Kn and w values following the procedure as outlined
    in Vanneste et al. (2013)

    :param clustering: clustering results as produced with
        :py:func:`_weighting`, can be by average linkage clustering or
        phylogenetic means (phyml, fasttree)
    :param pairwise_estimates: pairwise Ks, Kn and w estimates as produced by
        :py:meth:`codeml.Codeml.run_codeml`
    :return: a pandas data frame of weighted Ks, Kn and w values.
    """
    # None -> None
    if pairwise_estimates is None or clustering is None:
        return None

    # process the clustering structure to get weights
    leaves = pairwise_estimates['Ks'].shape[0]
    nodes = {i: [i] for i in range(leaves)}
    weights = {}
    out = set()

    for x in range(clustering.shape[0]):
        node_1, node_2, distance = clustering[x, 0], clustering[x, 1], \
                                   clustering[x, 2] * 2
        grouping_node = leaves + x
        nodes[grouping_node] = nodes[node_1] + nodes[node_2]

        for i in nodes[node_1]:
            for j in nodes[node_2]:
                if pairwise_distances:
                    distance = pairwise_distances[i][j]
                p1 = pairwise_estimates['Ks'].index[i]
                p2 = pairwise_estimates['Ks'].index[j]
                pair = "__".join(sorted([p1, p2]))
                weights[pair] = [
                    p1, p2, family_id.split("__")[-1],
                    pairwise_estimates['Ks'].iloc[i, j],
                    pairwise_estimates['Ka'].iloc[i, j],
                    pairwise_estimates['Omega'].iloc[i, j],
                    distance, grouping_node
                ]

                if pairwise_estimates['Ks'].iloc[i, j] > 5:
                    out.add(grouping_node)

    df = pd.DataFrame.from_dict(weights, orient='index')
    df.columns = ['Paralog1', 'Paralog2', 'Family',
                  'Ks', 'Ka', 'Omega', 'Distance', 'Node']

    return df


def add_alignment_stats(df, stats, l, l_stripped):
    """
    Add alignment statistics to the data frame

    :param df: pandas data frame
    :param stats: stats dict see :py:func:`alignment.pairwise_alignment_stats`
    :param l: alignment length
    :param l_stripped: stripped alignment length
    :return: data frame
    """
    identity = []
    coverage = []
    for row in df.index:
        paralogs = sorted([df.loc[row]['Paralog1'], df.loc[row]['Paralog2']])
        identity.append(stats[paralogs[0]][paralogs[1]][0])
        coverage.append(stats[paralogs[0]][paralogs[1]][1])
    df['AlignmentIdentity'] = pd.Series(identity, index=df.index)
    df['AlignmentCoverage'] = pd.Series(coverage, index=df.index)
    df['AlignmentLength'] = pd.Series([l] * len(df.index), index=df.index)
    df['AlignmentLengthStripped'] = pd.Series([l_stripped] * len(df.index),
                                              index=df.index)
    return df


def add_alignment_stats_(df, stats):
    st = pd.DataFrame.from_dict(stats, orient='index')
    df_ = pd.merge(df, st, left_index=True, right_index=True)
    return df_


# WRAP ANALYSE FAMILY FUNCTIONS IN TRY/EXCEPT ---------------------------------
def analyse_family_try_except(analysis_function, threshold_failed,
                              family_id, family, nucleotide, tmp='./',
                              codeml='codeml', preserve=False, times=1,
                              min_length=100, method='alc', aligner='muscle',
                              output_dir='./out', n_families=-1,
                              last_family_id=None):
    """
    Wrapping the analysis in a try/except checkpoint, so that if one family has
    a problem the whole loop is not interrupted and only that family will be
    ignored.
    
    :param analysis_function: the specific analysis function that has to be
        wrapped. It can be analyse_family or analyse_family_pairwise; they both
        accept the same argument list.
    """
    if len(glob.glob("*.exception")) / n_families > threshold_failed:
        logging.error(f"Too many gene family analyses failed, "
                      f"terminating threads...")
        raise GeneFamiliesError("error")

    is_last_family = (family_id == last_family_id)
    try:
        analysis_function(family_id, family, nucleotide, tmp, codeml, preserve,
                          times, min_length, method, aligner, output_dir,
                          n_families, is_last_family)
    except Exception:
        # this also prints the exception:
        logging.exception(f"Unexpected internal error during analysis of gene "
                          f"family {family_id}:")
        logging.error("Skipping gene family")
        # create a marker file and write exception to it as info
        failed_family_file = os.path.join(tmp, family_id + '.exception')
        with open(failed_family_file, 'w+') as o:
            traceback.print_exception(*sys.exc_info(), file=o)
    finally:
        if len(glob.glob("*.exception")) / n_families > threshold_failed:
            logging.error(f"Too many gene family analyses failed, "
                          f"terminating threads...")
            raise GeneFamiliesError("error")
        # once the last gene family finished start logging number of finished
        # gene families (even though the last gene family finished there could
        # still be earlier larger families running)
        last_family_ks_file = os.path.join(tmp, last_family_id + '.Ks')
        last_family_ks_exists = os.path.isfile(last_family_ks_file)
        last_family_codeml_file = os.path.join(tmp, last_family_id + '.codeml')
        last_family_codeml_exists = os.path.isfile(last_family_codeml_file)
        last_family_noks_file = os.path.join(tmp, last_family_id + '.noKs')
        last_family_noks_exists = os.path.isfile(last_family_noks_file)
        if (last_family_ks_exists and not last_family_codeml_exists) or \
                last_family_noks_exists:
            n_finished_families = len(glob.glob("*.Ks")) + \
                                  len(glob.glob("*.noKs")) + \
                                  len(glob.glob("*.exception"))
            # only print for the last family or the few last families to finish
            if is_last_family or n_families - n_finished_families <= 3:
                logging.info(f"Finished {n_finished_families}/{n_families} "
                             f"gene family analyses...")


# ANALYSE WHOLE FAMILY ---------------------------------------------------------
def analyse_family(
        family_id, family, nucleotide, tmp='./', codeml='codeml',
        preserve=False, times=1, min_length=100, method='alc', aligner='muscle',
        output_dir='./out', n_families=-1, is_last_family=False
):
    """
    Wrapper function for the analysis of one paralog family. Performs alignment
    with :py:meth:`alignment.MSA.run_aligner` and codeml analysis with
    :py:meth:`codeml.Codeml.run_codeml`. Subsequently also clustering with
    :py:func:`_average_linkage_clustering` is performed and weighted Ks, Kn and
    w values are calculated using :py:func:`_calculate_weighted_ks`.

    :param family_id: gene family id
    :param family: dictionary with sequences of paralogs
    :param nucleotide: nucleotide (CDS) sequences dictionary
    :param tmp: tmp directory
    :param codeml: codeml path
    :param preserve: preserve intermediate files
    :param times: number of times to perform ML estimation of Ks, Ka and omega
        values
    :param method: weighting method, from fast to slow: ``alc``, ``fasttree``,
        ``phyml``
    :param aligner: alignment program
    :param output_dir: output directory
    :return: ``csv`` file with results for the paralog family of interest
    """
    # pre-processing -----------------------------------------------------------
    if os.path.isfile(os.path.join(tmp, family_id + '.Ks')):
        logging.info('Found {}.Ks in tmp, will use this'.format(family_id))
        return
    if len(list(family.keys())) < 2:
        logging.debug("Skipping singleton gene family {}.".format(family_id))
        return

    logging.info('Performing analysis on gene family {} (size {})'
                 .format(family_id, len(list(family.keys()))))

    # if it is the last gene family start logging the number of finished gene
    # families (even though the last gene family finished there could still be
    # earlier larger families running)
    if is_last_family:
        n_finished_families = len(glob.glob("*.Ks")) + \
                              len(glob.glob("*.noKs")) + \
                              len(glob.glob("*.exception"))
        logging.info(f"Finished {n_finished_families}/{n_families} "
                     f"gene family analyses...")

    # multiple sequence alignment ----------------------------------------------
    if aligner == 'prank':
        # get nucleotide sequences instead of protein sequences
        logging.debug('Aligner is prank, will perform codon alignment')
        family = _get_nucleotide_sequences(family, nucleotide)
    logging.debug('Performing MSA ({0}) for {1}'.format(aligner, family_id))
    ff = write_fasta(family, os.path.join(tmp, family_id + '.fasta'))
    msa_path_protein = align(in_file=ff, out_file=ff + '.msa', aligner=aligner)
    msa_path, stats, successful = prepare_aln(msa_path_protein, nucleotide)
    if not successful:
        logging.warning("Failed to obtain codon alignment for gene family "
                        "{}".format(family_id))
        # TODO: why is the rest of the analysis not skipped in this case??

    # Calculate Ks values (codeml) ---------------------------------------------
    codeml = Codeml(codeml=codeml, tmp=tmp, id=family_id)
    logging.debug('Performing codeml analysis on {}'.format(family_id))
    results_dict, codeml_out = codeml.run_codeml(
            os.path.basename(msa_path), preserve=preserve, times=times)
    if not results_dict or len(results_dict['Ks'].index) == 0 \
            or results_dict['Ks'].empty:
        logging.warning('No codeml results for gene family {}'.format(family_id))
        # create a marker file and skip the rest of analysis
        with open(os.path.join(tmp, family_id + '.noKs'), 'w+') as o:
            o.write('No codeml results')
        return

    # Subdivide families -------------------------------------------------------

    # Calculate weights according to method ------------------------------------
    if len(list(family.keys())) == 2 and method == "phyml":
        # phyml breaks when only two genes are in a family
        logging.debug("PhyML breaks with only two genes, do ALC instead.")
        logging.debug("Distance will be in Ks units!")
        clustering, pairwise_distances, tree_path = _weighting(
                results_dict, msa=msa_path_protein, method="alc")
    else:
        clustering, pairwise_distances, tree_path = _weighting(
                results_dict, msa=msa_path_protein, method=method)
    if clustering is not None:
        out = _calculate_weighted_ks(
                clustering, results_dict, pairwise_distances, family_id
        )
        out = add_alignment_stats_(out, stats)
        logging.debug(out)
        out.to_csv(os.path.join(tmp, family_id + '.Ks'))
    else:
        logging.warning('No {} clustering/weighting results for '
                        'gene family {}'.format(method, family_id))
        # create a marker file and skip the rest of analysis
        with open(os.path.join(tmp, family_id + '.noKs'), 'w+') as o:
            o.write('No clustering/weighting results')
        return

    # preserve or remove data --------------------------------------------------
    if preserve:
        shutil.move(msa_path, os.path.join(output_dir, 'msa'))
        ### sp.run(['mv', msa_path, os.path.join(output_dir, 'msa')])
        shutil.move(codeml_out, os.path.join(output_dir, 'codeml'))
        ### sp.run(['mv', codeml_out, os.path.join(output_dir, 'codeml')])
        if tree_path:
            shutil.move(tree_path, os.path.join(output_dir, 'trees'))
            ### sp.run(['mv', tree_path, os.path.join(output_dir, 'trees')])
    else:
        os.remove(msa_path)
        ### sp.run(['rm', msa_path])
        os.remove(codeml_out)
        ### sp.run(['rm', codeml_out])
        if tree_path:
            os.remove(tree_path)
            ### sp.run(['rm', tree_path])


# PAIRWISE ANALYSIS ------------------------------------------------------------
def analyse_family_pairwise(
        family_id, family, nucleotide, tmp='./', codeml='codeml',
        preserve=False, times=1, min_length=100, method='alc',
        aligner='muscle', output_dir='./out', n_families=-1,
        is_last_family=False
):
    """
    Perform Ks analysis for one gene family using teh pairwise approach. The
    pairwise analysis approach is the one most commonly used in the literature,
    e.g. Vanneste et al. (2013), Vanneste et al. (2015), Barker et al. (2008),
    Li et al. (2015), Maere et al. (2005) and many more. Implementation details
    differ for many publications, and their are several implementation choices
    that are non-trivial. A main point of importance is the weighting of
    multiple pairwise Ks estimates based on a phylogeny (or proxy thereof), with
    some authors that do not weight (naive pairwise Ks distributions). Some
    authors use phylogenetic trees, oter use hierarchical clustering for the
    weighting. Lastly some authors perform reweighting after outlier removal,
    others don't. Here the approach of Vanneste et al. (2013) and later papers
    is largely followed:

    (1) A multiple sequence alignment at the protein level is constructed for a
    full paralogous family.

    (2) For every pair of sequences in the multiple sequence alignment (i.e.
    n*(n-1)/2 pairs for a family of n members), strip all gaps in the pairwise
    alignment.

    (3) Then back-translate the sequences to obtain a codon alignment. Discard
    the sequence pair if the alignment is too short (min_length, default=100).

    (4) Use codeml to perform ML estimation of Ks, Ka and Ka/Ks. Here the
    empirical codon frequencies of the pairwise alignment with the F3x4 method.

    (5) A matrix is constructed with all pairwise estimates. Entries that do not
    have a valid Ks estimate (because the alignment was too short) are removed
    from the matrix. This is done by iteratively removing the columns (and rows)
    with most NaN values. These sequences are also removed from the aligment.

    (6) A weighting method is applied. Either a phylogenetic tree is constructed
    for the family or the Pairwise Ks matrix is clustered using average linkage
    clustering. Every Ks estimate for a duplication node is then weighted by the
    number of estimates that are present for that node such that the different
    estimates for the same duplication event result in  weight of 1 in total.

    (7) Outliers are removed, where outliers are naively defined as pairs with
    Ks estimates > 5. Outlier removal essentially consists in putting the weight
    of these pairs to 0 and recalculating the weights for the remaining pairs in
    the family.

    :param family_id: gene family ID
    :param family: gene family sequence dictionary
    :param nucleotide: nucleotide sequences
    :param tmp: tmp directory
    :param codeml: codeml executable
    :param preserve: preserve codeml, alignment and tree results
    :param times: number of times to perform codeml estimation
    :param min_length: minimum gap-stripped alignment length to consider
    :param method: weighting method
    :param aligner: alignment method
    :param output_dir: output directory
    :return: nada
    """
    # pre-processing -----------------------------------------------------------
    if os.path.isfile(os.path.join(tmp, family_id + '.Ks')):
        logging.info('Found {}.Ks in tmp, will use this'.format(family_id))
        return
    if len(list(family.keys())) < 2:
        logging.debug("Skipping singleton gene family {}.".format(family_id))
        return
    logging.info('Performing analysis on gene family {}'.format(family_id))

    # multiple sequence alignment ----------------------------------------------
    if aligner == 'prank':
        # get nucleotide sequences instead of protein sequences
        logging.debug('Aligner is prank, will perform codon alignment')
        family = _get_nucleotide_sequences(family, nucleotide)
    logging.debug('Performing MSA ({0}) for {1}'.format(aligner, family_id))
    ff = write_fasta(family, os.path.join(tmp, family_id + '.fasta'))
    msa_path = align(in_file=ff, out_file=ff + '.msa', aligner=aligner)
    alns, stats = get_pairwise_alns(msa_path, nucleotide, min_length)

    # for every pair perform codeml analysis -----------------------------------
    codeml_out_string = ''
    family_dict = {}
    ks_mat = {g: {h: np.nan for h in family.keys()} for g in family.keys()}
    for i, pair in enumerate(alns):
        pid, seqs = pair
        g1, g2 = pid.split('__')
        pairwise_msa = write_fasta(seqs, msa_path + '.' + str(i))
        codeml_ = Codeml(codeml=codeml, tmp=tmp, id=family_id + '.' + str(i))
        logging.debug('Performing codeml analysis for {}'.format(family_id))
        results_dict, codeml_out = codeml_.run_codeml(
                os.path.basename(pairwise_msa), preserve=preserve, times=times)
        if preserve:
            with open(codeml_out, 'r') as f:
                codeml_out_string += f.read()

        # store results
        family_dict[pid] = {
            'Paralog1': g1, 'Paralog2': g2, 'Family': family_id,
            'Ks': results_dict['Ks'][g1][g2],
            'Ka': results_dict['Ka'][g1][g2],
            'Omega': results_dict['Omega'][g1][g2],
        }
        family_dict[pid].update(stats[pid])

        # also keep a Ks distance matrix
        ks_mat[g1][g2] = results_dict['Ks'][g1][g2]
        ks_mat[g2][g1] = results_dict['Ks'][g1][g2]
        ks_mat[g1][g1] = 0
        ks_mat[g2][g2] = 0

        os.remove(codeml_out)
        os.remove(pairwise_msa)

    ks_mat = pd.DataFrame.from_dict(ks_mat)
    logging.debug('Ks matrix: \n{}'.format(ks_mat))

    family_df = pd.DataFrame.from_dict(family_dict, orient='index')

    # filter the Ks matrix -----------------------------------------------------
    # As in the pairwise method, not every pair necessarily has a Ks estimate
    # (since filtering is on pairs within families, not entire families). Some
    # genes should be filtered out, as they might occur in pairs which have a
    # ks estimate and pairs which have not got one. In principle, this filtering
    # is not necessary, unless the weighting method is average linkage
    # clustering, in which case an incomplete distance matrix would be attained.
    # The strategy here is to remove iteratively the column/row with the most NA
    # entries until no more NA entries are in the Ks matrix. This will give the
    # most parsimonious filtered matrix
    while sum(ks_mat.isnull().sum()) != 0:
        logging.debug('Filtering matrix')
        ks_mat.drop(ks_mat.isnull().sum().idxmax(), axis=0, inplace=True)
        ks_mat.drop(ks_mat.isnull().sum().idxmax(), axis=1, inplace=True)

    logging.debug('Filtered Ks matrix: \n{}'.format(ks_mat))

    # filter the alignment -----------------------------------------------------
    # remove the sequences that were filtered out from the Ks matrix from the
    # alignment as well
    aln = read_fasta(msa_path)
    to_del = [k for k in aln.keys() if k not in ks_mat.index]
    for k in to_del: del aln[k]
    write_fasta(aln, msa_path)

    # weighting ----------------------------------------------------------------
    if len(list(family.keys())) == 2 and method == "phyml":
        # phyml breaks when only two genes are in a family
        logging.debug("PhyML breaks with only two genes, do ALC instead.")
        logging.debug("Distance will be in Ks units!")
        clustering, pairwise_distances, tree_path = _weighting(
                {'Ks': ks_mat}, msa=msa_path, method="alc")
    else:
        clustering, pairwise_distances, tree_path = _weighting(
                {'Ks': ks_mat}, msa=msa_path, method=method)
    if clustering is None:
        logging.warning('No Ks estimates for {}'.format(family_id))
        return

    weights = _calculate_weights(clustering, ks_mat, pairwise_distances)
    pd.merge(
            family_df, weights, left_index=True, right_index=True, how='outer'
    ).to_csv(os.path.join(tmp, family_id + '.Ks'))

    # preserve -----------------------------------------------------------------
    if preserve:
        base = os.path.basename(msa_path)
        os.rename(msa_path, os.path.join(output_dir, 'msa', base))
        if tree_path:
            base = os.path.basename(tree_path)
            os.rename(tree_path, os.path.join(output_dir, 'trees', base))
        codeml_path = os.path.join(output_dir, 'codeml', family_id + '.codeml')
        with open(codeml_path, 'w') as f:
            f.write(codeml_out_string)
    else:
        os.remove(msa_path)
        if tree_path:
            os.remove(tree_path)

    return


def ks_analysis_one_vs_one(
        nucleotide_sequences, protein_sequences, orthologs, tmp_dir='./tmp',
        output_dir='./ks.out', codeml_path='codeml', aligner='muscle',
        preserve=True, times=1, n_threads=4
):
    """
    Calculate a Ks distribution for one vs. one orthologs.

    :param nucleotide_sequences: sequence dictionary
    :param protein_sequences: protein sequence dictionary
    :param orthologs: file with ortholog families
    :param tmp_dir: tmp directory
    :param output_dir: output directory
    :param codeml_path: path to codeml executable
    :param preserve: preserve intermediate results (muscle, codeml)
    :param times: number of times to perform codeml analysis
    :param aligner: aligner to use (muscle|prank)
    :param n_threads: number of CPU cores to use
    :return: data frame
    """
    # Filter families with one vs one orthologs for the species pair. ---
    orthologs = process_gene_families(orthologs, ignore_prefix=False)
    protein = get_sequences(orthologs, protein_sequences)

    # preserve intermediate data if asked --------------------------------------
    if preserve:
        # msa and codeml results
        if not os.path.isdir(os.path.join(output_dir, 'msa')):
            os.mkdir(os.path.join(output_dir, 'msa'))
        if not os.path.isdir(os.path.join(output_dir, 'codeml')):
            os.mkdir(os.path.join(output_dir, 'codeml'))

    n_families = len(protein.keys())

    # start analysis -----------------------------------------------------------
    logging.info('Started analysis of {} ortholog gene families in parallel '
                 'using {} threads'.format(n_families, n_threads))

    last_family = list(protein)[-1]
    try:
        # if a gene family analysis fails with an error it will be skipped
        # if too many analyses fail the loop exits with a GeneFamiliesError
        Parallel(n_jobs=n_threads)(delayed(analyse_family_try_except)(
            analyse_family, _THRESHOLD_FAILED_ORTHOLOG_GF,
            family, protein[family], nucleotide_sequences, tmp_dir,
            codeml_path, preserve, times, 0, 'alc', aligner,
            output_dir, n_families, last_family
        ) for family in protein.keys())
    except GeneFamiliesError:
        _log_failed_families_threshold(_THRESHOLD_FAILED_ORTHOLOG_GF,
                                       n_families, tmp_dir)
        sys.exit(1)
    else:
        logging.info(f"Finished all gene family analyses")
        _log_no_ks_families()

    logging.info('Analysis done')

    logging.info('Making results data frame')
    results_frame = pd.DataFrame(
            columns=['Paralog1', 'Paralog2', 'Family', 'Ks', 'Ka', 'Omega'])

    # count the number of analyzed pairs ---------------------------------------
    counts = 0
    for f in os.listdir(tmp_dir):
        if f[-3:] == '.Ks':
            counts += 1
            df = pd.read_csv(os.path.join(tmp_dir, f), index_col=0)
            results_frame = pd.concat([results_frame, df], sort=True)
    results_frame.index = list(range(len(results_frame.index)))

    # rename the index of the data_frame to gene1_gene2 (alphabetically) -------
    new_index = results_frame[['Paralog1', 'Paralog2']].apply(
            lambda x: '__'.join(sorted([str(y) for y in x])), axis=1)
    results_frame.index = new_index

    # adding weights for completeness
    results_frame = compute_weights(results_frame)

    logging.info('Removing tmp directory')
    shutil.rmtree(tmp_dir)

    return results_frame


def ks_analysis_paranome(
        nucleotide_sequences, protein_sequences, paralogs,
        tmp_dir='./tmp', output_dir='./ks.out', codeml_path='codeml',
        preserve=True, times=1, ignore_prefixes=False, n_threads=4,
        min_length=100, method='alc', aligner='muscle',
        pairwise=False, max_gene_family_size=142
):
    """
    Calculate a Ks distribution for a whole paranome.

    :param nucleotide_sequences: sequence dictionary
    :param protein_sequences: protein sequence dictionary
    :param paralogs: file with paralog families
    :param tmp_dir: tmp directory
    :param output_dir: output directory
    :param codeml_path: path to codeml executable
    :param preserve: preserve intermediate results (muscle, codeml)
    :param times: number of times to perform codeml analysis
    :param ignore_prefixes: ignore prefixes in paralog/gene family file
        (e.g. in ath|AT1G45000, ath| will be ignored)
    :param min_length: minimum MSA length
    :param method: method to use, from fast to slow: ``alc``, ``fasttree``,
        ``phyml``
    :param aligner: alignment program to use (muscle|prank)
    :param n_threads: number of CPU cores to use:
    :param pairwise: perform pairwise (instead of gene family wise) analysis
    :param max_gene_family_size: maximum number of members a gene family
        may have
    :return: data frame
    """
    # ignore prefixes in gene families, since only one species -----------------
    paralogs = process_gene_families(paralogs, ignore_prefix=ignore_prefixes)
    protein = get_sequences(paralogs, protein_sequences)

    # preserve intermediate data if asked --------------------------------------
    if preserve:
        # msa and codeml results
        if not os.path.isdir(os.path.join(output_dir, 'msa')):
            os.mkdir(os.path.join(output_dir, 'msa'))
        if not os.path.isdir(os.path.join(output_dir, 'codeml')):
            os.mkdir(os.path.join(output_dir, 'codeml'))
        if method in ['fasttree', 'phyml']:
            if not os.path.isdir(os.path.join(output_dir, 'trees')):
                os.mkdir(os.path.join(output_dir, 'trees'))

    # sort family ids by family size -------------------------------------------
    # NOTE: I changed this so that filtering is also performed in the
    # non-pairwise analysis.
    sorted_families = sort_families_by_size(protein, True,
                                            max_gene_family_size)
    n_families = len(sorted_families)
    if n_families == 0:
        logging.error("No non-singleton gene families provided.")
        exit()

    # start analysis -----------------------------------------------------------
    logging.info('Started analysis of {} gene families in parallel using '
                 '{} threads'.format(n_families, n_threads))
    if pairwise:
        analysis_function = analyse_family_pairwise
    else:
        analysis_function = analyse_family

    try:
        # if a gene family analysis fails with an error it will be skipped
        # if too many analyses fail the loop exits with a GeneFamiliesError
        Parallel(n_jobs=n_threads)(delayed(analyse_family_try_except)(
            analysis_function, _THRESHOLD_FAILED_PARALOG_GF,
            family[0], protein[family[0]], nucleotide_sequences, tmp_dir,
            codeml_path, preserve, times, min_length, method, aligner,
            output_dir, n_families, sorted_families[n_families-1][0]
        ) for family in sorted_families)
    except GeneFamiliesError:
        _log_failed_families_threshold(_THRESHOLD_FAILED_PARALOG_GF,
                                       n_families, tmp_dir)
        sys.exit(1)
    else:
        logging.info(f"Finished all gene family analyses")
        _log_no_ks_families()

    logging.info('Analysis done')

    logging.info('Making results data frame')
    results_frame = pd.DataFrame(
            columns=['Paralog1', 'Paralog2', 'Family', 'Ks', 'Ka', 'Omega']
    )

    # count the number of analyzed pairs and get data frame --------------------
    counts = 0
    for f in os.listdir(tmp_dir):
        if f[-3:] == '.Ks':
            counts += 1
            df = pd.read_csv(os.path.join(tmp_dir, f), index_col=0)
            results_frame = pd.concat([results_frame, df], sort=True)
    results_frame.index = list(range(len(results_frame.index)))

    # rename the index of the data_frame to gene1_gene2 (alphabetically) -------
    new_index = results_frame[['Paralog1', 'Paralog2']].apply(
            lambda x: '__'.join(sorted([str(y) for y in x])), axis=1)
    results_frame.index = new_index
    logging.info("Computing weights, outlier cut-off at Ks > 20")
    results_frame = compute_weights(results_frame)

    logging.info('Removing tmp directory')
    shutil.rmtree(tmp_dir)

    return results_frame


def _log_no_ks_families():
    noks_families = glob.glob("*.noKs")
    n_noresults_families = len(noks_families)
    if n_noresults_families != 0:
        logging.warning("--")
        if n_noresults_families == 1:
            logging.warning(f"The analysis of {n_noresults_families} gene "
                            f"family gave no Ks results:")
        else:
            logging.warning(f"The analyses of {n_noresults_families} gene "
                            f"families gave no Ks results:")
        if n_noresults_families < 8:
            family_ids = [f.split(sep='.')[0] for f in sorted(noks_families)]
            logging.warning(f"{', '.join(family_ids)}")
        else:
            for family in sorted(noks_families):
                logging.warning(family.split(sep='.')[0])
    failed_families = glob.glob("*.exception")
    n_failed_families = len(failed_families)
    if n_failed_families != 0:
        # Warn the user that a few gene families failed
        # (not enough to justify an exit, though)
        logging.warning("--")
        if n_failed_families == 1:
            logging.warning(f"The analysis of {n_failed_families} gene family "
                            f"has failed due to unexpected internal errors")
            logging.warning(f"You may want to check the traceback above for "
                            f"the following gene family ID:")
        else:
            logging.warning(f"The analyses of {n_failed_families} gene "
                            f"families have failed due to unexpected internal "
                            f"errors")
            logging.warning(f"You may want to check the tracebacks above for "
                            f"the following gene family IDs:")
        if n_failed_families < 8:
            family_ids = [f.split(sep='.')[0] for f in sorted(failed_families)]
            logging.warning(f"{', '.join(family_ids)}")
        else:
            for family in sorted(failed_families):
                logging.warning(family.split(sep='.')[0])
    if n_noresults_families != 0 or n_failed_families:
        logging.warning("--")


def _log_failed_families_threshold(threshold, n_families, tmp_dir):
    failed_families = glob.glob("*.exception")
    n_failed_families = len(failed_families)
    logging.error("--")
    logging.error(f"The analyses of more than {int(threshold * 100)}% of "
                  f"gene families [{n_failed_families}/{n_families}] have "
                  f"failed due to unexpected internal errors")
    logging.error(f"Please check the nature of the error(s), remove the tmp "
                  f"directory [{tmp_dir}] and rerun the Ks analysis")
    logging.error(f"See the tracebacks above for the following gene family "
                  f"IDs:")
    if n_failed_families < 8:
        family_ids = [f.split(sep='.')[0] for f in sorted(failed_families)]
        logging.error(f"{', '.join(family_ids)}")
    else:
        for family in sorted(failed_families):
            logging.error(family.split(sep='.')[0])
    logging.error("Exiting")


def sort_families_by_size(families, pairwise=False, max_gene_family_size=142):
    """
    Returns a list of non-singleton gene families ordered by size and apply some
    filters

    :param families: nested gene family dictionary {family: {gene: sequence}}
    :param pairwise: pairwise analysis
    :param max_gene_family_size: maximum number of members a gene family may have
    :return: list of tuples [(family id, size)] sorted by size
    """
    sorted_families = []
    for k, v in families.items():
        if len(v.keys()) > 1:
            sorted_families.append((k, len(v.keys())))

    sorted_families = sorted(sorted_families, key=itemgetter(1), reverse=True)
    n_families = len(sorted_families)

    # filter for pairwise analysis
    if pairwise:
        sorted_families = [
            x for x in sorted_families if x[1] <= max_gene_family_size
        ]
        n_excluded_families = n_families - len(sorted_families)
        if n_excluded_families == 1:
            logging.warning('Filtered out the largest gene family because its '
                            'size is > {}'.format(max_gene_family_size))
            logging.warning("If you want to analyse this large family anyhow, "
                            "please raise the `max_gene_family_size` parameter")
        elif n_excluded_families > 1:
            logging.warning('Filtered out the {} largest gene families because '
                            'their size is > {}'.format(n_excluded_families,
                                                        max_gene_family_size))
            logging.warning("If you want to analyse these large families "
                            "anyhow, please raise the `max_gene_family_size` "
                            "parameter")

    return sorted_families


def compute_weights(df, min_ks=0.005, max_ks=20, aln_id=0, aln_len=300,
        aln_cov=0):
    """
    Modified from wgd.
    max_ks changed from 5 to 20 so to have weights computed also for Ks greater
    than 5 in the WeightOutliersExcluded column.
    """
    df = df[~df.index.duplicated()]  # for safety
    df["WeightOutliersIncluded"] = 1 / df.groupby(['Family', 'Node'])[
        'Ks'].transform('count')
    df_ = df[df["Ks"] <= max_ks]
    df_ = df_[df_["Ks"] >= min_ks]
    df_ = df_[df_["AlignmentCoverage"] >= aln_cov]
    df_ = df_[df_["AlignmentIdentity"] >= aln_id]
    df_ = df_[df_["AlignmentLength"] >= aln_len]
    df["WeightOutliersExcluded"] = np.zeros(len(df.index))
    df.loc[df_.index, "WeightOutliersExcluded"] = 1 / df_.groupby(
            ['Family', 'Node'])['Ks'].transform('count')
    return df
