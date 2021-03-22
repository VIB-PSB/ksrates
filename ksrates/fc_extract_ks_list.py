from pandas import read_csv
from numpy import float64, zeros


def ks_list_from_tsv(tsv_file, max_ks, data_type):
    """
    Extracts either whole-paranome Ks, anchor pairs Ks or ortholog Ks from a wgd Ks TSV file.
    If dealing with whole-paranome Ks values, it returns also a list of weights extracted from "WeightOutliersExcluded" column in the tsv_file.

    :param tsv_file: wgd output file containing either paranome, anchor pairs or ortholog Ks values (suffix formats: ".ks.tsv", "ks_anchors.tsv")
    :param max_ks: maximum Ks value to be accepted for the analysis
    :param data_type: specifies the nature of the Ks values in the tsv file (either "paralogs", "anchor pairs" or "orthologs")
    :return: ks_list_filtered, list of either paralog, anchors pairs or ortholog Ks values
    :return: weight_list_filtered, (returned only for paralog or ortholog data type) list of weights for Ks values        
    """
    with open(tsv_file, "r") as tsv_file:
        tsv = read_csv(tsv_file, sep="\t")

    if data_type == "paralogs" or data_type == "orthologs":
        filtered_tsv = tsv.loc[(tsv["Ks"].dtypes == float64) & (tsv["Ks"] <= max_ks)]
    elif data_type == "anchor pairs":
        filtered_tsv = tsv.loc[(tsv["Ks"].dtypes == float64) & (tsv["Ks"] >= 0.05) & (tsv["Ks"] <= max_ks)]

    ks_list_filtered = filtered_tsv["Ks"].to_list()
    weight_list_filtered = filtered_tsv["WeightOutliersExcluded"].to_list()

    if data_type == "paralogs" or data_type == "anchor pairs":
        return ks_list_filtered, weight_list_filtered
    if data_type == "orthologs":
        return ks_list_filtered


def compute_weights_anchor_pairs(df, min_ks=0.05, max_ks=20, aln_id=0, aln_len=300,
        aln_cov=0):
    """
    Modified from wgd.
    Computes the weights of anchor pair Ks estimates.
    
    :param min_ks: minimum Ks value considered (hard coded to 0.5 Ks)
    :param max_ks: maximum Ks value considered
    :param aln_id: minimum alignment identity considered
    :param aln_len: minimum alignment length (with gaps) considered
    :param aln_cov: minimum alignment coverage considered
    :return: dataframe with updated weights (outliers excluded)
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
