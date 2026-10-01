from pandas import read_csv
from numpy import float64
from wgd_ksrates.viz import filter_compute_weights


def ks_list_from_tsv(tsv_file, max_ks, data_type, min_ks=0.005):
    """
    Extracts paralog Ks (whole-paranome, anchor pairs or reciprocally retained) or ortholog Ks from a wgd Ks TSV file.
    If dealing with paralog Ks values, it returns also a list of weights extracted from "WeightOutliersExcluded" column in the tsv_file.

    :param tsv_file: wgd output file containing either paralog or ortholog Ks values (suffix formats: ".ks.tsv", "ks_anchors.tsv", "ks_recret_topX.tsv")
    :param max_ks: maximum Ks value to be accepted for the analysis
    :param data_type: specifies the nature of the Ks values in the tsv file (either "paralogs", "anchor pairs", "reciprocally retained" or "orthologs")
    :param min_ks: minimum Ks value to be accepted for the analysis
    :return: ks_list_filtered, list of either whole-paranome, anchors pairs, reciprocally retained or ortholog Ks values
    :return: weight_list_filtered, list of weights for Ks values (not returned for orthologs Ks)
    """
    with open(tsv_file, "r") as f:
        tsv = read_csv(f, sep="\t")
    return ks_list_from_df(tsv, max_ks, data_type, min_ks=min_ks)


def ks_list_from_df(df, max_ks, data_type, min_ks=0.005):
    """
    Same as ks_list_from_tsv, but takes an already-loaded DataFrame instead of reading a TSV file
    from disk. Useful when the caller already has the data in memory (e.g. reconstructed from the
    paralog Ks database) and doesn't want to pay for a redundant read.

    :param df: DataFrame with the same columns as a wgd Ks TSV file (at minimum: Family, Node, Ks,
               AlignmentCoverage, AlignmentIdentity, AlignmentLength)
    :param max_ks: maximum Ks value to be accepted for the analysis
    :param data_type: specifies the nature of the Ks values in df (either "paralogs", "anchor pairs", "reciprocally retained" or "orthologs")
    :param min_ks: minimum Ks value to be accepted for the analysis
    :return: ks_list_filtered, list of either whole-paranome, anchors pairs, reciprocally retained or ortholog Ks values
    :return: weight_list_filtered, list of weights for Ks values (not returned for orthologs Ks)
    """
    # Re-calculate the weights when excluding the outliers (i.e. Ks > max_ks)
    # Return the Ks data with the updated WeightOutliersExcluded column
    filtered_tsv_updated_weights = filter_compute_weights(df, min_ks, max_ks)

    ks_list_filtered = filtered_tsv_updated_weights["Ks"].to_list()
    weight_list_filtered = filtered_tsv_updated_weights["WeightOutliersExcluded"].to_list()

    if data_type == "paralogs" or data_type == "anchor pairs" or data_type == "reciprocally retained":
        return ks_list_filtered, weight_list_filtered
    if data_type == "orthologs":
        return ks_list_filtered
