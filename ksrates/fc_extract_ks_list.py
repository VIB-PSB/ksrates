from pandas import read_csv
from numpy import float64


def ks_list_from_tsv(tsv_file, max_ks, data_type):
    """
    Extracts paralog Ks (whole-paranome, anchor pairs or reciprocally retained) or ortholog Ks from a wgd Ks TSV file.
    If dealing with paralog Ks values, it returns also a list of weights extracted from "WeightOutliersExcluded" column in the tsv_file.

    :param tsv_file: wgd output file containing either paralog or ortholog Ks values (suffix formats: ".ks.tsv", "ks_anchors.tsv", "ks_recret_topX.tsv")
    :param max_ks: maximum Ks value to be accepted for the analysis
    :param data_type: specifies the nature of the Ks values in the tsv file (either "paralogs", "anchor pairs", "reciprocally retained" or "orthologs")
    :return: ks_list_filtered, list of either whole-paranome, anchors pairs, reciprocally retained or ortholog Ks values
    :return: weight_list_filtered, list of weights for Ks values (not returned for orthologs Ks)
    """
    with open(tsv_file, "r") as tsv_file:
        tsv = read_csv(tsv_file, sep="\t")

    filtered_tsv = tsv.loc[(tsv["Ks"].dtypes == float64) & (tsv["Ks"] <= max_ks)]

    ks_list_filtered = filtered_tsv["Ks"].to_list()
    weight_list_filtered = filtered_tsv["WeightOutliersExcluded"].to_list()

    if data_type == "paralogs" or data_type == "anchor pairs" or data_type == "reciprocally retained":
        return ks_list_filtered, weight_list_filtered
    if data_type == "orthologs":
        return ks_list_filtered
