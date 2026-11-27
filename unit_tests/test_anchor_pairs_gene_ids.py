"""
Unit tests for the handling of the gene IDs in the Ks tables of the anchor pairs. The gene IDs of some genomes are made only
of digits, and an ID can also start with zeros: they have to be handled as text everywhere, like the names of the pairs in the
Ks table, which are the two IDs sorted as text and joined with "__". Read by pandas as numbers, the IDs would crash the join,
be sorted as numbers (so the pair names would not match those of the Ks table) and lose their leading zeros.

Covers fc_wgd._write_anchor_pairs_ks (matches the anchor pairs of i-ADHoRe to their Ks values),
fc_cluster_anchors.parse_ks_anchors_from_df (the lookup of the anchor Ks values) and
fc_ks_database_paralog._extract_columns_from_tsv (the Ks tables that go into the database).

Run with: python -m pytest unit_tests/test_anchor_pairs_gene_ids.py
"""
import os
import sys
import tempfile
import unittest

import pandas

# Makes the repo root importable, since this file lives one level down in unit_tests/.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import ksrates.fc_wgd as fcWgd
import ksrates.fc_cluster_anchors as fcCluster
import ksrates.fc_ks_database_paralog as fcKsDbPara

ANCHORPOINTS_HEADER = "id\tmultiplicon\tbasecluster\tgene_x\tgene_y\tcoord_x\tcoord_y\tis_real_anchorpoint\n"


def pair_name(gene_1, gene_2):
    """Name of a pair in the Ks table: the two IDs sorted as text."""
    return "__".join(sorted([gene_1, gene_2]))


def ks_table(pairs):
    """Ks table as the paranome step writes it, with the IDs as text; pairs is a list of (gene_1, gene_2)."""
    return pandas.DataFrame({
        "AlignmentCoverage": 1.0, "AlignmentIdentity": 1.0, "AlignmentLength": 600.0,
        "Family": "GF_000001", "Ks": [0.5 + 0.1 * i for i in range(len(pairs))],
        "Node": 2.0, "Paralog1": [p[0] for p in pairs], "Paralog2": [p[1] for p in pairs],
        "WeightOutliersIncluded": 1.0,
    }, index=[pair_name(*p) for p in pairs])


class WriteAnchorPairsKsTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.anchorpoints = os.path.join(self.tmp.name, "anchorpoints.txt")
        self.ks_file = os.path.join(self.tmp.name, "ks.tsv")
        self.out_file = os.path.join(self.tmp.name, "ks_anchors.tsv")

    def tearDown(self):
        self.tmp.cleanup()

    def write_inputs(self, anchor_pairs, ks_pairs):
        """anchor_pairs: (gene_x, gene_y) as written by i-ADHoRe; ks_pairs: (gene_1, gene_2) of the rows of the Ks table."""
        with open(self.anchorpoints, "w") as f:
            f.write(ANCHORPOINTS_HEADER)
            for i, (gene_x, gene_y) in enumerate(anchor_pairs, start=1):
                f.write(f"{i}\t1\t1\t{gene_x}\t{gene_y}\t{i}\t{i}\t-1\n")
        ks_table(ks_pairs).to_csv(self.ks_file, sep="\t")

    def output(self):
        return pandas.read_table(self.out_file, index_col=0, dtype={"Paralog1": str, "Paralog2": str})

    def run_function(self):
        fcWgd._write_anchor_pairs_ks(self.anchorpoints, self.ks_file, self.out_file, 0.05)

    def test_gene_ids_with_letters(self):
        self.write_inputs([("AT1G01020", "AT1G01010")], [("AT1G01010", "AT1G01020"), ("AT1G09999", "AT1G09998")])
        self.run_function()
        self.assertEqual(sorted(self.output().index), ["AT1G01010__AT1G01020"])

    def test_gene_ids_made_only_of_digits(self):
        # The pair names of the Ks table are sorted as text ("105005" before "37570"), not as numbers
        self.write_inputs([(37570, 105005), (404282, 404126)], [("37570", "105005"), ("404126", "404282"), ("1", "2")])
        self.run_function()
        self.assertEqual(sorted(self.output().index), ["105005__37570", "404126__404282"])

    def test_leading_zeros_are_kept_in_the_pair_names(self):
        self.write_inputs([("00123", "00045")], [("00123", "00045")])
        self.run_function()
        self.assertEqual(sorted(self.output().index), ["00045__00123"])

    def test_ids_in_the_output_table_keep_their_text(self):
        self.write_inputs([("00123", "00045"), ("7", "100")], [("00123", "00045"), ("7", "100")])
        self.run_function()
        table = self.output()
        self.assertEqual(sorted(zip(table["Paralog1"], table["Paralog2"])), [("00123", "00045"), ("7", "100")])

    def test_anchor_pairs_without_a_ks_value_are_left_out(self):
        self.write_inputs([(1, 2), (3, 4)], [("1", "2")])
        self.run_function()
        self.assertEqual(sorted(self.output().index), ["1__2"])

    def test_repeated_anchor_pair_is_written_once(self):
        self.write_inputs([(5, 6), (6, 5)], [("5", "6")])
        self.run_function()
        self.assertEqual(sorted(self.output().index), ["5__6"])


class ParseKsAnchorsFromDfTest(unittest.TestCase):

    def test_text_ids(self):
        df = pandas.DataFrame({"Paralog1": ["b", "a"], "Paralog2": ["a", "c"], "Ks": [0.5, 1.5]})
        self.assertEqual(fcCluster.parse_ks_anchors_from_df(df), {("a", "b"): "0.5", ("a", "c"): "1.5"})

    def test_ids_stored_as_numbers_are_looked_up_as_text(self):
        # The anchor pairs of anchorpoints.txt are tuples of text sorted as text
        df = pandas.DataFrame({"Paralog1": [37570], "Paralog2": [105005], "Ks": [0.5]})
        self.assertEqual(fcCluster.parse_ks_anchors_from_df(df), {("105005", "37570"): "0.5"})

    def test_leading_zeros_are_kept(self):
        df = pandas.DataFrame({"Paralog1": ["00123"], "Paralog2": ["00045"], "Ks": [0.5]})
        self.assertEqual(fcCluster.parse_ks_anchors_from_df(df), {("00045", "00123"): "0.5"})


class ExtractColumnsFromTsvTest(unittest.TestCase):

    def test_ids_keep_their_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ks.tsv")
            ks_table([("00123", "00045"), ("37570", "105005")]).to_csv(path, sep="\t")
            columns = fcKsDbPara._extract_columns_from_tsv(path)
        self.assertEqual(columns["Paralog1"], ["00123", "37570"])
        self.assertEqual(columns["Paralog2"], ["00045", "105005"])
        self.assertEqual(len(columns["Ks"]), 2)


if __name__ == "__main__":
    unittest.main()
