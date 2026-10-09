"""
Unit tests for the ortholog part of the unified populate-ks-db / inspect-ks-db commands: legacy ortholog TSVs
-> database -> TSV dump must round-trip. The sqld server is replaced by a local SQLite file (see fake_sqld.py).

Run with: python -m pytest unit_tests/test_populate_inspect_ks_db.py
"""
import glob
import os
import sys
import tempfile
import unittest

import pandas
from click.testing import CliRunner

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

import ksrates.fc_ks_database_utilities as fcKsDbUtil
import ksrates.fc_ks_database_ortholog as fcKsDbOrtho
from ksrates_cli import cli
from fake_sqld import install_fake_sqld


class PopulateInspectOrthologTest(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.db, self.unreachable_db = install_fake_sqld(self, self.tmpdir.name)
        self.peak_tsv = os.path.join(self.tmpdir.name, "peak.tsv")
        self.ks_tsv = os.path.join(self.tmpdir.name, "ks.tsv")
        with open(self.peak_tsv, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tMode\tMode_SD\n")
            handle.write("Aaa bbb_Ccc ddd\tAaa bbb\tCcc ddd\t1.5\t0.1\n")
            handle.write("Aaa bbb_Eee fff\tAaa bbb\tEee fff\t2.5\t0.2\n")
        with open(self.ks_tsv, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tKs_Values\n")
            handle.write("Aaa bbb_Ccc ddd\tAaa bbb\tCcc ddd\t[0.1, 0.2, 0.3]\n")
        self.runner = CliRunner()

    def _populate(self, *extra, database=None):
        return self.runner.invoke(cli, ["populate-ks-db", "--orthologs", self.peak_tsv, self.ks_tsv,
                                        "--database", database or self.db, *extra], catch_exceptions=False)

    def test_populate_copies_values_forward(self):
        result = self._populate()
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(fcKsDbOrtho.read_peak(self.db, "Aaa bbb", "Ccc ddd"), (1.5, 0.1))
        self.assertEqual(fcKsDbOrtho.read_ks_list(self.db, "Ccc ddd", "Aaa bbb"), [0.1, 0.2, 0.3])
        self.assertEqual(fcKsDbOrtho.read_peak(self.db, "Aaa bbb", "Eee fff"), (2.5, 0.2))
        self.assertIsNone(fcKsDbOrtho.read_ks_list(self.db, "Aaa bbb", "Eee fff"))  # peak-only pair stays peak-only

    def test_populate_never_mixes_halves_from_different_sources(self):
        fcKsDbUtil.initialize_database(self.db)
        # database has only a peak for a pair that the TSVs have complete: replaced as a consistent whole
        fcKsDbOrtho.write_to_ortholog_db(self.db, "Aaa bbb", "Ccc ddd", peak=(9.9, 0.9))
        # database has a Ks list for a pair of which the TSVs only have the peak: not mixed, skipped
        fcKsDbOrtho.write_to_ortholog_db(self.db, "Aaa bbb", "Eee fff", ks_list=[7.0, 8.0])
        self._populate()
        self.assertEqual(fcKsDbOrtho.read_peak(self.db, "Aaa bbb", "Ccc ddd"), (1.5, 0.1))
        self.assertEqual(fcKsDbOrtho.read_ks_list(self.db, "Aaa bbb", "Ccc ddd"), [0.1, 0.2, 0.3])
        self.assertIsNone(fcKsDbOrtho.read_peak(self.db, "Aaa bbb", "Eee fff"))
        self.assertEqual(fcKsDbOrtho.read_ks_list(self.db, "Aaa bbb", "Eee fff"), [7.0, 8.0])

    def test_populate_leaves_complete_pairs_alone_unless_forced(self):
        fcKsDbUtil.initialize_database(self.db)
        fcKsDbOrtho.write_to_ortholog_db(self.db, "Aaa bbb", "Ccc ddd", peak=(9.9, 0.9), ks_list=[5.0])
        self._populate()
        self.assertEqual(fcKsDbOrtho.read_peak(self.db, "Aaa bbb", "Ccc ddd"), (9.9, 0.9))
        self.assertEqual(fcKsDbOrtho.read_ks_list(self.db, "Aaa bbb", "Ccc ddd"), [5.0])
        self._populate("--force")
        self.assertEqual(fcKsDbOrtho.read_peak(self.db, "Aaa bbb", "Ccc ddd"), (1.5, 0.1))
        self.assertEqual(fcKsDbOrtho.read_ks_list(self.db, "Aaa bbb", "Ccc ddd"), [0.1, 0.2, 0.3])

    def test_populate_requires_at_least_one_kind_of_data(self):
        result = self.runner.invoke(cli, ["populate-ks-db", "--database", self.db])
        self.assertNotEqual(result.exit_code, 0)

    def test_populate_with_unreachable_database_fails(self):
        result = self._populate(database=self.unreachable_db)
        self.assertNotEqual(result.exit_code, 0)

    def test_inspect_dump_roundtrips_to_legacy_tsv_format(self):
        self._populate()
        result = self.runner.invoke(cli, ["inspect-ks-db", self.db, "--type", "ortholog"], catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        peaks_dump = glob.glob(os.path.join(self.tmpdir.name, "ks_database_server_address_ortholog_peaks_*.tsv"))[0]
        ks_dump = glob.glob(os.path.join(self.tmpdir.name, "ks_database_server_address_ortholog_ks_lists_*.tsv"))[0]
        pandas.testing.assert_frame_equal(pandas.read_csv(peaks_dump, sep="\t", index_col=0),
                                          pandas.read_csv(self.peak_tsv, sep="\t", index_col=0))
        # Ks lists: same legacy format, readable by the same code that reads the legacy file
        original = pandas.read_csv(self.ks_tsv, sep="\t", index_col=0)
        dumped = pandas.read_csv(ks_dump, sep="\t", index_col=0)
        pandas.testing.assert_frame_equal(dumped, original)

    def test_inspect_filter_list_and_delete(self):
        self._populate()
        result = self.runner.invoke(cli, ["inspect-ks-db", self.db, "--type", "ortholog", "--list"], catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        listing = glob.glob(os.path.join(self.tmpdir.name, "*ortholog_pairs_list_*.tsv"))[0]
        with open(listing) as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0], "pair_key\tpeak\tks_list")
        self.assertIn("Aaa bbb_Eee fff\tTrue\tFalse", lines)
        # delete needs an explicit type and a filter, and asks for confirmation
        self.assertNotEqual(self.runner.invoke(cli, ["inspect-ks-db", self.db, "Eee", "--delete"]).exit_code, 0)
        self.assertNotEqual(self.runner.invoke(cli, ["inspect-ks-db", self.db, "--type", "ortholog", "--delete"]).exit_code, 0)
        result = self.runner.invoke(cli, ["inspect-ks-db", self.db, "Eee", "--type", "ortholog", "--delete"], input="y\n", catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(fcKsDbUtil.list_keys(self.db, fcKsDbOrtho.ORTHOLOG_TABLE), ["Aaa bbb_Ccc ddd"])


if __name__ == "__main__":
    unittest.main()
