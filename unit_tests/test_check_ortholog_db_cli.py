"""
Unit tests for the "ksrates check-ortholog-db" CLI command (ksrates_cli.py), which
setOrthologAnalysis (main.nf) uses to decide whether a species pair's wgd/peak computation can be
skipped entirely because it's already in the shared peak_database_path/ks_list_database_path
TSVs - mirrors "check-paralog-db"'s existing exit-code/stdout contract.

Run with: python -m pytest unit_tests/test_check_ortholog_db_cli.py
"""
import os
import sys
import tempfile
import unittest

from click.testing import CliRunner

# Makes the repo root importable, since this file lives one level down in unit_tests/.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ksrates_cli import cli


# Minimal fake ksrates config file, just enough for Configuration() to parse without errors.
# {peak_db}/{ks_list_db} get filled in per-test with paths inside that test's own temp dir.
CONFIG_TEMPLATE = """\
[SPECIES]
focal_species = sp1
newick_tree = ((sp1, sp2), sp3);
latin_names = sp1: Species one, sp2: Species two, sp3: Species three
peak_database_path = {peak_db}
ks_list_database_path = {ks_list_db}
"""


class CheckOrthologDbTest(unittest.TestCase):

    def setUp(self):
        # Runs before every test method below: gives each test its own throwaway temp dir and a
        # ready-to-use fake config file pointing at (not-yet-created) peak/Ks-list TSV paths
        # inside it. Individual tests then create those TSV files themselves (or don't, to
        # simulate "not in the database yet").
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.peak_db_path = os.path.join(self.tmpdir.name, "peak_db.tsv")
        self.ks_list_db_path = os.path.join(self.tmpdir.name, "ks_list_db.tsv")
        self.config_path = os.path.join(self.tmpdir.name, "config.txt")
        with open(self.config_path, "w") as handle:
            handle.write(CONFIG_TEMPLATE.format(
                peak_db=self.peak_db_path, ks_list_db=self.ks_list_db_path,
            ))
        self.runner = CliRunner()

    def _invoke(self, species1, species2):
        """Shared helper: runs "ksrates check-ortholog-db config.txt species1 species2"."""
        return self.runner.invoke(
            cli, ["check-ortholog-db", self.config_path, species1, species2],
            catch_exceptions=False,
        )

    def test_pair_missing_from_both_databases_reports_no(self):
        # Neither TSV exists yet at all - matches a brand new setup.
        result = self._invoke("sp1", "sp2")

        self.assertEqual(result.exit_code, 1)  # exit 1 = "not present" (main.nf treats this as "needs computing")
        self.assertEqual(result.output.strip(), "no")

    def test_pair_present_in_both_databases_reports_yes(self):
        # Arrange: manually write both TSVs with a row for the sp1/sp2 pair, as if a previous
        # run had already computed it.
        with open(self.peak_db_path, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tMode\tMode_SD\n")
            handle.write("Species one_Species two\tSpecies one\tSpecies two\t1.23\t0.1\n")
        with open(self.ks_list_db_path, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tKs_Values\n")
            handle.write("Species one_Species two\tSpecies one\tSpecies two\t[1.1, 1.2, 1.3]\n")

        result = self._invoke("sp1", "sp2")

        self.assertEqual(result.exit_code, 0)  # exit 0 = "present" (safe to skip recomputation)
        self.assertEqual(result.output.strip(), "yes")

    def test_pair_present_in_peak_db_only_reports_no(self):
        # Test problematic case: a Ks list can be stored even though peak estimation failed
        # and was skipped - this must not be treated as "fully present". Only the peak TSV is
        # written here; the Ks-list TSV is deliberately left missing.
        with open(self.peak_db_path, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tMode\tMode_SD\n")
            handle.write("Species one_Species two\tSpecies one\tSpecies two\t1.23\t0.1\n")
        # ks_list_db_path intentionally left missing.

        result = self._invoke("sp1", "sp2")

        # Must still say "no" - BOTH files need the pair, having just one isn't enough.
        self.assertEqual(result.exit_code, 1)
        self.assertEqual(result.output.strip(), "no")

    def test_pair_order_is_order_independent(self):
        # Same fixture as the "present" test above - only the argument order changes below.
        with open(self.peak_db_path, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tMode\tMode_SD\n")
            handle.write("Species one_Species two\tSpecies one\tSpecies two\t1.23\t0.1\n")
        with open(self.ks_list_db_path, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tKs_Values\n")
            handle.write("Species one_Species two\tSpecies one\tSpecies two\t[1.1, 1.2, 1.3]\n")

        # sp2/sp1 instead of sp1/sp2 on the command line - same pair, reversed argument order.
        result = self._invoke("sp2", "sp1")

        self.assertEqual(result.exit_code, 0)
        self.assertEqual(result.output.strip(), "yes")

    def test_different_pair_not_in_database_reports_no(self):
        # Same sp1/sp2 rows as the "present" test above, but this time we ask about sp1/sp3 -
        # a pair that was never stored. Guards against a too-loose match (e.g. "any row
        # mentioning sp1 counts") rather than an exact pair-key match.
        with open(self.peak_db_path, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tMode\tMode_SD\n")
            handle.write("Species one_Species two\tSpecies one\tSpecies two\t1.23\t0.1\n")
        with open(self.ks_list_db_path, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tKs_Values\n")
            handle.write("Species one_Species two\tSpecies one\tSpecies two\t[1.1, 1.2, 1.3]\n")

        result = self._invoke("sp1", "sp3")

        self.assertEqual(result.exit_code, 1)
        self.assertEqual(result.output.strip(), "no")


if __name__ == "__main__":
    unittest.main()
