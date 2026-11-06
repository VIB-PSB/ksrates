"""
Unit tests for the removal of the large intermediate files of the reciprocal retention pipeline
(fc_rec_ret_orthomcl.remove_orthomcl_intermediate_files and the expert configuration field
remove_reciprocal_retention_orthomcl_files).

Run with: python -m pytest unit_tests/test_rec_ret_cleanup.py
"""
import os
import sys
import tempfile
import unittest
from configparser import ConfigParser

# Makes the repo root importable, since this file lives one level down in unit_tests/.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import ksrates.fc_configfile as fcConf
import ksrates.fc_rec_ret_orthomcl as fcRecRet

SPECIES = "Sp1"


def write(path, text="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


class RemoveOrthomclIntermediateFilesTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = os.path.join(self.tmp.name, "reciprocal_retention", "orthomcl")
        self.focal_aa = os.path.join(self.dir, "Sp1-translated.faa")
        self.orthomcl_out = os.path.join(self.dir, f"{SPECIES}_original_37.orthomcl.out")
        # Files to keep
        write(self.orthomcl_out, "ORTHOMCL0(2 genes, 2 taxa):\tg1(sp) g2(sp)\n")
        write(os.path.join(self.dir, f"{SPECIES}_matched_top_2000_gfs_max_15_extra.tsv"))
        write(os.path.join(self.dir, f"{SPECIES}_matched_top_4000_gfs_max_15_extra.tsv"))
        write(os.path.join(self.dir, f"{SPECIES}_rec_ret_top_2000_ids.mcl.tsv"))
        write(os.path.join(self.dir, f"{SPECIES}_rec_ret_top_4000_ids.mcl.tsv"))
        # Files to remove
        for name in (f"{SPECIES}_original_37.dmd.tsv", f"{SPECIES}_original_37.dmnd", f"{SPECIES}_original_37.faa",
                     "all.gg", "Sp1-translated.faa"):
            write(os.path.join(self.dir, name), "y" * 1000)
        write(os.path.join(self.dir, "Sep_1", "parameter.log"), "MCL Inflation  1.5\n")
        write(os.path.join(self.dir, "Sep_1", "all_orthomcl.out"), "y" * 1000)
        write(os.path.join(self.dir, "Sep_1", "tmp", "all.bpo"), "y" * 1000)

    def tearDown(self):
        self.tmp.cleanup()

    def clean(self, **kwargs):
        return fcRecRet.remove_orthomcl_intermediate_files(SPECIES, self.dir, self.focal_aa, self.orthomcl_out, **kwargs)

    def names(self):
        return sorted(os.listdir(self.dir))

    def test_keeps_only_the_results_and_the_parameter_log(self):
        self.clean()
        self.assertEqual(self.names(), sorted([
            f"{SPECIES}_matched_top_2000_gfs_max_15_extra.tsv",
            f"{SPECIES}_matched_top_4000_gfs_max_15_extra.tsv",
            f"{SPECIES}_orthomcl_parameter.log",
            f"{SPECIES}_original_37.orthomcl.out",
            f"{SPECIES}_rec_ret_top_2000_ids.mcl.tsv",
            f"{SPECIES}_rec_ret_top_4000_ids.mcl.tsv",
        ]))

    def test_parameter_log_keeps_its_content(self):
        self.clean()
        with open(os.path.join(self.dir, f"{SPECIES}_orthomcl_parameter.log")) as f:
            self.assertEqual(f.read(), "MCL Inflation  1.5\n")

    def test_returns_the_freed_bytes(self):
        # 5 files in the directory and 2 files in the working directory, 1000 bytes each; parameter.log is moved, not counted
        self.assertEqual(self.clean(), 7000)

    def test_nothing_is_removed_without_the_orthogroups_file(self):
        os.remove(self.orthomcl_out)
        before = self.names()
        self.assertEqual(self.clean(), 0)
        self.assertEqual(self.names(), before)

    def test_empty_orthogroups_file_counts_as_missing(self):
        write(self.orthomcl_out, "")
        before = self.names()
        self.assertEqual(self.clean(), 0)
        self.assertEqual(self.names(), before)

    def test_protected_path_is_kept(self):
        protected = os.path.join(self.dir, f"{SPECIES}_original_37.dmd.tsv")
        self.clean(protected_paths=[protected, ""])
        self.assertIn(f"{SPECIES}_original_37.dmd.tsv", self.names())
        self.assertNotIn("all.gg", self.names())

    def test_logs_what_was_done(self):
        with self.assertLogs(level="INFO") as first:
            self.clean()
        self.assertTrue(any("Removed the OrthoMCL intermediate files" in line for line in first.output))
        with self.assertLogs(level="INFO") as second:
            self.clean()
        self.assertTrue(any("No OrthoMCL intermediate files to remove" in line for line in second.output))
        self.assertFalse(any("Removed the OrthoMCL" in line for line in second.output))

    def test_second_call_is_harmless(self):
        self.clean()
        before = self.names()
        self.assertEqual(self.clean(), 0)
        self.assertEqual(self.names(), before)

    def test_existing_parameter_log_is_not_overwritten(self):
        write(os.path.join(self.dir, f"{SPECIES}_orthomcl_parameter.log"), "first run\n")
        self.clean()
        with open(os.path.join(self.dir, f"{SPECIES}_orthomcl_parameter.log")) as f:
            self.assertEqual(f.read(), "first run\n")
        self.assertIn(f"{SPECIES}_orthomcl_parameter_Sep_1.log", self.names())

    def test_other_directories_are_left_alone(self):
        write(os.path.join(self.dir, "notes", "keep.txt"))              # name is not a date
        write(os.path.join(self.dir, "Oct_2", "something.txt"))         # date name but not an OrthoMCL working directory
        self.clean()
        self.assertTrue(os.path.isfile(os.path.join(self.dir, "notes", "keep.txt")))
        self.assertTrue(os.path.isfile(os.path.join(self.dir, "Oct_2", "something.txt")))

    def test_several_working_directories(self):
        write(os.path.join(self.dir, "Sep_1_2", "parameter.log"), "second run\n")
        write(os.path.join(self.dir, "Sep_1_2", "tmp", "all.bpo"), "y" * 1000)
        self.clean()
        self.assertFalse(os.path.exists(os.path.join(self.dir, "Sep_1")))
        self.assertFalse(os.path.exists(os.path.join(self.dir, "Sep_1_2")))
        self.assertIn(f"{SPECIES}_orthomcl_parameter.log", self.names())
        self.assertIn(f"{SPECIES}_orthomcl_parameter_Sep_1_2.log", self.names())

    def test_focal_fasta_outside_the_orthomcl_directory_is_not_removed(self):
        outside = os.path.join(self.tmp.name, "input.faa")
        write(outside)
        fcRecRet.remove_orthomcl_intermediate_files(SPECIES, self.dir, outside, self.orthomcl_out)
        self.assertTrue(os.path.isfile(outside))


class RemoveOrthomclFilesConfigTest(unittest.TestCase):

    @staticmethod
    def config(expert_text):
        """A Configuration without reading any file, with only the expert configuration set."""
        config = fcConf.Configuration.__new__(fcConf.Configuration)
        if expert_text is None:
            config.expert_config = None
        else:
            config.expert_config = ConfigParser()
            config.expert_config.read_string(expert_text)
        return config

    def test_default_is_to_keep_the_files(self):
        self.assertFalse(self.config("[EXPERT PARAMETERS]\n").get_remove_reciprocal_retention_orthomcl_files(True))
        self.assertFalse(self.config(None).get_remove_reciprocal_retention_orthomcl_files(True))

    def test_yes_and_no(self):
        yes = self.config("[EXPERT PARAMETERS]\nremove_reciprocal_retention_orthomcl_files = Yes\n")
        no = self.config("[EXPERT PARAMETERS]\nremove_reciprocal_retention_orthomcl_files = no\n")
        self.assertTrue(yes.get_remove_reciprocal_retention_orthomcl_files(True))
        self.assertFalse(no.get_remove_reciprocal_retention_orthomcl_files(True))

    def test_unrecognised_value_keeps_the_files(self):
        config = self.config("[EXPERT PARAMETERS]\nremove_reciprocal_retention_orthomcl_files = maybe\n")
        self.assertFalse(config.get_remove_reciprocal_retention_orthomcl_files(True))

    def test_not_used_without_reciprocal_retention(self):
        config = self.config("[EXPERT PARAMETERS]\nremove_reciprocal_retention_orthomcl_files = yes\n")
        self.assertFalse(config.get_remove_reciprocal_retention_orthomcl_files(False))


if __name__ == "__main__":
    unittest.main()
