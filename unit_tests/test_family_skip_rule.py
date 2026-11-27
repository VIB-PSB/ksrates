"""
Unit tests for the rule that skips the codeml analysis of paralog gene families whose alignment has too few
gap-free codon columns (count_gapfree_codons, FamilySkipRule, family_skip_reason, and the expert configuration
fields min_gapfree_columns, big_family_size and min_gapfree_columns_big_family).

Run with: python -m pytest unit_tests/test_family_skip_rule.py
"""
import os
import sys
import unittest
from configparser import ConfigParser

# Makes the repo root importable, since this file lives one level down in unit_tests/.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import ksrates.fc_configfile as fcConf
from wgd_ksrates.alignment import count_gapfree_codons
from wgd_ksrates.ks_distribution import FamilySkipRule, family_skip_reason, noks_skipped_summary


class CountGapfreeCodonsTest(unittest.TestCase):

    def test_all_columns_without_gaps(self):
        self.assertEqual(count_gapfree_codons({"a": "ATGAAACCC", "b": "ATGAAGCCC"}), 3)

    def test_a_gap_in_one_sequence_removes_the_column(self):
        self.assertEqual(count_gapfree_codons({"a": "ATGAAA---CCC", "b": "ATGAAATTTCCC"}), 3)

    def test_ambiguity_characters_remove_the_column(self):
        self.assertEqual(count_gapfree_codons({"a": "ATGNNNCCC", "b": "ATGAAACCC"}), 2)

    def test_lowercase_bases_count(self):
        self.assertEqual(count_gapfree_codons({"a": "atgaaa", "b": "ATGAAA"}), 2)

    def test_incomplete_last_codon_is_ignored(self):
        self.assertEqual(count_gapfree_codons({"a": "ATGA", "b": "ATGA"}), 1)

    def test_no_column_without_gaps(self):
        self.assertEqual(count_gapfree_codons({"a": "ATG---", "b": "---AAA"}), 0)

    def test_empty_alignment(self):
        self.assertEqual(count_gapfree_codons({}), 0)


class FamilySkipReasonTest(unittest.TestCase):

    def test_no_rule_skips_nothing(self):
        self.assertIsNone(family_skip_reason(200, 0, None))

    def test_rule_with_zero_limits_is_off(self):
        self.assertIsNone(family_skip_reason(200, 0, FamilySkipRule(0, 100, 0)))

    def test_min_columns_skips_families_below_the_limit(self):
        rule = FamilySkipRule(2, 100, 0)
        self.assertIsNotNone(family_skip_reason(10, 0, rule))
        self.assertIsNotNone(family_skip_reason(10, 1, rule))
        self.assertIsNone(family_skip_reason(10, 2, rule))
        self.assertIsNone(family_skip_reason(500, 2, rule))

    def test_big_family_limit_only_applies_to_big_families(self):
        rule = FamilySkipRule(0, 100, 4)
        self.assertIsNone(family_skip_reason(99, 3, rule))
        self.assertIsNotNone(family_skip_reason(100, 3, rule))
        self.assertIsNone(family_skip_reason(100, 4, rule))

    def test_both_limits_together(self):
        rule = FamilySkipRule(2, 100, 4)
        self.assertIn("fewer than 2", family_skip_reason(20, 1, rule))
        self.assertIn("family of 150 genes", family_skip_reason(150, 3, rule))
        self.assertIsNone(family_skip_reason(20, 3, rule))
        self.assertIsNone(family_skip_reason(150, 4, rule))


class NoKsSummaryTest(unittest.TestCase):

    def test_no_skipped_family(self):
        texts = {"GF_000001": "No codeml results", "GF_000002": "No clustering/weighting results"}
        self.assertEqual(noks_skipped_summary(texts), [])

    def test_skipped_families_are_counted_and_listed(self):
        texts = {
            "GF_000001": "No codeml results",
            "GF_000002": "Skipped before codeml: 0 gap-free codon columns, fewer than 2 [genes=137 columns=0]",
            "GF_000003": "Skipped before codeml: 0 gap-free codon columns, fewer than 2 [genes=81 columns=0]",
            "GF_000004": "Skipped before codeml: 1 gap-free codon columns, fewer than 2 [genes=57 columns=1]",
            "GF_000005": "Skipped before codeml: 3 gap-free codon columns, fewer than 4, in a family of 191 genes "
                         "(at least 100) [genes=191 columns=3]",
        }
        lines = noks_skipped_summary(texts)
        self.assertIn("4 gene families were skipped before codeml", lines[0])
        self.assertIn(" - 2 with no gap-free column at all", lines)
        self.assertIn(" - 2 with some columns, but too few:", lines)
        self.assertIn("     GF_000004: 57 genes, 1 gap-free codon columns", lines)
        self.assertIn("     GF_000005: 191 genes, 3 gap-free codon columns", lines)
        self.assertFalse(any("GF_000001" in line for line in lines))


class ExpertFieldsTest(unittest.TestCase):

    @staticmethod
    def configuration(expert_text=None):
        """A Configuration without reading any file, with only the expert configuration set."""
        config = fcConf.Configuration.__new__(fcConf.Configuration)
        if expert_text is None:
            config.expert_config = None
        else:
            config.expert_config = ConfigParser()
            config.expert_config.read_string(expert_text)
        return config

    def test_defaults_without_expert_configuration(self):
        config = self.configuration(None)
        self.assertEqual(config.get_min_gapfree_columns(), 0)
        self.assertEqual(config.get_big_family_size(), 100)
        self.assertEqual(config.get_min_gapfree_columns_big_family(), 0)

    def test_defaults_when_the_fields_are_missing(self):
        config = self.configuration("[EXPERT PARAMETERS]\nmax_gene_family_size = 200\n")
        self.assertEqual(config.get_min_gapfree_columns(), 0)
        self.assertEqual(config.get_big_family_size(), 100)
        self.assertEqual(config.get_min_gapfree_columns_big_family(), 0)

    def test_values_are_read(self):
        config = self.configuration(
            "[EXPERT PARAMETERS]\nmin_gapfree_columns = 2\nbig_family_size = 80\nmin_gapfree_columns_big_family = 4\n")
        self.assertEqual(config.get_min_gapfree_columns(), 2)
        self.assertEqual(config.get_big_family_size(), 80)
        self.assertEqual(config.get_min_gapfree_columns_big_family(), 4)

    def test_invalid_values_give_the_default(self):
        config = self.configuration(
            "[EXPERT PARAMETERS]\nmin_gapfree_columns = -1\nbig_family_size = many\nmin_gapfree_columns_big_family = 2.5\n")
        self.assertEqual(config.get_min_gapfree_columns(), 0)
        self.assertEqual(config.get_big_family_size(), 100)
        self.assertEqual(config.get_min_gapfree_columns_big_family(), 0)


if __name__ == "__main__":
    unittest.main()
