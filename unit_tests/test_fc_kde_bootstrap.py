"""
Unit tests for fc_kde_bootstrap.estimate_peak()'s handling of ortholog pairs that have too few
(or numerically degenerate) Ks values below max_ks_ortho to support a bootstrap peak estimate.

Reproduces the real-world crash this guards against: a highly divergent species pair whose orthologs
are nearly all Ks-saturated (>10Ks), leaving only a handful of points below it - too few
for scipy.stats.gaussian_kde to build a non-singular covariance matrix. Before this guard, that
raised an uncaught numpy.linalg.LinAlgError that crashed the whole orthologs-analysis/wgdOrthologs
task; now it's treated as an expected outcome (skip the peak, log a warning).

Run with: python -m pytest unit_tests/test_fc_kde_bootstrap.py
"""
import os
import sys
import tempfile
import unittest

import numpy as np

# Makes the repo root importable as "ksrates.xxx", since this file lives one level down in unit_tests/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ksrates.fc_kde_bootstrap import estimate_peak, MIN_KS_VALUES_FOR_PEAK_ESTIMATE


def _write_ks_tsv(path, species1, species2, ks_values):
    """
    Writes a minimal wgd-style ortholog .ks.tsv fixture with just the columns
    fc_extract_ks_list.ks_list_from_tsv()/filter_compute_weights() actually read: Family, Node,
    Ks, AlignmentCoverage, AlignmentIdentity, AlignmentLength. AlignmentLength must clear
    filter_compute_weights()'s default 300bp minimum or the row gets silently filtered out.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write("\tAlignmentCoverage\tAlignmentIdentity\tAlignmentLength\tFamily\tKs\tNode\n")
        for i, ks in enumerate(ks_values):
            pair_id = f"{species1}-{i}__{species2}-{i}"
            handle.write(f"{pair_id}\t1.0\t0.5\t500\tGF_{i:06d}\t{ks}\t2.0\n")


class EstimatePeakTest(unittest.TestCase):

    def setUp(self):
        # Runs before every test method below. Gives each test its own throwaway directory and
        # drops us into it, since estimate_peak() reads/writes paths relative to the current
        # working directory (e.g. "ortholog_distributions/...", "peak_db.tsv")
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)  # delete the temp dir once the test ends
        self._original_cwd = os.getcwd()
        os.chdir(self.tmpdir.name)
        self.addCleanup(os.chdir, self._original_cwd)  # restore cwd even if the test fails
        self.peak_db_path = "peak_db.tsv"
        self.ks_list_db_path = "ks_list_db.tsv"

    def _estimate(self, ks_values, flag_not_in_peak_db=True, flag_not_in_ks_db=True):
        """Shared helper: writes a fake .ks.tsv fixture, then calls the real function under test."""
        ks_tsv = os.path.join("ortholog_distributions", "wgd_SP1_SP2", "SP1_SP2.ks.tsv")
        _write_ks_tsv(ks_tsv, "SP1", "SP2", ks_values)  # arrange: fake input data
        return estimate_peak(  # act: this is the actual function being tested
            "SP1", "SP2", "Species one", "Species two",
            max_ks_ortho=10, n_iter=20, x_lim_ortho=5, bin_width_ortho=0.1,
            ks_list_db_path=self.ks_list_db_path, db_path=self.peak_db_path,
            flag_not_in_peak_db=flag_not_in_peak_db, flag_not_in_ks_db=flag_not_in_ks_db,
        )

    def test_too_few_ks_values_skips_gracefully_without_raising(self):
        # Test probelamtic case: it has 3 Ks values, which will be an insufficient amount to proceed
        ks_values = [7.57, 7.80, 8.71]
        self.assertLess(len(ks_values), MIN_KS_VALUES_FOR_PEAK_ESTIMATE)

        failed = self._estimate(ks_values)  # this used to crash with LinAlgError before the fix

        # "failed" here just means "no peak computed" (a normal, expected outcome) - not a crash
        self.assertTrue(failed, "a too-few-values pair should report failure, not raise")
        peak_db_is_empty = (
            not os.path.exists(self.peak_db_path) or os.path.getsize(self.peak_db_path) == 0
        )
        self.assertTrue(peak_db_is_empty, "no peak should have been written for a failed estimate")

    def test_too_few_ks_values_still_records_the_raw_ks_list(self):
        # Even if the Ks peak is unavailable, the raw (already-extracted) Ks list still get stored
        self._estimate([7.57, 7.80, 8.71])

        with open(self.ks_list_db_path) as handle:
            content = handle.read()
        self.assertIn("Species one_Species two", content)  # the row key estimate_peak() writes

    def test_enough_ks_values_computes_a_real_peak(self):
        # Test normal case: many Ks values, spread widely enough to never produce a degenerate
        # bootstrap covariance matrix
        ks_values = list(np.linspace(0.5, 2.5, 20))

        failed = self._estimate(ks_values)

        self.assertFalse(failed)  # a real peak should be computed, so "failed" must be False
        with open(self.peak_db_path) as handle:
            content = handle.read()
        self.assertIn("Species one_Species two", content)

    def test_already_in_peak_and_ks_db_skips_recomputation_entirely(self):
        # flag_not_in_peak_db=False, flag_not_in_ks_db=False means "this pair is already fully
        # stored": nothing should be computed or written at all, regardless of how many Ks
        # values are available (even a too-small list here must NOT trigger a fresh computation).
        failed = self._estimate([1.0, 1.1, 1.2], flag_not_in_peak_db=False, flag_not_in_ks_db=False)

        self.assertFalse(failed)
        # Neither database file should have been touched - confirms no needless recomputation.
        self.assertFalse(os.path.exists(self.peak_db_path) and os.path.getsize(self.peak_db_path) > 0)
        self.assertFalse(os.path.exists(self.ks_list_db_path) and os.path.getsize(self.ks_list_db_path) > 0)


if __name__ == "__main__":
    unittest.main()
