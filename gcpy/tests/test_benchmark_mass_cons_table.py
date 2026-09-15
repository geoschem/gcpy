"""
Unit tests for gcpy/benchmark/modules/benchmark_mass_cons_table.py.

Covers the difference statistics reported beneath the monthly mass
table, where the End mass row's Abs Diff and % Diff columns were
computed from the start-of-run mass instead of the end-of-run mass.
"""
import numpy as np

from gcpy.benchmark.modules.benchmark_mass_cons_table import \
    compute_diff_statistics, compute_statistics


# Ref and Dev global PassiveTracer masses [Tg], as reported for
# gcc.14.8.0-rc.0 vs gchp.14.8.0-rc.0.  Ref drifts in its last few
# digits; Dev is constant to the precision np.sum can resolve.
REF_MASSES = np.array([
    17.656309337013, 17.656309337164, 17.656309337116, 17.656309336992,
    17.656309337007, 17.656309337124, 17.656309337151, 17.656309337125,
    17.656309337052, 17.656309337093, 17.656309337046, 17.656309337140,
])
DEV_MASSES = np.full(12, 17.656423649381)


def test_compute_statistics_start_and_end_are_distinct():
    stats = compute_statistics(REF_MASSES)
    assert stats["start_mass"] == REF_MASSES[0]
    assert stats["end_mass"] == REF_MASSES[-1]
    assert stats["start_mass"] != stats["end_mass"]


def test_end_mass_diff_uses_end_mass_not_start_mass():
    # Regression test: compute_diff_statistics computed its end_mass
    # entry from the "start_mass" key, so the End mass row of the
    # table duplicated the Start mass row's Abs Diff and % Diff.
    ref = compute_statistics(REF_MASSES)
    dev = compute_statistics(DEV_MASSES)
    diff = compute_diff_statistics(ref, dev)

    assert diff["end_mass__absdiff"] == (
        dev["end_mass"] - ref["end_mass"]
    )
    assert diff["end_mass__absdiff"] != diff["start_mass__absdiff"], (
        "End mass diff is being taken from the start-of-run mass"
    )


def test_start_mass_diff_still_uses_start_mass():
    ref = compute_statistics(REF_MASSES)
    dev = compute_statistics(DEV_MASSES)
    diff = compute_diff_statistics(ref, dev)

    assert diff["start_mass__absdiff"] == (
        dev["start_mass"] - ref["start_mass"]
    )


def test_subgram_conservation_survives_the_abs_diff_format():
    # The Abs diff [g] rows used to run through np.int64, so a spread
    # under one gram printed as an exact 0 and read as perfect
    # conservation.  GCHP's real spread is a fraction of a gram.
    masses = np.array([17.656423649381, 17.656423649381 + 3.6e-13])
    stats = compute_statistics(masses)
    grams = stats["minmax_absdiff_g"]

    assert 0.0 < grams < 1.0, "expected a sub-gram spread"
    assert np.int64(grams) == 0, "the old format really did print 0"
    assert float(f"{grams:.4e}") > 0.0, "the new format must not vanish"
