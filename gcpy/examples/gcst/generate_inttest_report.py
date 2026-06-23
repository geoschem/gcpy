#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_inttest_report.py
==========================
Parse GEOS-Chem Classic or GCHP integration-test execution and diff logs
and produce a concise plain-text report summarising:

  1. Which execution tests passed / failed / are still pending.
  2. Which tests produce output that differs from a reference version.

The model type (GCClassic vs GCHP) is auto-detected from the execution log.

Any annotation of *why* tests differ (e.g. "due to GFAS changes") or
flagging of known issues is intentionally left to the user; this script
reports facts only.

Usage
-----
    python generate_inttest_report.py \\
        --exec-log  <execution_log.txt>  \\
        --diff-log  <diff_log.txt>        \\
        --ref-label <gcc.14.8.0-alpha.11> \\
        [--output   <report.txt>]

Arguments
---------
--exec-log    Path to the execution test log produced by the integration
              test suite (GCClassic or GCHP format both accepted).
--diff-log    Path to the diff log produced by diffTest.sh.
--ref-label   Human-readable reference-version label that appears in
              diff-log file paths (e.g. gcc.14.8.0-alpha.11).
--output      Optional output file.  Default: stdout.

Notes
-----
* A test is flagged as "differing" if any non-zero number of differences
  is reported for either its OutputDir or Restarts, or if a raw ``diff``
  invocation produced output within that test's section of the diff log.
* A "perfect" run (all execution tests passed AND all zero-diff) produces
  only the banner and a single confirmation line.
* Only standard-library modules are used (argparse, pathlib, re, sys),
  so this script is compatible with any Python >= 3.9 environment,
  including the GCPy conda environment.
"""

import argparse
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Formatting constants
# ---------------------------------------------------------------------------

_SEP_HEAVY = '=' * 78
_SEP_LIGHT = '-' * 78
_BANNER_LINES = [
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
    '%%%  All execution tests passed!  %%%',
    '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
]


# ---------------------------------------------------------------------------
# Shared helper
# ---------------------------------------------------------------------------

def _find_int(pattern: str, text: str, default: int = 0) -> int:
    """Return the first integer captured by *pattern* in *text*, or *default*."""
    match = re.search(pattern, text)
    return int(match.group(1)) if match else default


# ---------------------------------------------------------------------------
# Execution-log parser
# ---------------------------------------------------------------------------

def parse_exec_log(log_path: str) -> dict:
    """
    Parse a GCClassic or GCHP execution test log.

    Parameters
    ----------
    log_path : str
        Path to the execution test log file.

    Returns
    -------
    dict with the following keys:

    ``model``
        ``'GCClassic'`` or ``'GCHP'``.
    ``header_block``
        ``list[str]`` – verbatim lines from the opening ``====...`` line
        through the closing ``====...`` line (inclusive).  Reproduced as-is
        in the report so that commit hashes, SLURM job IDs, etc. are
        preserved exactly.
    ``test_lines_raw``
        ``list[str]`` – verbatim lines from the *Execution tests:*
        subsection (blank lines included).  Reproduced verbatim so the
        dot-padded alignment of the original log is kept.
    ``tests``
        ``list[(name: str, result: str)]``.  *result* is ``'PASS'``,
        ``'FAIL'``, or ``'OTHER'`` (running / not-yet-completed).
    ``n_total``, ``n_passed``, ``n_failed``, ``n_pending``
        Counts parsed directly from the log.

    Raises
    ------
    ValueError
        If the model type cannot be detected or required section headers
        are absent.
    OSError
        If *log_path* cannot be read.
    """
    text  = Path(log_path).read_text(encoding='utf-8', errors='replace')
    lines = text.splitlines()

    # ------------------------------------------------------------------ model
    if 'GEOS-Chem Classic: Execution Test Results' in text:
        model = 'GCClassic'
    elif 'GCHP: Execution Test Results' in text:
        model = 'GCHP'
    else:
        raise ValueError(
            f"Cannot detect model type in '{log_path}'.\n"
            "Expected 'GEOS-Chem Classic: Execution Test Results' "
            "or 'GCHP: Execution Test Results'."
        )

    # ----------------------------------------------- header block (verbatim)
    # The header sits between the first and second '=====...' separator lines.
    sep_idx = [i for i, ln in enumerate(lines) if ln.startswith('=====')]
    if len(sep_idx) < 2:
        raise ValueError(
            f"Could not find two '===...' separator lines in '{log_path}'."
        )
    header_block = lines[sep_idx[0]: sep_idx[1] + 1]

    # Total test count – parsed from the header so it is available even
    # when the test-results section is absent (e.g. log captured mid-run).
    n_total = _find_int(
        r'Number of execution tests:\s*(\d+)',
        '\n'.join(header_block),
    )

    # -------------------------------------------- locate section boundaries
    exec_start = next(
        (i for i, ln in enumerate(lines) if ln.strip() == 'Execution tests:'),
        None,
    )
    summary_start = next(
        (i for i, ln in enumerate(lines)
         if ln.strip() == 'Summary of test results:'),
        None,
    )

    # ------------------------------- raw test lines (verbatim, for output)
    # Layout:
    #   exec_start + 0 : "Execution tests:"
    #   exec_start + 1 : "---...---" separator  (skip)
    #   exec_start + 2 : first test result line  ← begin slice here
    #   summary_start  : "Summary of test results:" ← end slice here
    test_lines_raw: list[str] = []
    if exec_start is not None and summary_start is not None:
        test_lines_raw = lines[exec_start + 2: summary_start]

    # ----------------------------------------------- parse individual tests
    # The non-greedy group before '\.{3,}' lets us handle test names that
    # contain a single embedded dot, e.g. gc_2x25_ModelE2.1_fullchem,
    # without confusing them with the dot-padding separator (which is always
    # three or more consecutive dots).
    _TEST_PAT = re.compile(
        r'^(.+?)\.{3,}Execute Simulation\.{3,}(\S+)\s*$',
        re.IGNORECASE,
    )
    tests: list[tuple[str, str]] = []
    for ln in test_lines_raw:
        m = _TEST_PAT.match(ln.strip())
        if m:
            result = m.group(2).strip().upper()
            if result not in ('PASS', 'FAIL'):
                result = 'OTHER'
            tests.append((m.group(1).strip(), result))

    # -------------------------------------------------------- summary counts
    summary_text = (
        '\n'.join(lines[summary_start:]) if summary_start is not None else ''
    )
    n_passed  = _find_int(r'Execution tests passed:\s*(\d+)',            summary_text)
    n_failed  = _find_int(r'Execution tests failed:\s*(\d+)',            summary_text)
    n_pending = _find_int(r'Execution tests not yet completed:\s*(\d+)', summary_text)

    # Fall back to computed total if the header line was absent.
    if n_total == 0:
        n_total = n_passed + n_failed + n_pending

    return dict(
        model          = model,
        header_block   = header_block,
        test_lines_raw = test_lines_raw,
        tests          = tests,
        n_total        = n_total,
        n_passed       = n_passed,
        n_failed       = n_failed,
        n_pending      = n_pending,
    )


# ---------------------------------------------------------------------------
# Diff-log parser
# ---------------------------------------------------------------------------

def parse_diff_log(log_path: str) -> dict:
    """
    Parse the output of ``diffTest.sh``.

    A test is flagged as **differing** when either:

    * A ``-> N difference(s) found in …`` line is present (nc4 / nc files), or
    * A raw ``diff …`` invocation line is present (text-file comparison,
      e.g. ``plane.log``).

    Lines matching ``-> No differences in …`` and file-path ``*`` lines are
    silently ignored.

    Parameters
    ----------
    log_path : str
        Path to the diff log file.

    Returns
    -------
    dict[str, bool]
        Maps each test name to ``True`` (has at least one difference) or
        ``False`` (zero-diff or not present in the log).

    Raises
    ------
    OSError
        If *log_path* cannot be read.
    """
    lines = Path(log_path).read_text(encoding='utf-8', errors='replace').splitlines()

    # "Checking <test_name>" at the start of a (possibly stripped) line.
    # \S+ rather than \w+ so that names with dots (ModelE2.1) are captured
    # in full.
    _CHECKING = re.compile(r'^Checking\s+(\S+)\s*$')

    # "-> N difference(s) found in OutputDir/Restarts"
    _HAS_DIFF = re.compile(r'->\s*\d+\s+difference', re.IGNORECASE)

    # Raw diff invocation: "diff -r file1 file2" or similar.
    # The word boundary \b prevents false matches on "differences".
    _RAW_DIFF = re.compile(r'^\s*diff\b')

    results : dict[str, bool] = {}
    current : str | None      = None
    has_diff: bool            = False

    for line in lines:
        # A new "Checking …" line starts a new test block; flush the previous.
        m = _CHECKING.match(line.strip())
        if m:
            if current is not None:
                results[current] = has_diff
            current  = m.group(1).strip()
            has_diff = False
            continue

        if current is not None:
            if _HAS_DIFF.search(line) or _RAW_DIFF.match(line):
                has_diff = True
            # "-> No differences in …" and "   * path" lines leave has_diff
            # unchanged (False unless already set True by an earlier line).

    # Flush the final test in the file.
    if current is not None:
        results[current] = has_diff

    return results


# ---------------------------------------------------------------------------
# Report helpers
# ---------------------------------------------------------------------------

def _exec_summary_sentence(exec_data: dict, model_long: str) -> str:
    """
    Return a single human-readable sentence summarising execution results.

    Examples
    --------
    "All GEOS-Chem Classic execution tests passed."
    "All GEOS-Chem Classic execution tests passed except gc_4x5_merra2_fullchem_APM."
    "GEOS-Chem Classic: 38 tests total – 2 failed, 1 not yet completed."
    """
    n_failed  = exec_data['n_failed']
    n_pending = exec_data['n_pending']
    n_total   = exec_data['n_total']
    failed    = [n for n, r in exec_data['tests'] if r == 'FAIL']

    if n_failed == 0 and n_pending == 0:
        return f'All {model_long} execution tests passed.'

    if n_failed == 1 and n_pending == 0 and failed:
        return (
            f'All {model_long} execution tests passed '
            f'except {failed[0]}.'
        )

    parts = []
    if n_failed:
        parts.append(f'{n_failed} failed')
    if n_pending:
        parts.append(f'{n_pending} not yet completed')
    return f'{model_long}: {n_total} tests total – ' + ', '.join(parts) + '.'


# ---------------------------------------------------------------------------
# Report generator
# ---------------------------------------------------------------------------

def generate_report(exec_data: dict, diff_data: dict, ref_label: str) -> str:
    """
    Build and return the full plain-text integration-test report.

    Parameters
    ----------
    exec_data : dict
        Output of :func:`parse_exec_log`.
    diff_data : dict
        Output of :func:`parse_diff_log`.
    ref_label : str
        Reference-version label (e.g. ``gcc.14.8.0-alpha.11``).

    Returns
    -------
    str
        The complete report as a single string (lines joined by ``\\n``,
        with a trailing newline).
    """
    model_long = (
        'GEOS-Chem Classic' if exec_data['model'] == 'GCClassic' else 'GCHP'
    )

    all_exec_passed = exec_data['n_failed'] == 0 and exec_data['n_pending'] == 0
    any_diffs       = any(diff_data.values())
    failed_tests    = [n for n, r in exec_data['tests'] if r == 'FAIL']
    other_tests     = [n for n, r in exec_data['tests'] if r == 'OTHER']
    diffing_tests   = sorted(k for k, v in diff_data.items() if v)

    out: list[str] = []

    # ------------------------------------------------------------------
    # Perfect run – compact output only
    # ------------------------------------------------------------------
    if all_exec_passed and not any_diffs:
        out.extend(_BANNER_LINES)
        out.append('')
        out.append(f'All {model_long} tests were zero-diff w/r/t {ref_label}.')
        return '\n'.join(out) + '\n'

    # ------------------------------------------------------------------
    # Section 1 – Execution test results
    # ------------------------------------------------------------------

    # Opening one-liner (e.g. "All … passed except …")
    out.append(_exec_summary_sentence(exec_data, model_long))
    out.append('')

    # Header block verbatim from the log (commit hashes, SLURM ID, etc.)
    out.extend(exec_data['header_block'])
    out.append('')

    # Test-result table verbatim from the log (preserves dot-padding alignment)
    out.append('Execution tests:')
    out.append(_SEP_LIGHT)
    out.extend(exec_data['test_lines_raw'])
    out.append('')

    # Summary counts
    out.append('Summary of test results:')
    out.append(_SEP_LIGHT)
    out.append(f'Execution tests passed:            {exec_data["n_passed"]}')
    out.append(f'Execution tests failed:            {exec_data["n_failed"]}')
    out.append(f'Execution tests not yet completed: {exec_data["n_pending"]}')

    # Show banner when execution is clean but diffs remain (common case after
    # an intentional science change that does not break any simulations).
    if all_exec_passed:
        out.append('')
        out.extend(_BANNER_LINES)

    # Explicit list of failed / not-yet-completed tests
    if failed_tests:
        out.append('')
        out.append('Tests that failed execution:')
        for t in failed_tests:
            out.append(f'  - {t}')

    if other_tests:
        out.append('')
        out.append('Tests with unknown / incomplete execution status:')
        for t in other_tests:
            out.append(f'  - {t}')

    # ------------------------------------------------------------------
    # Section 2 – Diff results
    # ------------------------------------------------------------------
    out.append('')
    out.append(_SEP_HEAVY)
    out.append(f'Diff results w/r/t {ref_label}:')
    out.append(_SEP_LIGHT)

    if not any_diffs:
        out.append(
            f'All {model_long} tests were zero-diff w/r/t {ref_label}.'
        )
    else:
        out.append(
            f'{model_long} tests with differences w/r/t {ref_label}:'
        )
        for t in diffing_tests:
            out.append(f'  - {t}')

    out.append('')
    return '\n'.join(out) + '\n'


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        prog='generate_inttest_report.py',
        description='Generate a plain-text GEOS-Chem integration-test report.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '\n'
            '  # GEOS-Chem Classic\n'
            '  python generate_inttest_report.py \\\n'
            '      --exec-log  gcc_exec_log.txt \\\n'
            '      --diff-log  gcc_diff_log.txt \\\n'
            '      --ref-label gcc.14.8.0-alpha.11 \\\n'
            '      --output    gcc_report.txt\n'
            '\n'
            '  # GCHP (separate invocation)\n'
            '  python generate_inttest_report.py \\\n'
            '      --exec-log  gchp_exec_log.txt \\\n'
            '      --diff-log  gchp_diff_log.txt \\\n'
            '      --ref-label gchp.14.8.0-alpha.11 \\\n'
            '      --output    gchp_report.txt\n'
        ),
    )
    parser.add_argument(
        '--exec-log',
        required=True,
        metavar='FILE',
        help='Execution test log produced by the integration test suite.',
    )
    parser.add_argument(
        '--diff-log',
        required=True,
        metavar='FILE',
        help='Diff log produced by diffTest.sh.',
    )
    parser.add_argument(
        '--ref-label',
        required=True,
        metavar='LABEL',
        help=(
            'Reference-version label as it appears in diff-log file paths '
            '(e.g. gcc.14.8.0-alpha.11 or gchp.14.8.0-alpha.11).'
        ),
    )
    parser.add_argument(
        '--output',
        metavar='FILE',
        default=None,
        help='Write the report to FILE instead of stdout.',
    )

    args = parser.parse_args()

    try:
        exec_data = parse_exec_log(args.exec_log)
        diff_data = parse_diff_log(args.diff_log)
        report    = generate_report(exec_data, diff_data, args.ref_label)
    except (ValueError, OSError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)

    if args.output:
        try:
            Path(args.output).write_text(report, encoding='utf-8')
            print(f'Report written to {args.output}', file=sys.stderr)
        except OSError as exc:
            print(f'ERROR writing to {args.output!r}: {exc}', file=sys.stderr)
            sys.exit(1)
    else:
        sys.stdout.write(report)


if __name__ == '__main__':
    main()
