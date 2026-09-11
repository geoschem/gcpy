#!/usr/bin/env python3
r"""
Script to scrape statistics from a 1-month GCHP benchmark run,
which can then be placed in the "GEOS-Chem 1-month Benchmark Stats"
Google spreadsheet.

Unlike GEOS-Chem Classic, GCHP does not wrap its executable in
``/usr/bin/time -v``, so there is no "Elapsed wall clock time" or
"Maximum resident set size" line in the public run log.  Peak memory
is instead read from the "Mem/Swap Used (MB)" lines that MAPL prints
to the run log throughout the simulation (the maximum such value,
found near the end of the run, is taken as the memory used).  Timing
info is not read from ``allPEs.log`` directly (it is only archived
inside a private tar.gz, not a public artifact); instead this script
reads the "Benchmark_Timers_<ref>_vs_<dev>.txt" table that gcpy's own
:mod:`gcpy.benchmark.modules.benchmark_scrape_gchp_timers` module already produces from
``allPEs.log`` during the benchmark pipeline, and publishes as a
public diff-table artifact. Because GCHP's MAPL timers are coarser
than GEOS-Chem Classic's (e.g. no separate Gas-phase chem/Photolysis/
Aerosol chem/Linearized chem/Unit conversions timers), only the
timers that GCHP actually reports are scraped and printed.

Examples
--------

.. code-block:: console

   $ conda activate gcpy_env
   $ python -m gcpy.benchmark.modules.benchmark_gchp_stats \
     14.5.0-alpha.5 \
     14.5.0-alpha.6
"""
import re
import sys
import requests
from gcpy.util import replace_whitespace, verify_variable_type

# ----------------------------------------------------------------------
# Global variables
# ----------------------------------------------------------------------

ROOT = "https://s3.amazonaws.com/benchmarks-cloud"

RUN_LOG_TEMPLATE = f"{ROOT}/benchmarks/1Mon/gchp/ID/RunGCHP.txt"

METRICS_TEMPLATE = \
    f"{ROOT}/diff-plots/1Mon/ID/BenchmarkResults/" + \
    "GCHP_version_comparison/Tables/OH_metrics.txt"

TIMERS_TEMPLATE = \
    f"{ROOT}/diff-plots/1Mon/ID/BenchmarkResults/" + \
    "GCHP_version_comparison/Tables/Benchmark_Timers_REF_vs_DEV.txt"

# Maps a column label to the corresponding GCHP MAPL timer, given
# as either an unqualified name (if unique in the timer hierarchy,
# e.g. "GC_CONV") or a fully-qualified dotted key (if the same name
# occurs more than once, e.g. "All.Run.GCHP.DYNAMICS"). See
# `flatten_timers` for how these dotted keys are built.
#
# Order follows the Summary section's own nesting: All.Run's
# top-level siblings first (SetService, Initialize, Finalize handled
# separately below), then Run's children (EXTDATA, HIST, GCHP's
# child DYNAMICS), then the GCHPchem chem-timer-table total, then
# its own individual (GC_*) timers -- GCHP does not report chemistry
# at Classic's granularity (no separate Gas-phase chem/Photolysis/
# Aerosol chem/Linearized chem/Unit conversions timers), so those
# have no entry and are not printed.
TIMERS = {
    "Run": "All.Run",
    "SetService": "All.SetService",
    "Initialize": "All.Initialize",
    "EXTDATA": "All.Run.EXTDATA",
    "HIST": "All.Run.HIST",
    "Transport": "All.Run.GCHP.DYNAMICS",
    "GEOS-Chem": "GCHPchem",
    "Convection": "GC_CONV",
    "Dry deposition": "GC_DRYDEP",
    "HEMCO": "GC_EMIS",
    "Boundary layer mixing": "GC_TURB",
    "All chemistry": "GC_CHEM",
    "Wet deposition": "GC_WETDEP",
    "Diagnostics": "GC_DIAGN",
    "Finalize": "All.Finalize",
}

# ----------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------

def print_stats(stats):
    """
    Prints OH metrics and timing statistics.

    Parameters
    ----------
    stats : dict
        Dictionary with statistics to print.
    """
    # Time and memory
    line = f"{stats['Wall Time']},,,{stats['Memory']},"

    # OH metrics
    line += f"{stats['Mean OH']},,{stats['CH3CCl3']},{stats['CH4']},,"

    # Timers
    for label in TIMERS:
        line += f"{stats[label]},"

    print(line)


def seconds_to_hms(seconds):
    """
    Converts a number of seconds to a "h:mm:ss" string, matching the
    format GEOS-Chem Classic gets for free from ``/usr/bin/time -v``.

    Parameters
    ----------
    seconds : float
        Number of seconds to convert.

    Returns
    -------
    result : str
        Elapsed time in "h:mm:ss" format.
    """
    total = round(seconds)
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours}:{minutes:02d}:{secs:02d}"


def scrape_memory(text):
    """
    Extracts the peak memory usage (in GB) from the "Mem/Swap Used
    (MB)" lines that MAPL prints to the GCHP run log throughout the
    simulation.  The largest such value (typically found near the
    end of the run) is taken as the memory used.

    Parameters
    ----------
    text : str
        Text scraped from the GCHP run log file.

    Returns
    -------
    memory : str
        Peak memory usage in GB.
    """
    matches = re.findall(
        r"Mem/Swap Used \(MB\) at[^=]+=\s*([-\d.Ee+]+)\s+[-\d.Ee+]+",
        text,
    )
    if not matches:
        raise ValueError("Could not find any 'Mem/Swap Used (MB)' lines!")

    peak_mb = max(float(value) for value in matches)
    return str(peak_mb / 1.0e3)


def scrape_oh_metrics(text):
    """
    Extracts the OH metrics (mean OH, CH3CCl3 and CH4 lifetimes)
    from the given text.

    Parameters
    ----------
    text : str
        Text scraped from the "OH_metrics.txt" file.

    Returns
    -------
    stats : dict
        Dictionary with the OH metrics.
    """
    stats = {}
    for line_count, line in enumerate(reversed(text.splitlines())):
        if line_count == 2 and "Dev" in line:
            stats["CH4"] = line.split(":")[1].strip()
        if line_count == 10 and "Dev" in line:
            stats["CH3CCl3"] = line.split(":")[1].strip()
        if line_count == 18 and "Dev" in line:
            stats["Mean OH"] = line.split(":")[1].strip()

    return stats


def flatten_timers(text):
    """
    Parses the "Benchmark_Timers_<ref>_vs_<dev>.txt" table (which
    prints GCHP's nested MAPL timers using "--"-prefixed labels to
    denote nesting depth) into a flat dictionary keyed by the
    dot-separated ancestry chain of each timer, e.g.
    ``"All.Run.GCHP.DYNAMICS"``.

    Parameters
    ----------
    text : str
        Text scraped from the "Benchmark_Timers_<ref>_vs_<dev>.txt"
        file.

    Returns
    -------
    timers : dict
        Dictionary mapping the dotted timer name to a
        ``(ref_seconds, dev_seconds)`` tuple.
    """
    hdr = [""] * 12
    timers = {}

    for line in text.splitlines():
        stripped = line.strip()

        # Skip blank lines, separator lines ("---...") and headers
        if not stripped or stripped.strip("-") == "" or "Ref [s]" in stripped:
            continue

        fields = stripped.split()
        if len(fields) < 3:
            continue

        label = fields[0]
        try:
            ref_seconds = float(fields[1])
            dev_seconds = float(fields[2])
        except ValueError:
            continue

        # The number of prefixed "-" characters denotes nesting depth
        name = label.lstrip("-")
        depth = (len(label) - len(name)) // 2
        hdr[depth] = name
        key = ".".join(hdr[:depth + 1])
        timers[key] = (ref_seconds, dev_seconds)

    return timers


def get_timer(timers, name):
    """
    Looks up a single GCHP timer by its unqualified name, e.g.
    ``"GC_CONV"``.  Raises an error if the name is ambiguous (found
    at more than one place in the timer hierarchy); callers that
    need a specific occurrence of an ambiguous name (e.g.
    "DYNAMICS", which appears once per top-level Init/Run/Finalize
    phase) should look it up by its full dotted key instead.

    Parameters
    ----------
    timers : dict
        Dictionary as returned by `flatten_timers`.
    name : str
        Unqualified timer name to look up.

    Returns
    -------
    dev_seconds : float
        The "Dev" column value (in seconds) for this timer.
    """
    # Exact key match (e.g. top-level "GCHPchem" or "All")
    if name in timers:
        return timers[name][1]

    # Otherwise look up by the last segment of the dotted key
    matches = [key for key in timers if key.split(".")[-1] == name]
    if len(matches) != 1:
        err_msg = f"Found {len(matches)} matches for timer '{name}'!"
        raise ValueError(err_msg)

    return timers[matches[0]][1]


def scrape_timers(text):
    """
    Extracts the wall time and the GCHP timers listed in the
    module-level `TIMERS` dict from the given text.

    Parameters
    ----------
    text : str
        Text scraped from the "Benchmark_Timers_<ref>_vs_<dev>.txt"
        file.

    Returns
    -------
    stats : dict
        Dictionary with the wall time and timer statistics.
    """
    timers = flatten_timers(text)

    stats = {"Wall Time": seconds_to_hms(timers["All"][1])}

    for label, name in TIMERS.items():
        stats[label] = str(get_timer(timers, name))

    return stats


def get_text_from_web(url):
    """
    Returns the text from a file located on the web.

    Parameters
    ----------
    url : str
        URL of the file to be parsed.
    """
    try:
        text = requests.get(url, timeout=10).text
    except FileNotFoundError as exc:
        err_msg = f"Could not download {url} from AWS!"
        raise FileNotFoundError(err_msg) from exc

    return text


def main(ref_label, dev_label):
    """
    Main program.  Given the labels from two benchmark simulations
    (ref and dev), downloads the relevant files from AWS and scrapes
    statistics from them.

    Parameters
    ----------
    ref_label : str
        Label for the Ref version.
    dev_label : str
        Label for the Dev version.
    """
    verify_variable_type(ref_label, str)
    verify_variable_type(dev_label, str)

    # Replace whitespace in the ref and dev labels
    ref_label = replace_whitespace(ref_label)
    dev_label = replace_whitespace(dev_label)

    # Run IDs for the Dev run log and the Ref-vs-Dev diff artifacts
    dev_id = f"gchp-c24-1Mon-{dev_label}"
    diff_id = f"diff-gchp-c24-1Mon-{ref_label}-gchp-c24-1Mon-{dev_label}"

    # Wall time & memory come from the Dev run's own console log
    run_log_text = get_text_from_web(RUN_LOG_TEMPLATE.replace("ID", dev_id))
    stats = {"Memory": scrape_memory(run_log_text)}

    # OH metrics come from the Ref-vs-Dev diff table
    metrics_text = get_text_from_web(
        METRICS_TEMPLATE.replace("ID", diff_id)
    )
    stats.update(scrape_oh_metrics(metrics_text))

    # Wall time & timers come from the Ref-vs-Dev timing diff table
    timers_url = TIMERS_TEMPLATE.replace("ID", diff_id)
    timers_url = timers_url.replace("REF", f"gchp-c24-1Mon-{ref_label}")
    timers_url = timers_url.replace("DEV", dev_id)
    timers_text = get_text_from_web(timers_url)
    stats.update(scrape_timers(timers_text))

    print_stats(stats)

# ----------------------------------------------------------------------
# For use from the command line
# ----------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) != 3:
        ERR_MSG = "Usage: benchmark_gchp_stats.py REF-LABEL DEV-LABEL"
        raise ValueError(ERR_MSG)

    main(sys.argv[1], sys.argv[2])
