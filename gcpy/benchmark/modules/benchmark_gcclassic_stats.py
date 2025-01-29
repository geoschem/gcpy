#!/usr/bin/env python3
"""
Script to scrape statistics from a 1-month GEOS-Chem Classic benchmark run,
which can then be placed in the "GEOS-Chem 1-month Benchmark Stats"
Google spreadsheet.

Calling sequence:
$ python -m gcpy.benchmark.modules.benchmark_scrape_gcclassic_stats 14.5.0-alpha.5 14.5.0-alpha.6
"""
import sys
import requests
from gcpy.util import replace_whitespace, verify_variable_type

# ----------------------------------------------------------------------
# Global variables
# ----------------------------------------------------------------------

ROOT = "https://s3.amazonaws.com/benchmarks-cloud"

LOG_TEMPLATE = f"{ROOT}/benchmarks/1Mon/gcc/ID/RunGCC.txt"

METRICS_TEMPLATE = f"{ROOT}/diff-plots/1Mon/ID/BenchmarkResults/Tables/OH_metrics.txt"

TIMERS = [
    "GEOS-Chem                     :",
    "HEMCO                         :",
    "=> Gas-phase chem             :",
    "=> Photolysis                 :",
    "=> Aerosol chem               :",
    "=> Linearized chem            :",
    "Transport                     :",
    "Convection                    :",
    "Boundary layer mixing         :",
    "Dry deposition                :",
    "Wet deposition                :",
    "Diagnostics                   :",
    "Unit conversions              :",
]

# ----------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------

def print_stats(stats):
    """
    Prints OH metrics and timing statistics.

    Args
    stats (dict) : Dictionary with statistics to print
    """
    # Time and memory
    line = f"{stats['Wall Time']},,,{stats['Memory']},"

    # OH metrics
    line += f"{stats['Mean OH']},,{stats['CH3CCl3']},{stats['CH4']},,"

    # Timers
    timers = TIMERS
    for timer in timers:
        timer = format_timer(timer.split(":", maxsplit=1)[0])
        line += f"{stats[timer]},"

    print(line)


def format_timer(timer):
    """
    Strips spaces and preceding "=>" characters from a
    GEOS-Chem Classic timer name
    """
    return timer.strip().replace("=> ", "").replace(":", "")


def parse_timer(timer):
    """
    Extracts the timer name and time in seconds from the given text.

    Args
    timer (str) : Line of text with GEOS-Chem Classic timing output
    """
    sub_strings = timer.split(":")
    timer = format_timer(sub_strings[0])
    seconds = sub_strings[3].split()[1].strip()
    return timer, seconds


def scrape_stats(text):
    """
    Extracts timing statistics and OH metrics from the given text.

    Args
    text (str) : Text scraped from the log file and metrics file.
    """
    # Copy global variable to local for efficiency
    timers = TIMERS

    # Define empty dictionary for output and a counter
    stats = {}
    line_count = 0

    # Read the text backwards since the timers and OH are at the end
    for line in reversed(text.splitlines()):

        # Skip reading the rest of the file once we have
        # found the start of the timers section
        if "G E O S - C H E M   T I M E R S" in line:
            break

        # Look for the various metrics
        if line_count == 2 and "Dev" in line:
            stats["CH4"] = line.split(":")[1].strip()
        if line_count == 10 and "Dev" in line:
            stats["CH3CCl3"] = line.split(":")[1].strip()
        if line_count == 18 and "Dev" in line:
            stats["Mean OH"] = line.split(":")[1].strip()

        # Skip commands
        if "++ sed" in line:
            line_count += 1
            continue

        # Wall time
        if "wall clock" in line:
            stats["Wall Time"] = line.split("m:ss):")[1].strip()

        # Memory (GB)
        if "Maximum resident set size" in line:
            stats["Memory"] = str(float(line.split(":")[1]) / 1.0e6).strip()

        # GEOS-Chem Classic timers
        for timer in timers:
            if timer in line:
                timer, seconds = parse_timer(line)
                stats[timer] = str(round(float(seconds)))

        # Increment counter
        line_count += 1

    return stats


def get_text_from_web(url):
    """
    Returns the text from a file located on the web.

    Args
    url (str) : URL of the file to be parsed.
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
    (ref and dev), downloads the relevant files from AWS and passes
    the text to function "scrape_info" where it will be analyzed.

    Args
    ref_label (str) : Label for the Ref version
    dev_label (str) : Label for the Dev version
    """
    verify_variable_type(ref_label, str)
    verify_variable_type(dev_label, str)

    # Replace whitespace in the ref and dev labels
    ref_label = replace_whitespace(ref_label)
    dev_label = replace_whitespace(dev_label)

    # Scrape the log file text into a variable
    bmk_id = f"gcc-4x5-1Mon-{dev_label}"
    text = get_text_from_web(LOG_TEMPLATE.replace("ID", bmk_id))

    # Append the metrics file text
    bmk_id = f"diff-gcc-4x5-1Mon-{ref_label}-gcc-4x5-1Mon-{dev_label}"
    text += get_text_from_web(METRICS_TEMPLATE.replace("ID", bmk_id))

    # Scrape the relevant statistics from the text and print to stdout
    stats = scrape_stats(text)
    print_stats(stats)

# ----------------------------------------------------------------------
# For use from the command line
# ----------------------------------------------------------------------

if __name__ == '__main__':

    if len(sys.argv) != 3:
        ERR_MSG = "Usage: stats.py REF-LABEL DEV-LABEL"
        raise ValueError(ERR_MSG)

    main(sys.argv[1], sys.argv[2])
