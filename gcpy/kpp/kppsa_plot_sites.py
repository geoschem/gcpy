#!/usr/bin/env python3
"""
Script to visualize output from the KPP standalone box model.
"""

# Imports
import argparse
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from gcpy.util import verify_variable_type
from gcpy.kpp.kppsa_utils import \
    kppsa_get_file_list,  kppsa_get_unique_site_names,\
    kppsa_plot_one_page, kppsa_read_csv_files


def kppsa_plot_species_at_sites(
        ref_file_list,
        ref_label,
        dev_file_list,
        dev_label,
        species,
        pdfname,
):
    """
    Creates vertical profile plots of a given species
    from KPP-Standalone box model output.

    Args
    ref_file_list : list : KPP-Standalone log files for "Ref" version
    ref_label     : str  : Label for the "Ref" version
    dev_file_list : list : KPP-Standalone log files for "Dev" version
    dev_label     : str  : Label for the "Dev" version
    species       : str  : Name of the species to plot
    pdfname       : str  : Name of the output PDF file
    """
    verify_variable_type(ref_file_list, list)
    verify_variable_type(ref_label, str)
    verify_variable_type(dev_file_list, list)
    verify_variable_type(dev_label, str)
    verify_variable_type(species, str)
    verify_variable_type(pdfname, str)

    # Read data
    ref_data = kppsa_read_csv_files(ref_file_list)
    dev_data = kppsa_read_csv_files(dev_file_list)

    # Get a list of site names sorted from N to S
    site_names = kppsa_get_unique_site_names(ref_data)

    # Figure setup
    plt.style.use("seaborn-v0_8-darkgrid")
    rows_per_page = 3
    cols_per_page = 2
    plots_per_page = rows_per_page * cols_per_page

    # Open the plot as a PDF document
    if ".pdf" not in pdfname:
        pdfname += ".pdf"
    pdf = PdfPages(f"{pdfname}")

    # Loop over the number of obs sites that fit on a page
    for start in range(0, len(site_names), plots_per_page):
        end = start + plots_per_page - 1
        kppsa_plot_one_page(
            pdf,
            site_names[start:end+1],
            ref_data,
            ref_label,
            dev_data,
            dev_label,
            species,
            rows_per_page,
            cols_per_page,
            font_scale=1.0,
        )

    # Close the PDF file
    pdf.close()

    # Reset the plot style (this prevents the seaborn style from
    # being applied to other model vs. obs plotting scripts)
    plt.style.use("default")


def main():
    """
    Parses arguments from the command line and calls
    kppsa_plot_species_at_sites.

    Command-line arguments
    --refdir   : Folder with KPP-Standalone output from Ref model
    --reflabel : Plot label for the Ref model
    --devdir   : Folder with KPP-Standalone output from Dev model
    --devlabel : Plot label for the Dev model
    --pattern  : Look for filenames matching this pattern
    --species  : Species to plot
    --pdfname  : Name of the PDF file to be created
    """
    # Tell the parser which arguments to look for
    parser = argparse.ArgumentParser(
        description="Single-panel plotting example program"
    )
    parser.add_argument(
        "--refdir",
        metavar="REFDIR",
        type=str,
        required=True,
        help="Directory w/ KPP-Standalone log files (Ref version)"
    )
    parser.add_argument(
        "--reflabel",
        metavar="REFLABEL",
        type=str,
        required=False,
        help="Descriptive label for the Ref data",
        default="Ref"
    )
    parser.add_argument(
        "--devdir",
        metavar="DEVDIR",
        type=str,
        required=True,
        help="Directory w/ KPP-Standalone log files (Dev version)"
    )
    parser.add_argument(
        "--devlabel",
        metavar="DEVLABEL",
        type=str,
        required=False,
        help="Descriptive label for the Ref data",
        default="Dev"
    )
    parser.add_argument(
        "--pattern",
        metavar="PATTERN",
        type=str,
        required=False,
        help="Search for file names matching this pattern",
    )
    parser.add_argument(
        "--species",
        metavar="SPECIES",
        type=str,
        required=True,
        help="Species to plot"
    )
    parser.add_argument(
        "--pdfname",
        metavar="PDF-FILE-NAME",
        type=str,
        required=False,
        help="Name of the PDF file to be created",
        default="kppsa_output.pdf"
    )

    # Parse command-line arguments
    args = parser.parse_args()

    # Get a list of KPP-Standalone files matching the criteria (Ref)
    ref_file_list = kppsa_get_file_list(
        args.refdir,
        args.pattern,
    )
    if len(ref_file_list) == 0:
        msg = "Could not find any files matching {pattern} for Ref!"
        raise ValueError(msg)

    dev_file_list = kppsa_get_file_list(
        args.devdir,
        args.pattern,
    )
    if len(ref_file_list) == 0:
        msg = "Could not find any files matching {pattern} for Dev!"
        raise ValueError(msg)

    # Plot data
    kppsa_plot_species_at_sites(
        ref_file_list,
        args.reflabel,
        dev_file_list,
        args.devlabel,
        args.species,
        args.pdfname,
    )

if __name__ == '__main__':
    main()
