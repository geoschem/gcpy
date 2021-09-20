"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.
"""

import os
import warnings
import itertools
from distutils.version import LooseVersion
import yaml
import numpy as np
import pandas as pd
import xarray as xr
from joblib import Parallel, delayed
from tabulate import tabulate
import gcpy.util as util
from gcpy.plot import compare_single_level, compare_zonal_mean
from gcpy.regrid import create_regridders
from gcpy.grid import get_troposphere_mask
from gcpy.units import convert_units
import gcpy.constants as gcon

# Save warnings format to undo overwriting built into PyPDF2
warning_format = warnings.showwarning

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")

# YAML files
aod_spc = "aod_species.yml"
spc_categories = "benchmark_categories.yml"
emission_spc = "emission_species.yml"
emission_inv = "emission_inventories.yml"


def create_total_emissions_table(
        refdata,
        refstr,
        devdata,
        devstr,
        species,
        outfilename,
        ref_interval=[2678400.0],
        dev_interval=[2678400.0],
        template="Emis{}_",
        refmetdata=None,
        devmetdata=None,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates a table of emissions totals (by sector and by inventory)
    for a list of species in contained in two data sets.  The data sets,
    which typically represent output from two differnet model versions,
    are usually contained in netCDF data files.

    Args:
        refdata: xarray Dataset
            The first data set to be compared (aka "Reference" or "Ref").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        devdata: xarray Dataset
            The second data set to be compared (aka "Development" or "Dev").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        species: dict
            Dictionary containing the name of each species and the target
            unit that emissions will be converted to. The format of
            species is as follows:

                { species_name: target_unit", etc. }

            where "species_name" and "target_unit" are strs.
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.

    Keyword Args (optional):
        ref_interval: float
            The length of the ref data interval in seconds. By default, interval
            is set to the number of seconds in a 31-day month (86400 * 31),
            which corresponds to typical benchmark simulation output.
            Default value: [2678400.0]
        dev_interval: float
            The length of the dev data interval in seconds. By default, interval
            is set to the number of seconds in a 31-day month (86400 * 31),
            which corresponds to typical benchmark simulation output.
            Default value: [2678400.0]
        template: str
            Template for the diagnostic names that are contained both
            "Reference" and "Development" data sets.  If not specified,
            template will be set to "Emis{}", where {} will be replaced
            by the species name.
            Default value: "Emis{}_"
        ref_area_varname: str
            Name of the variable containing the grid box surface areas
            (in m2) in the ref dataset.
            Default value: 'AREA'
        dev_area_varname: str
            Name of the variable containing the grid box surface areas
            (in m2) in the dev dataset.
            Default value: 'AREA'
        refmetdata: xarray dataset
            Dataset containing ref meteorology and area
            Default value: None
        devmetdata: xarray dataset
            Dataset containing dev meteorology and area
            Default value: None
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species properties (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Make sure refdata and devdata are both xarray Dataset objects
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata argument must be an xarray Dataset!")
    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # Get ref area [m2]
    if "AREA" in refdata.data_vars.keys():
        refarea = refdata["AREA"]
    elif refmetdata is not None:
        refarea = refmetdata["Met_AREAM2"]
    else:
        msg = "AREA variable is not in the ref Dataset and " + \
              "optional dataset containing Met_AREAM2 is not passed!"
        raise ValueError(msg)

    # Get dev area [m2]
    if "AREA" in devdata.data_vars.keys():
        devarea = devdata["AREA"]
    elif devmetdata is not None:
        devarea = devmetdata["Met_AREAM2"]
    else:
        msg = "AREA variable is not in the dev Dataset and optional " + \
              "dataset containing Met_AREAM2 is not passed!"
        raise ValueError(msg)

    # Load a YAML file containing species properties (such as
    # molecular weights), which we will need for unit conversions.
    # This is located in the "data" subfolder of this folder where
    # this benchmark.py file is found.
    properties_path = os.path.join(spcdb_dir, "species_database.yml")
    properties = yaml.load(open(properties_path), Loader=yaml.FullLoader)

    # ==================================================================
    # Get the list of emission variables for which we will print totals
    # ==================================================================

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refdata, devdata] = util.add_missing_variables(refdata, devdata)

    # Find all common variables between the two datasets
    # and get the lists of variables only in Ref and only in Dev,
    vardict = util.compare_varnames(refdata, devdata, quiet=True)
    cvars = vardict["commonvars"]
    refonly = vardict["refonly"]
    devonly = vardict["devonly"]

    # =================================================================
    # Open the file for output
    # =================================================================
    try:
        f = open(outfilename, "w")
    except FileNotFoundError:
        msg = "Could not open {} for writing!".format(outfilename)
        raise FileNotFoundError(msg)

    # =================================================================
    # Loop through all of the species are in species_dict
    # =================================================================
    for species_name, target_units in species.items():

        # Get a list of emission variable names for each species
        diagnostic_template = template.format(species_name)
        varnames = util.get_emissions_varnames(cvars, diagnostic_template)

        # Also add variables that might be in either Ref or Dev
        # but not the other.  This will allow us to print totals
        # for all species (and print NaN for the missing ones).
        if len(refonly) > 0:
            matching = [v for v in refonly if diagnostic_template in v]
            varnames = varnames + matching
        if len(devonly) > 0:
            matching = [v for v in devonly if diagnostic_template in v]
            varnames = varnames + matching

        # Sort the list again  to account for new variables added above
        varnames.sort()

        # If no emissions are found, then skip to next species
        if len(varnames) == 0:
            msg = "No emissions found for {} ... skippping"
            print(msg.format(species_name))
            continue

        # Check if there is a total emissions variable in the list
        vartot = [v for v in varnames if "_TOTAL" in v.upper()]

        # Push the total variable to the last list element
        # so that it will be printed last of all
        if len(vartot) == 1:
            varnames.append(varnames.pop(varnames.index(vartot[0])))

        # Title strings
        if "Inv" in template:
            print("Computing inventory totals for {}".format(species_name))
            title1 = "### Emissions totals for inventory {} [Tg]".format(
                species_name)
        else:
            print("Computing emissions totals for {}".format(species_name))
            title1 = "### Emissions totals for species {} [Tg]".format(species_name)

        title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

        # Print header to file
        print("#" * 83, file=f)
        print("{}{}".format(title1.ljust(80), "###"), file=f)
        print("{}{}".format(title2.ljust(80), "###"), file=f)
        print("#" * 83, file=f)
        print(
            "{}{}{}{}{}".format(
                " ".ljust(19),
                "Ref".rjust(20),
                "Dev".rjust(20),
                "Dev - Ref".rjust(14),
                "% diff".rjust(10),
            ),
            file=f)

        # =============================================================
        # Loop over all emissions variables corresponding to this
        # species and print their totals in Ref and Dev to the file.
        # =============================================================
        for v in varnames:

            if "Inv" in template:
                spc_name = v.split("_")[1]
            else:
                spc_name = species_name

            # Get a list of properties for the given species
            species_properties = properties.get(spc_name)

            # If no properties are found, then skip to next species
            if species_properties is None:
                print("No properties found for {} ... skippping".format(spc_name))
                continue

            # Convert units of Ref and Dev and save to numpy ndarray objects
            # (or set to NaN if the variable is not found in Ref or Dev)
            if v in refonly and v not in devonly:

                # Convert units of Ref
                refarray = convert_units(
                    refdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval=ref_interval,
                    area_m2=refarea,
                )

                # Set Dev to NaN (missing values) everywhere
                devarray = util.create_dataarray_of_nan(
                    name=refdata[v].name,
                    sizes=devdata.sizes,
                    coords=devdata.coords,
                    attrs=refdata[v].attrs,
                )

            elif v in devonly and v not in refonly:

                # Convert units of Dev
                devarray = convert_units(
                    devdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval=dev_interval,
                    area_m2=devarea,
                )

                # Set Ref to NaN (missing values) everywhere
                refarray = util.create_dataarray_of_nan(
                    name=devdata[v].name,
                    sizes=refdata.sizes,
                    coords=refdata.coords,
                    attrs=devdata[v].attrs,
                )

            else:

                # Convert units of both Ref and Dev
                refarray = convert_units(
                    refdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval=ref_interval,
                    area_m2=refarea,
                )
                devarray = convert_units(
                    devdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval=dev_interval,
                    area_m2=devarea,
                )

            # ==========================================================
            # Print emission totals for Ref and Dev
            # ==========================================================
            util.print_totals(refarray, devarray, f)

        # Add newlines before going to the next species
        print(file=f)
        print(file=f)

    # =================================================================
    # Close file
    # =================================================================
    f.close()


def create_global_mass_table(
        refdata,
        refstr,
        devdata,
        devstr,
        varlist,
        met_and_masks,
        label,
        trop_only=False,
        outfilename="GlobalMass_TropStrat.txt",
        verbose=False,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates a table of global masses for a list of species in contained in
    two data sets.  The data sets,  which typically represent output from two
    different model versions, are usually contained in netCDF data files.

    Args:
        refdata: xarray Dataset
            The first data set to be compared (aka "Reference").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        devdata: xarray Dataset
            The second data set to be compared (aka "Development").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        varlist: list of strings
            List of species concentation variable names to include
            in the list of global totals.
        met_and_masks: dict of xarray DataArray
            Dictionary containing the meterological variables and
            masks for the Ref and Dev datasets.
        label: str
            Label to go in the header string.  Can be used to
            pass the month & year.

    Keyword Args (optional):
        trop_only: bool
            Set this switch to True if you wish to print totals
            only for the troposphere.
            Default value: False (i.e. print whole-atmosphere totals).
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "GlobalMass_TropStrat.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species properties (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Make sure refdata and devdata are xarray Dataset objects
    if not isinstance(refdata, xr.Dataset):
        raise TypeError("The refdata argument must be an xarray Dataset!")
    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # Make sure required arguments are passed
    if varlist is None:
        raise ValueError('The "varlist" argument was not passed!')
    if met_and_masks is None:
        raise ValueError('The "met_and_masks" argument was not passed!')

    # Load a YAML file containing species properties (such as
    # molecular weights), which we will need for unit conversions.
    # This is located in the "data" subfolder of this current directory.2
    properties_path = os.path.join(spcdb_dir, "species_database.yml")
    properties = yaml.load(open(properties_path), Loader=yaml.FullLoader)

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Create file
    try:
        f = open(outfilename, "w")
    except FileNotFoundError:
        msg = "Could not open {} for writing!".format(outfilename)
        raise FileNotFoundError(msg)

    # Title strings
    if trop_only:
        title1 = "### Global mass (Gg) {} (Trop only)".format(label)
    else:
        title1 = "### Global mass (Gg) {} (Trop + Strat)".format(label)
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

    # Print header to file
    print("#" * 83, file=f)
    print("{}{}".format(title1.ljust(80), "###"), file=f)
    print("{}{}".format(title2.ljust(80), "###"), file=f)
    print("#" * 83, file=f)
    print(
        "{}{}{}{}{}".format(
            " ".ljust(19),
            "Ref".rjust(20),
            "Dev".rjust(20),
            "Dev - Ref".rjust(14),
            "% diff".rjust(10),
        ),
        file=f,
    )

    # ==================================================================
    # Print global masses for all species
    #
    # NOTE: By this point, all species will be in both Ref and Dev'
    # because we have added them in the calling routine
    # ==================================================================
    for v in varlist:

        # Get the species name
        spc_name = v.split("_")[1]

        # Get a list of properties for the given species
        species_properties = properties.get(spc_name)

        # If no properties are found, then skip to next species
        if species_properties is None:
            if verbose:
                msg = "No properties found for {} ... skippping"
                print(msg.format(spc_name))
            continue

        # Specify target units
        target_units = "Gg"
        mol_wt_g = species_properties.get("MW_g")
        if mol_wt_g is None:
            if verbose:
                msg = "No molecular weight found for {} ... skippping"
                print(msg.format(spc_name))
            continue

        # ==============================================================
        # Convert units of Ref and save to a DataArray
        # (or skip if Ref contains NaNs everywhere)
        # ==============================================================
        refarray = refdata[v]
        if not np.isnan(refdata[v].values).all():
            refarray = convert_units(
                refarray,
                spc_name,
                species_properties,
                target_units,
                area_m2=met_and_masks["Ref_Area"],
                delta_p=met_and_masks["Ref_Delta_P"],
                box_height=met_and_masks["Ref_BxHeight"],
            )

        # ==============================================================
        # Convert units of Dev and save to a DataArray
        # (or skip if Dev contains NaNs everywhere)
        # ==============================================================
        devarray = devdata[v]
        if not np.isnan(devdata[v].values).all():
            devarray = convert_units(
                devarray,
                spc_name,
                species_properties,
                target_units,
                area_m2=met_and_masks["Dev_Area"],
                delta_p=met_and_masks["Dev_Delta_P"],
                box_height=met_and_masks["Dev_BxHeight"],
            )

        # ==============================================================
        # Print global masses for Ref and Dev
        # (we will mask out tropospheric boxes in util.print_totals)
        # ==============================================================
        if trop_only:
            util.print_totals(
                refarray,
                devarray,
                f,
                masks=met_and_masks,
            )
        else:
            util.print_totals(
                refarray,
                devarray,
                f,
            )

    # ==================================================================
    # Close files
    # ==================================================================
    f.close()


def make_benchmark_conc_plots(
        ref,
        refstr,
        dev,
        devstr,
        dst="./benchmark",
        subdst=None,
        overwrite=False,
        verbose=False,
        collection="SpeciesConc",
        benchmark_type="FullChemBenchmark",
        plot_by_spc_cat=True,
        restrict_cats=[],
        plots=["sfc", "500hpa", "zonalmean"],
        use_cmap_RdBu=False,
        log_color_scale=False,
        sigdiff_files=None,
        normalize_by_area=False,
        cats_in_ugm3=["Aerosols", "Secondary_Organic_Aerosols"],
        areas=None,
        refmet=None,
        devmet=None,
        weightsdir='.',
        n_job=-1,
        second_ref=None,
        second_dev=None,
        time_mean=False,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates PDF files containing plots of species concentration
    for model benchmarking purposes.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str OR list of str
            A string to describe ref (e.g. version number)
            OR list containing [ref1str, ref2str] for diff-of-diffs plots
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str OR list of str
            A string to describe dev (e.g. version number)
            OR list containing [dev1str, dev2str] for diff-of-diffs plots

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where a PDF
            file containing plots will be written.
            Default value: ./benchmark
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
                collection: str
                        Name of collection to use for plotting.
                        Default value: "SpeciesConc"
                benchmark_type: str
                        A string denoting the type of benchmark output to plot,
                        either FullChemBenchmark or TransportTracersBenchmark.
                        Default value: "FullChemBenchmark"
        plot_by_spc_cat: logical
            Set this flag to False to send plots to one file rather
            than separate file per category.
            Default value: True
        restrict_cats: list of strings
            List of benchmark categories in benchmark_categories.yml to make
            plots for. If empty, plots are made for all categories.
            Default value: empty
        plots: list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']
        log_color_scale: bool
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False
        normalize_by_area: bool
            Set this flag to true to enable normalization of data
            by surfacea area (i.e. kg s-1 --> kg s-1 m-2).
                        Default value: False
        cats_in_ugm3: list of str
            List of benchmark categories to to convert to ug/m3
            Default value: ["Aerosols", "Secondary_Organic_Aerosols"]
        areas: dict of xarray DataArray:
            Grid box surface areas in m2 on Ref and Dev grids.
            Default value: None
        refmet: str
            Path name for ref meteorology
            Default value: None
        devmet: str
            Path name for dev meteorology
            Default value: None
        sigdiff_files: list of str
            Filenames that will contain the lists of species having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        second_ref: str
            Path name for a second "Ref" (aka "Reference") data set for
            diff-of-diffs plotting. This dataset should have the same model
            type and grid as ref.
            Default value: None
        second_dev: str
            Path name for a second "Ref" (aka "Reference") data set for
            diff-of-diffs plotting. This dataset should have the same model
            type and grid as ref.
            Default value: None
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
    """

    # NOTE: this function could use some refactoring;
    # abstract processing per category?

    # ==================================================================
    # Initialization and data read
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        os.mkdir(dst)

    # Define extra title text (usually a date string)
    # for the top-title of the plot
    if subdst is not None:
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Pick the proper function to read the data
    reader = util.dataset_reader(time_mean)

    # Open datasets
    refds = reader(ref, drop_variables=gcon.skip_these_vars)
    devds = reader(dev, drop_variables=gcon.skip_these_vars)

    # Open met datasets if passed as arguments
    refmetds = None
    devmetds = None
    if refmet:
        refmetds = reader(refmet, drop_variables=gcon.skip_these_vars)
    if devmet:
        devmetds = reader(devmet, drop_variables=gcon.skip_these_vars)

    # Determine if doing diff-of-diffs
    if second_ref is not None and second_dev is not None:
        diff_of_diffs = True
    else:
        diff_of_diffs = False

    # Open second datasets if passed as arguments (used for diff of diffs)
    if diff_of_diffs:
        second_refds = reader(second_ref, drop_variables=gcon.skip_these_vars)
        second_devds = reader(second_dev, drop_variables=gcon.skip_these_vars)
    else:
        second_refds = None
        second_devds = None

    # Compute the annual mean of datasets (if necessary)
    if time_mean:
        refds = util.dataset_mean(refds)
        devds = util.dataset_mean(devds)
        refmetds = util.dataset_mean(refmetds)
        devmetds = util.dataset_mean(devmetds)
        second_ref = util.dataset_mean(second_refds)
        second_dev = util.dataset_mean(second_devds)

    # Create regridding files if necessary while not in parallel loop
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # If we are normalizing by area, then merge the surface areas
    # on the Ref & Dev grids into the Ref & Dev datasets, but only
    # if they are not already there. The area variables should both
    # be called 'AREA' and be in units of m2.
    if normalize_by_area:
        if areas is not None:
            if 'AREA' not in refds.data_vars:
                refds = xr.merge([refds, areas["Ref"]])
            if 'AREA' not in devds.data_vars:
                devds = xr.merge([devds, areas["Dev"]])
            if diff_of_diffs:
                if 'AREA' not in second_refds.data_vars:
                    second_refds = xr.merge([second_refds, areas["Ref"]])
                if 'AREA' not in second_devds.data_vars:
                    second_devds = xr.merge([second_devds, areas["Dev"]])
        else:
            msg = "ERROR: normalize_by_area = True but " \
                + "the 'areas' argument was not passed!"
            raise ValueError(msg)

    # ==================================================================
    # If sending plots to one file then do all plots here and return.
    # Keep original units for plots.
    # ==================================================================
    if not plot_by_spc_cat:
        [refds, devds] = util.add_missing_variables(refds, devds)
        var_prefix = 'SpeciesConc_'
        varlist = [k for k in refds.data_vars.keys() if var_prefix in k]
        varlist.sort()

        # Surface
        pdfname = os.path.join(dst, 'SpeciesConc_Sfc.pdf')
        compare_single_level(refds, refstr, devds, devstr,
                             varlist=varlist,
                             pdfname=pdfname,
                             use_cmap_RdBu=use_cmap_RdBu,
                             log_color_scale=log_color_scale,
                             extra_title_txt=extra_title_txt,
                             normalize_by_area=normalize_by_area,
                             weightsdir=weightsdir,
                             second_ref=second_refds,
                             second_dev=second_devds,
                             spcdb_dir=spcdb_dir)

        util.add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                                  verbose=verbose)
        # 500 hPa
        pdfname = os.path.join(dst, 'SpeciesConc_500hPa.pdf')
        compare_single_level(refds, refstr, devds, devstr,
                             ilev=22,
                             varlist=varlist,
                             pdfname=pdfname,
                             use_cmap_RdBu=use_cmap_RdBu,
                             log_color_scale=log_color_scale,
                             normalize_by_area=normalize_by_area,
                             extra_title_txt=extra_title_txt,
                             weightsdir=weightsdir,
                             second_ref=second_refds,
                             second_dev=second_devds,
                             spcdb_dir=spcdb_dir)

        util.add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                                  verbose=verbose)
        # Zonal mean
        pdfname = os.path.join(dst, 'SpeciesConc_ZnlMn.pdf')
        compare_zonal_mean(refds, refstr, devds, devstr,
                           varlist=varlist,
                           pdfname=pdfname,
                           use_cmap_RdBu=use_cmap_RdBu,
                           log_color_scale=log_color_scale,
                           normalize_by_area=normalize_by_area,
                           extra_title_txt=extra_title_txt,
                           weightsdir=weightsdir,
                           second_ref=second_refds,
                           second_dev=second_devds,
                           spcdb_dir=spcdb_dir)

        util.add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                                  verbose=verbose)
        return

    # ==================================================================
    # Otherwise plot by categories. Convert units to ug/m3 for
    # aerosol categories: Aerosols and Secondary Organic Aerosols.
    # ==================================================================

    # FullChemBenchmark has lumped species (TransportTracers does not)
    if "FullChem" in benchmark_type:
        print("\nAdding lumped species to ref dataset")
        refds = util.add_lumped_species_to_dataset(refds)
        print("\nAdding lumped species to dev dataset")
        devds = util.add_lumped_species_to_dataset(devds)
        if diff_of_diffs:
            second_refds = util.add_lumped_species_to_dataset(second_refds)
            second_devds = util.add_lumped_species_to_dataset(second_devds)
        util.archive_lumped_species_definitions(dst)

    # Get the list of species categories
    catdict = util.get_species_categories(benchmark_type)
    util.archive_species_categories(dst)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = util.add_missing_variables(refds, devds)

    if diff_of_diffs:
        [refds, second_refds] = util.add_missing_variables(refds, second_refds)
        [devds, second_devds] = util.add_missing_variables(devds, second_devds)

    # Collection prefix
    coll_prefix = collection.strip() + "_"

    # ==================================================================
    # Create the plots!
    # ==================================================================

    # Use dictionaries to maintain order of significant difference categories
    dict_sfc = {}
    dict_500 = {}
    dict_zm = {}

    def createplots(filecat):
        cat_diff_dict = {'sfc': [], '500': [], 'zm': []}

        # Plots units ug/m3 for certain species categories
        convert_to_ugm3 = False
        if cats_in_ugm3 is not None and filecat in cats_in_ugm3:
            convert_to_ugm3 = True

        # Suppress harmless run-time warnings from all threads
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)

        # If restrict_cats list is passed,
        # skip all categories except those in the list
        if restrict_cats and filecat not in restrict_cats:
            return {filecat: cat_diff_dict}

        # Create a directory for each category.
        # If subdst is passed, then create a subdirectory in each
        # category directory (e.g. as for the 1-year benchmark).
        catdir = os.path.join(dst, filecat)
        if not os.path.isdir(catdir):
            os.mkdir(catdir)
        if subdst is not None:
            catdir = os.path.join(catdir, subdst)
            if not os.path.isdir(catdir):
                os.mkdir(catdir)

        varlist = []
        warninglist = []
        for subcat in catdict[filecat]:
            for spc in catdict[filecat][subcat]:
                varname = coll_prefix + spc
                if varname not in refds.data_vars or \
                   varname not in devds.data_vars:
                    warninglist.append(varname)
                    continue
                varlist.append(varname)
        if warninglist != []:
            msg = "\n\nWarning: variables in {} category not in dataset: {}"
            print(msg.format(filecat, warninglist))

        # -----------------------
        # Surface plots
        # -----------------------
        if "sfc" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_Surface_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(
                    catdir, "{}_Surface.pdf".format(filecat))

            diff_sfc = []
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                ilev=0,
                refmet=refmetds,
                devmet=devmetds,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                normalize_by_area=normalize_by_area,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_sfc,
                weightsdir=weightsdir,
                convert_to_ugm3=convert_to_ugm3,
                second_ref=second_refds,
                second_dev=second_devds,
                n_job=n_job,
                spcdb_dir=spcdb_dir
            )
            diff_sfc[:] = [v.replace(coll_prefix, "") for v in diff_sfc]
            cat_diff_dict['sfc'] = diff_sfc
            util.add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )

        # -----------------------
        # 500 hPa plots
        # -----------------------
        if "500hpa" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_500hPa_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(catdir, "{}_500hPa.pdf".format(filecat))

            diff_500 = []
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                ilev=22,
                refmet=refmetds,
                devmet=devmetds,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                normalize_by_area=normalize_by_area,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_500,
                weightsdir=weightsdir,
                convert_to_ugm3=convert_to_ugm3,
                second_ref=second_refds,
                second_dev=second_devds,
                n_job=n_job,
                spcdb_dir=spcdb_dir
            )
            diff_500[:] = [v.replace(coll_prefix, "") for v in diff_500]
            #dict_500[filecat] = diff_500
            cat_diff_dict['500'] = diff_500
            util.add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )

        # -----------------------
        # Zonal mean plots
        # -----------------------
        if "zonalmean" in plots or "zm" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_FullColumn_ZonalMean_{}.pdf".format(
                        filecat, subdst)
                )
            else:
                pdfname = os.path.join(
                    catdir, "{}_FullColumn_ZonalMean.pdf".format(filecat)
                )

            diff_zm = []
            compare_zonal_mean(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                refmet=refmetds,
                devmet=devmetds,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                normalize_by_area=normalize_by_area,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_zm,
                weightsdir=weightsdir,
                convert_to_ugm3=convert_to_ugm3,
                second_ref=second_refds,
                second_dev=second_devds,
                n_job=n_job,
                spcdb_dir=spcdb_dir
            )
            diff_zm[:] = [v.replace(coll_prefix, "") for v in diff_zm]
            #dict_zm = diff_zm
            cat_diff_dict['zm'] = diff_zm
            util.add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )

            # Strat_ZonalMean plots will use a log-pressure Y-axis, with
            # a range of 1..100 hPa, as per GCSC request. (bmy, 8/13/19)
            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_Strat_ZonalMean_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(catdir, "{}_Strat_ZonalMean.pdf".format(
                    filecat))

            compare_zonal_mean(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                refmet=refmetds,
                devmet=devmetds,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                pres_range=[1, 100],
                log_yaxis=True,
                extra_title_txt=extra_title_txt,
                log_color_scale=log_color_scale,
                normalize_by_area=normalize_by_area,
                convert_to_ugm3=convert_to_ugm3,
                weightsdir=weightsdir,
                second_ref=second_refds,
                second_dev=second_devds,
                n_job=n_job,
                spcdb_dir=spcdb_dir
            )
            util.add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )
        return {filecat: cat_diff_dict}

    # Create the plots in parallel
    results = Parallel(n_jobs=n_job)(
        delayed(createplots)(filecat) for _, filecat in enumerate(catdict)
    )

    dict_sfc = {list(result.keys())[0]: result[list(
        result.keys())[0]]['sfc'] for result in results}
    dict_500 = {list(result.keys())[0]: result[list(
        result.keys())[0]]['500'] for result in results}
    dict_zm = {list(result.keys())[0]: result[list(
        result.keys())[0]]['zm'] for result in results}

    # ==============================================================
    # Write the list of species having significant differences,
    # which we need to fill out the benchmark approval forms.
    # ==============================================================
    if sigdiff_files is not None:
        for filename in sigdiff_files:
            if "sfc" in plots:
                if "sfc" in filename:
                    with open(filename, "a+") as f:
                        for c, diff_list in dict_sfc.items():
                            print("* {}: ".format(c), file=f, end="")
                            for v in diff_list:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                        f.close()

            if "500hpa" in plots:
                if "500hpa" in filename:
                    with open(filename, "a+") as f:
                        for c, diff_list in dict_500.items():
                            print("* {}: ".format(c), file=f, end="")
                            for v in diff_list:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                        f.close()

            if "zonalmean" in plots or "zm" in plots:
                if "zonalmean" in filename or "zm" in filename:
                    with open(filename, "a+") as f:
                        for c, diff_list in dict_zm.items():
                            print("* {}: ".format(c), file=f, end="")
                            for v in diff_list:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                        f.close()

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    refds = xr.Dataset()
    devds = xr.Dataset()
    refmetds = xr.Dataset()
    devmetds = xr.Dataset()
    second_ref = xr.Dataset()
    second_dev = xr.Dataset()


def make_benchmark_emis_plots(
        ref,
        refstr,
        dev,
        devstr,
        dst="./benchmark",
        subdst=None,
        plot_by_spc_cat=False,
        plot_by_hco_cat=False,
        overwrite=False,
        verbose=False,
        flip_ref=False,
        flip_dev=False,
        log_color_scale=False,
        sigdiff_files=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates PDF files containing plots of emissions for model
    benchmarking purposes. This function is compatible with benchmark
    simulation output only. It is not compatible with transport tracers
    emissions diagnostics.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where
            PDF files containing plots will be written.
            Default value: './benchmark
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        plot_by_spc_cat: bool
            Set this flag to True to separate plots into PDF files
            according to the benchmark species categories (e.g. Oxidants,
            Aerosols, Nitrogen, etc.)  These categories are specified
            in the YAML file benchmark_species.yml.
            Default value: False
        plot_by_hco_cat: bool
            Set this flag to True to separate plots into PDF files
            according to HEMCO emissions categories (e.g. Anthro,
            Aircraft, Bioburn, etc.)
            Default value: False
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
        flip_ref: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        flip_dev: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        log_color_scale: bool
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False
         sigdiff_files: list of str
            Filenames that will contain the lists of species having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting.
            Value of -1 allows the application to decide.
            Default value: -1
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False

    Remarks:
        (1) If both plot_by_spc_cat and plot_by_hco_cat are
            False, then all emission plots will be placed into the
            same PDF file.

        (2) Emissions that are 3-dimensional will be plotted as
            column sums.
    """
    # =================================================================
    # Initialization and data read
    # =================================================================

    # Create destination folder if it does not exist
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite "\
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Create the "Emissions" category folder.  If subdst is passed,
    # then create a sub-folder (needed for the 1-year benchmarks).
    emisdir = os.path.join(dst, "Emissions")
    if not os.path.isdir(emisdir):
        os.mkdir(emisdir)
    if subdst is not None:
        emisdir = os.path.join(emisdir, subdst)
        if not os.path.isdir(emisdir):
            os.mkdir(emisdir)
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Get the function that will read the dataset
    reader = util.dataset_reader(time_mean)

    # Ref dataset
    try:
        refds = reader(ref, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Dev dataset
    try:
        devds = reader(dev, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = util.dataset_mean(refds)
        devds = util.dataset_mean(devds)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = util.add_missing_variables(refds, devds)

    # Create regridding files if necessary while not in parallel loop
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # Combine 2D and 3D variables into an overall list
    quiet = not verbose
    vardict = util.compare_varnames(refds, devds, quiet=quiet)
    vars2D = vardict["commonvars2D"]
    vars3D = vardict["commonvars3D"]
    varlist = vars2D + vars3D

    # ==================================================================
    # Compute column sums for 3D emissions
    # Make sure not to clobber the DataArray attributes
    # ==================================================================
    with xr.set_options(keep_attrs=True):
        for v in vars3D:
            if "lev" in refds[v].dims:
                refds[v] = refds[v].sum(dim="lev")
            if "lev" in devds[v].dims:
                devds[v] = devds[v].sum(dim="lev")

    # ==================================================================
    # If inputs plot_by* are both false, plot all emissions in same file
    # ==================================================================
    if not plot_by_spc_cat and not plot_by_hco_cat:
        if subdst is not None:
            pdfname = os.path.join(emisdir, "Emissions_{}.pdf".format(subdst))
        else:
            pdfname = os.path.join(emisdir, "Emissions.pdf")

        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        util.add_bookmarks_to_pdf(pdfname, varlist,
                                  remove_prefix="Emis", verbose=verbose)
        return

    # Get emissions variables (non-inventory), categories, and species
    emis_vars = [v for v in varlist if v[:4] == "Emis"]
    emis_cats = sorted(set([v.split("_")[1] for v in emis_vars]))
    emis_spc = sorted(set([v.split("_")[0][4:] for v in emis_vars]))

    # This is fixed in 12.3.2, comment out for now (bmy, 5/1/19)
    #    # Handle Bioburn and BioBurn as same categories (temporary until 12.3.1)
    #    emis_cats.remove('BioBurn')

    # Sort alphabetically (assume English characters)
    emis_vars.sort(key=str.lower)

    # ==================================================================
    # if plot_by_hco_cat is true, make a file for each HEMCO emissions
    # category that is in the diagnostics file
    #
    # Also write the list of emission quantities that have significant
    # diffs.  We'll need that to fill out the benchmark forms.
    # ==================================================================

    if plot_by_hco_cat:
        emisspcdir = os.path.join(dst, "Emissions")
        if not os.path.isdir(emisspcdir):
            os.mkdir(emisspcdir)
        if subdst is not None:
            emisspcdir = os.path.join(emisspcdir, subdst)
            if not os.path.isdir(emisspcdir):
                os.mkdir(emisspcdir)

        # for c in emis_cats:
        def createfile_hco_cat(c):
            # Handle cases of bioburn and bioBurn (temporary until 12.3.1)
            if c == "Bioburn":
                varnames = [k for k in emis_vars
                            if any(b in k for b in ["Bioburn", "BioBurn"])
                            ]
            else:
                varnames = [k for k in emis_vars if c in k]

            # Create the PDF name.  If subdst is passed, then also add
            # subdst to the file name (e.g. as for 1-year benchmarks).
            if subdst is not None:
                pdfname = os.path.join(
                    emisspcdir, "{}_Emissions_{}.pdf".format(c, subdst)
                )
            else:
                pdfname = os.path.join(
                    emisspcdir, "{}_Emissions.pdf".format(c))
            diff_dict = {}
            diff_emis = []
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varnames,
                ilev=0,
                pdfname=pdfname,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_emis,
                weightsdir=weightsdir,
                n_job=n_job,
                spcdb_dir=spcdb_dir
            )

            util.add_bookmarks_to_pdf(
                pdfname, varnames, remove_prefix="Emis", verbose=verbose
            )
            # Save the list of quantities with significant differences for
            # this category into the diff_dict dictionary for use below
            diff_emis[:] = [v.replace("Emis", "") for v in diff_emis]
            diff_emis[:] = [v.replace("_" + c, "") for v in diff_emis]
            diff_dict[c] = diff_emis
            return diff_dict

        results = Parallel(n_jobs=n_job)(delayed(createfile_hco_cat)(c)
                                         for c in emis_cats)

        dict_emis = {list(result.keys())[0]: result[list(result.keys())[0]]
                     for result in results}

        # =============================================================
        # Write the list of species having significant differences,
        # which we need to fill out the benchmark approval forms.
        # =============================================================
        if sigdiff_files is not None:
            for filename in sigdiff_files:
                if "emis" in filename:
                    with open(filename, "w+") as f:
                        for c, diff_list in dict_emis.items():
                            print("* {}: ".format(c), file=f, end="")
                            for v in diff_list:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                        f.close()

    # ==================================================================
    # if plot_by_spc_cat is true, make a file for each benchmark
    # species category with emissions in the diagnostics file
    # ==================================================================
    if plot_by_spc_cat:

        catdict = util.get_species_categories()
        # in case any emissions are skipped (for use in nested pdf bookmarks)
        warninglist = ([])
        # for checking if emissions species not defined in benchmark category
        # file
        allcatspc = ([])
        emisdict = {}  # used for nested pdf bookmarks
        # for i, filecat in enumerate(catdict):

        def createfile_bench_cat(filecat):
            # Get emissions for species in this benchmark category
            varlist = []
            emisdict[filecat] = {}
            catspc = []
            for subcat in catdict[filecat]:
                for spc in catdict[filecat][subcat]:
                    catspc.append(spc)
                    if spc in emis_spc:
                        emisdict[filecat][spc] = []
                        emisvars = [v for v in emis_vars
                                    if spc == v.split("_")[0][4:]]
                        for var in emisvars:
                            emisdict[filecat][spc].append(
                                var.replace("Emis", ""))
                            varlist.append(var)
            if not varlist:
                print(
                    "\nWarning: no emissions species in benchmark species" + \
                    "category {}".format(
                        filecat
                    )
                )
                return catspc

            # Use same directory structure as for concentration plots
            catdir = os.path.join(dst, filecat)
            if not os.path.isdir(catdir):
                os.mkdir(catdir)
            if subdst is not None:
                catdir = os.path.join(catdir, subdst)
                if not os.path.isdir(catdir):
                    os.mkdir(catdir)

            # Create emissions file for this benchmark species category
            # If subdst is passed, add it to the pdf name (e.g. as
            # is needed for the 1-year benchmarks).
            if subdst is not None:
                pdfname = os.path.join(
                    catdir, "{}_Emissions_{}.pdf".format(filecat, subdst)
                )
            else:
                pdfname = os.path.join(catdir, "{}_Emissions.pdf".format(
                    filecat))
            # Create the PDF
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                ilev=0,
                pdfname=pdfname,
                flip_ref=flip_ref,
                flip_dev=flip_dev,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
                weightsdir=weightsdir,
                n_job=n_job,
                spcdb_dir=spcdb_dir
            )
            util.add_nested_bookmarks_to_pdf(
                pdfname, filecat, emisdict, warninglist)
            return catspc
        results = Parallel(n_jobs=n_job)(
            delayed(createfile_bench_cat)(filecat)
            for i, filecat in enumerate(catdict)
        )

        allcatspc = [spc for result in results for spc in result]
        # Give warning if emissions species is not assigned a benchmark
        # category
        for spc in emis_spc:
            if spc not in allcatspc:
                print("Warning: species {} has emissions diagnostics but is not"
                      " in benchmark_categories.yml".format(spc))

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    refds = xr.Dataset()
    devds = xr.Dataset()


def make_benchmark_emis_tables(
        reflist,
        refstr,
        devlist,
        devstr,
        dst="./benchmark",
        refmet=None,
        devmet=None,
        overwrite=False,
        ref_interval=[2678400.0],
        dev_interval=[2678400.0],
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates a text file containing emission totals by species and
    category for benchmarking purposes.

    Args:
        reflist: list of str
             List with the path names of the emissions file or files
             (multiple months) that will constitute the "Ref"
             (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        devlist: list of str
             List with the path names of the emissions file or files
             (multiple months) that will constitute the "Dev"
             (aka "Development") data set
        devstr: str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./benchmark
        refmet: str
            Path name for ref meteorology
            Default value: None
        devmet: str
            Path name for dev meteorology
            Default value: None
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        ref_interval: list of float
            The length of the ref data interval in seconds. By default, interval
            is set to [2678400.0], which is the number of seconds in July
            (our 1-month benchmarking month).
            Default value: [2678400.0]
        dev_interval: list of float
            The length of the dev data interval in seconds. By default, interval
            is set to [2678400.0], which is the number of seconds in July
            (our 1-month benchmarking month).
            Default value: [2678400.0]
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository

    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Create destination folder
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        os.mkdir(dst)

    # Create the "Tables" category folder if it does not exist
    emisdir = os.path.join(dst, "Tables")
    if not os.path.isdir(emisdir):
        os.mkdir(emisdir)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Read the input datasets
    # Also read the meteorology datasets if passed. These are optional since it
    # the refds and devds have variable AREA already (always true) and
    # unit conversions do not require any meteorology.
    if len(reflist) == 1:
        reflist = [reflist]
    if len(devlist) == 1:
        devlist = [devlist]
    refmetds = None
    devmetds = None

    if LooseVersion(xr.__version__) < LooseVersion("0.15.0"):
        refds = xr.open_mfdataset(reflist, drop_variables=gcon.skip_these_vars)
        devds = xr.open_mfdataset(devlist, drop_variables=gcon.skip_these_vars)
        if refmet is not None:
            refmetds = xr.open_mfdataset(
                refmet, drop_variables=gcon.skip_these_vars)
        if devmet is not None:
            devmetds = xr.open_mfdataset(
                devmet, drop_variables=gcon.skip_these_vars)
    else:
        # , combine="nested", concat_dim="time")
        refds = xr.open_mfdataset(reflist, drop_variables=gcon.skip_these_vars)
        # , combine="nested", concat_dim="time")
        devds = xr.open_mfdataset(devlist, drop_variables=gcon.skip_these_vars)
        if refmet is not None:
            # , combine="nested", concat_dim="time")
            refmetds = xr.open_mfdataset(
                refmet, drop_variables=gcon.skip_these_vars)
        if devmet is not None:
            # , combine="nested", concat_dim="time")
            devmetds = xr.open_mfdataset(
                devmet, drop_variables=gcon.skip_these_vars)

    # ==================================================================
    # Create table of emissions
    # ==================================================================

    # Emissions species dictionary
    species = yaml.load(
        open(os.path.join(os.path.dirname(__file__), emission_spc)),
        Loader=yaml.FullLoader
    )
    inventories = yaml.load(
        open(os.path.join(os.path.dirname(__file__), emission_inv)),
        Loader=yaml.FullLoader
    )

    # Destination files
    file_emis_totals = os.path.join(emisdir, "Emission_totals.txt")
    file_inv_totals = os.path.join(emisdir, "Inventory_totals.txt")

    # Create table of emissions by species
    create_total_emissions_table(
        refds,
        refstr,
        devds,
        devstr,
        species,
        file_emis_totals,
        ref_interval,
        dev_interval,
        template="Emis{}_",
        refmetdata=refmetds,
        devmetdata=devmetds,
        spcdb_dir=spcdb_dir
    )

    # Create table of emissions by inventory
    create_total_emissions_table(
        refds,
        refstr,
        devds,
        devstr,
        inventories,
        file_inv_totals,
        ref_interval,
        dev_interval,
        template="Inv{}_",
        refmetdata=refmetds,
        devmetdata=devmetds,
        spcdb_dir=spcdb_dir
    )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    refds = xr.Dataset()
    devds = xr.Dataset()
    refmetds = xr.Dataset()
    devmetds = xr.Dataset()


def make_benchmark_jvalue_plots(
        ref,
        refstr,
        dev,
        devstr,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        local_noon_jvalues=False,
        plots=["sfc", "500hpa", "zonalmean"],
        overwrite=False,
        verbose=False,
        flip_ref=False,
        flip_dev=False,
        log_color_scale=False,
        sigdiff_files=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates PDF files containing plots of J-values for model
    benchmarking purposes.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
        varlist: list of str
            List of J-value variables to plot.  If not passed,
            then all J-value variables common to both dev
            and ref will be plotted.  The varlist argument can be
            a useful way of restricting the number of variables
            plotted to the pdf file when debugging.
            Default value: None
        dst: str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./benchmark.
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        local_noon_jvalues: bool
            Set this flag to plot local noon J-values.  This will
            divide all J-value variables by the JNoonFrac counter,
            which is the fraction of the time that it was local noon
            at each location.
            Default value: False
        plots: list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
        flip_ref: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        flip_dev: bool
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False
        log_color_scale: bool
            Set this flag to True if you wish to enable plotting data
            (not diffs) on a log color scale.
            Default value: False
        sigdiff_files: list of str
            Filenames that will contain the lists of J-values having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False

    Remarks:
         Will create 4 files containing J-value plots:
            (1 ) Surface values
            (2 ) 500 hPa values
            (3a) Full-column zonal mean values.
            (3b) Stratospheric zonal mean values
         These can be toggled on/off with the plots keyword argument.

         At present, we do not yet have the capability to split the
         plots up into separate files per category (e.g. Oxidants,
         Aerosols, etc.).  This is primarily due to the fact that
         we archive J-values from GEOS-Chem for individual species
         but not family species.  We could attempt to add this
         functionality later if there is sufficient demand.
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Create the destination folder if it does not exist
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in tht directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        os.mkdir(dst)

    # Get the function that will read file(s) into a Dataset
    reader = util.dataset_reader(time_mean)

    # Ref dataset
    try:
        refds = reader(ref, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Dev dataset
    try:
        devds = reader(dev, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = util.dataset_mean(refds)
        devds = util.dataset_mean(devds)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = util.add_missing_variables(refds, devds)

    # Create regridding files if necessary
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # Get a list of the 3D variables in both datasets
    if varlist is None:
        quiet = not verbose
        vardict = util.compare_varnames(refds, devds, quiet=quiet)
        cmn = vardict["commonvars3D"]

    # ==================================================================
    # Local noon or continuously-averaged J-values?
    # ==================================================================
    if local_noon_jvalues:

        # Get a list of local noon J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = "JNoon_"
        if varlist is None:
            varlist = [v for v in cmn if prefix in v]

        # Make sure JNoonFrac (fraction of times it was local noon
        # in each column) is present in both Ref and Dev datasets
        if "JNoonFrac" not in cmn:
            msg = "JNoonFrac is not common to Ref and Dev datasets!"
            raise ValueError(msg)

        # JNoon_* are cumulative sums of local noon J-values; we need
        # to divide these by JNoonFrac to get the average value
        refds = util.divide_dataset_by_dataarray(refds,
                                                 refds["JNoonFrac"], varlist)
        devds = util.divide_dataset_by_dataarray(devds,
                                                 devds["JNoonFrac"], varlist)

        # Subfolder of dst where PDF files will be printed
        catdir = "JValuesLocalNoon"

    else:

        # Get a list of continuously averaged J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = "Jval"
        if varlist is None:
            varlist = [v for v in cmn if prefix in v]

        # Subfolder of dst where PDF files will be printed
        catdir = "JValues"

    # ==================================================================
    # Create the plots
    # ==================================================================

    # Make the output folder if it doesn't exist.  If subdst is passed,
    # then create a sub-folder of this directory (e.g. which is needed
    # for the 1-year benchmarks)
    jvdir = os.path.join(dst, catdir)
    if not os.path.isdir(jvdir):
        os.mkdir(jvdir)
    if subdst is not None:
        jvdir = os.path.join(jvdir, subdst)
        if not os.path.isdir(jvdir):
            os.mkdir(jvdir)
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Surface plots
    if "sfc" in plots:
        if subdst is not None:
            pdfname = os.path.join(jvdir, "{}_Surface_{}.pdf".format(
                prefix, subdst))
        else:
            pdfname = os.path.join(jvdir, "{}_Surface.pdf".format(prefix))

        diff_sfc = []
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            ilev=0,
            pdfname=pdfname,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            sigdiff_list=diff_sfc,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        diff_sfc[:] = [v.replace(prefix, "") for v in diff_sfc]
        util.add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=prefix,
            verbose=verbose)

    # 500hPa plots
    if "500hpa" in plots:
        if subdst is not None:
            pdfname = os.path.join(jvdir, "{}_500hPa_{}.pdf".format(
                prefix, subdst))
        else:
            pdfname = os.path.join(jvdir, "{}_500hPa.pdf".format(prefix))

        diff_500 = []
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            ilev=22,
            pdfname=pdfname,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            sigdiff_list=diff_500,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        diff_500[:] = [v.replace(prefix, "") for v in diff_500]
        util.add_bookmarks_to_pdf(pdfname, varlist,
                                  remove_prefix=prefix, verbose=verbose)
    # Full-column zonal mean plots
    if "zonalmean" in plots:
        if subdst is not None:
            pdfname = os.path.join(
                jvdir, "{}_FullColumn_ZonalMean_{}.pdf".format(prefix, subdst)
            )
        else:
            pdfname = os.path.join(jvdir, "{}_FullColumn_ZonalMean.pdf".format(
                prefix))

        diff_zm = []
        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            sigdiff_list=diff_zm,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        diff_zm[:] = [v.replace(prefix, "") for v in diff_zm]
        util.add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=prefix,
            verbose=verbose)

        # Strat_ZonalMean plots will use a log-pressure Y-axis, with
        # a range of 1..100 hPa, as per GCSC request. (bmy, 8/13/19)
        if subdst is not None:
            pdfname = os.path.join(
                jvdir, "{}_Strat_ZonalMean_{}.pdf".format(prefix, subdst)
            )
        else:
            pdfname = os.path.join(
                jvdir, "{}_Strat_ZonalMean.pdf".format(prefix))

        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            pres_range=[1, 100],
            log_yaxis=True,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            extra_title_txt=extra_title_txt,
            log_color_scale=log_color_scale,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        util.add_bookmarks_to_pdf(pdfname, varlist,
                                  remove_prefix=prefix, verbose=verbose)

        # ==============================================================
        # Write the lists of J-values that have significant differences,
        # which we need to fill out the benchmark approval forms.
        # ==============================================================
        if sigdiff_files is not None:
            for filename in sigdiff_files:
                if "sfc" in plots:
                    if "sfc" in filename:
                        with open(filename, "a+") as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_sfc:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                            f.close()

                if "500" in plots:
                    if "500" in filename:
                        with open(filename, "a+") as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_500:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                            f.close()

                if "zonalmean" in plots or "zm" in plots:
                    if "zonalmean" in filename or "zm" in filename:
                        with open(filename, "a+") as f:
                            print("* J-Values: ", file=f, end="")
                            for v in diff_zm:
                                print("{} ".format(v), file=f, end="")
                            print(file=f)
                            f.close()

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    refds = xr.Dataset()
    devds = xr.Dataset()


def make_benchmark_aod_plots(
        ref,
        refstr,
        dev,
        devstr,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        overwrite=False,
        verbose=False,
        log_color_scale=False,
        sigdiff_files=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates PDF files containing plots of column aerosol optical
    depths (AODs) for model benchmarking purposes.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
        varlist: list of str
            List of AOD variables to plot.  If not passed, then all
            AOD variables common to both Dev and Ref will be plotted.
            Use the varlist argument to restrict the number of
            variables plotted to the pdf file when debugging.
            Default value: None
        dst: str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./benchmark.
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False
        log_color_scale: bool
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False
        sigdiff_files: list of str
            Filenames that will contain the list of quantities having
            having significant differences in the column AOD plots.
            These lists are needed in order to fill out the benchmark
            approval forms.
            Default value: None
        weightsdir: str
            Directory in which to place (and possibly reuse) xESMF regridder
            netCDF files.
            Default value: '.'
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
    """
    # ==================================================================
    # Initialization and also read data
    # ==================================================================

    # Create the destination directory if it does not exist
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        os.mkdir(dst)

    # Create the "Aerosols" directory as a subfolder of dst.
    # If subdst is passed, then create a subdirectory of the "Aerosols"
    # directory (e.g. which is needed for the 1-year benchmarks).
    aoddir = os.path.join(dst, "Aerosols")
    if not os.path.isdir(aoddir):
        os.mkdir(aoddir)
    if subdst is not None:
        aoddir = os.path.join(aoddir, subdst)
        if not os.path.isdir(aoddir):
            os.mkdir(aoddir)
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Get the function that will read file(s) into a dataset
    reader = util.dataset_reader(time_mean)

    # Read the Ref dataset
    try:
        refds = reader(ref, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Read the Dev dataset
    try:
        devds = reader(dev, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = util.dataset_mean(refds)
        devds = util.dataset_mean(devds)

    # Create regridding files if necessary
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # NOTE: GCHP diagnostic variable exports are defined before the
    # input.geos file is read.  This means "WL1" will not have been
    # replaced with "550nm" in the variable names.  Do this string
    # replace operation here, so that we can compare GCC and GCHP
    # data directly. (bmy, 4/29/19)
    with xr.set_options(keep_attrs=True):

        # Rename variables in the Ref dataset
        old2new = {}
        for v in refds.data_vars.keys():
            if "WL1" in v:
                newname = v.replace("WL1", "550nm")
                old2new[v] = newname
        refds = refds.rename(old2new)

        # Rename variables in the Dev dataset
        old2new = {}
        for v in devds.data_vars.keys():
            if "WL1" in v:
                newname = v.replace("WL1", "550nm")
                old2new[v] = newname
        devds = devds.rename(old2new)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = util.add_missing_variables(refds, devds)

    # Find common AOD variables in both datasets
    # (or use the varlist passed via keyword argument)
    if varlist is None:
        quiet = not verbose
        vardict = util.compare_varnames(refds, devds, quiet)
        cmn3D = vardict["commonvars3D"]
        varlist = [v for v in cmn3D if "AOD" in v and "_bin" not in v]

    # Dictionary and list for new display names
    newvars = yaml.load(
        open(os.path.join(os.path.dirname(__file__), aod_spc)),
        Loader=yaml.FullLoader
    )
    newvarlist = []

    # ==================================================================
    # Compute the total AOD by summing over the constituent members
    # ==================================================================

    # Take one of the variables so we can use its dims, coords,
    # attrs to create the DataArray object for total AOD
    v = varlist[0]

    # Create a DataArray to hold total column AOD
    # This is the same shape as the DataArray objects in refds
    reftot = xr.DataArray(
        np.zeros(refds[v].values.shape),
        name="AODTotal",
        dims=refds[v].dims,
        coords=refds[v].coords,
        attrs=refds[v].attrs,
    )

    # Create a DataArray to hold total column AOD
    # This is the same shape as the DataArray objects in devds
    devtot = xr.DataArray(
        np.zeros(devds[v].values.shape),
        name="AODTotal",
        dims=devds[v].dims,
        coords=devds[v].coords,
        attrs=devds[v].attrs,
    )

    # Save the variable attributes so that we can reattach them
    refattrs = reftot.attrs
    devattrs = devtot.attrs

    # Compute the sum of all AOD variables
    # Avoid double-counting SOA from aqueous isoprene, which is
    # already accounted for in AODHyg550nm_OCPI.  Also see
    # Github issue: https://github.com/geoschem/gcpy/issues/65
    for v in varlist:
        if "AODSOAfromAqIsoprene550nm" not in v:
            reftot = reftot + refds[v]
            devtot = devtot + devds[v]

    # Reattach the variable attributes
    reftot.name = "AODTotal"
    reftot.attrs = refattrs
    reftot.attrs["long_name"] = "Total aerosol optical depth"
    devtot.name = "AODTotal"
    devtot.attrs = devattrs
    devtot.attrs["long_name"] = "Total aerosol optical depth"

    # Merge these variables back into the dataset
    refds = xr.merge([refds, reftot])
    devds = xr.merge([devds, devtot])

    # Also add AODTotal to the list
    varlist.append("AODTotal")

    # ==================================================================
    # Compute column AODs
    # Create a new DataArray for each column AOD variable,
    # using the new display name, and preserving attributes.
    # Merge the new DataArrays back into the DataSets.
    # ==================================================================
    for v in varlist:

        # Get the new name for each AOD variable (it's easier to display)
        if v in newvars:
            newname = newvars[v]
            newvarlist.append(newname)
        else:
            raise ValueError("Could not find a display name for {}".format(v))

        # Don't clobber existing DataArray and Dataset attributes
        with xr.set_options(keep_attrs=True):

            # Add column AOD of newname to Ref
            array = refds[v].sum(dim="lev")
            array.name = newname
            array.attrs["units"] = "1"
            refds = xr.merge([refds, array])

            # Add column AOD of newname to Dev
            array = devds[v].sum(dim="lev")
            array.name = newname
            array.attrs["units"] = "1"
            devds = xr.merge([devds, array])

    # ==================================================================
    # Create the plots
    # ==================================================================
    if subdst is not None:
        pdfname = os.path.join(aoddir, "Aerosols_ColumnOptDepth_{}.pdf".format(
            subdst))
    else:
        pdfname = os.path.join(aoddir, "Aerosols_ColumnOptDepth.pdf")

    diff_aod = []
    compare_single_level(
        refds,
        refstr,
        devds,
        devstr,
        varlist=newvarlist,
        ilev=0,
        pdfname=pdfname,
        log_color_scale=log_color_scale,
        extra_title_txt=extra_title_txt,
        sigdiff_list=diff_aod,
        weightsdir=weightsdir,
        n_job=n_job,
        spcdb_dir=spcdb_dir
    )
    diff_aod[:] = [v.replace("Column_AOD_", "") for v in diff_aod]
    util.add_bookmarks_to_pdf(
        pdfname, newvarlist, remove_prefix="Column_AOD_", verbose=verbose
    )

    # ==================================================================
    # Write the list of AOD quantities having significant differences,
    # which we will need to fill out the benchmark forms.
    # ==================================================================
    if sigdiff_files is not None:
        for filename in sigdiff_files:
            if "sfc" in filename:
                with open(filename, "a+") as f:
                    print("* Column AOD: ", file=f, end="")
                    for v in diff_aod:
                        print("{} ".format(v), file=f, end="")
                    print(file=f)
                    f.close()

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    refds = xr.Dataset()
    devds = xr.Dataset()


def make_benchmark_mass_tables(
        ref,
        refstr,
        dev,
        devstr,
        varlist=None,
        dst="./benchmark",
        subdst=None,
        overwrite=False,
        verbose=False,
        label="at end of simulation",
        spcdb_dir=os.path.dirname(__file__),
        ref_met_extra='',
        dev_met_extra=''
):
    """
    Creates a text file containing global mass totals by species and
    category for benchmarking purposes.

    Args:
        reflist: str
            Pathname that will constitute
            the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: list of str
            Pathname that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.
        devstr: str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
        varlist: list of str
            List of variables to include in the list of totals.
            If omitted, then all variables that are found in either
            "Ref" or "Dev" will be included.  The varlist argument
            can be a useful way of reducing the number of
            variables during debugging and testing.
            Default value: None
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./benchmark
        subdst: str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False.
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository
        ref_met_extra: str
            Path to ref Met file containing area data for use with restart files
            which do not contain the Area variable.
            Default value: ''
        dev_met_extra: str
            Path to dev Met file containing area data for use with restart files
            which do not contain the Area variable.
            Default value: ''
    """

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        try:
            os.makedirs(dst)
        except FileExistsError:
            pass

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Read data
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=xr.SerializationWarning)
        refds = xr.open_dataset(ref, drop_variables=gcon.skip_these_vars)
        devds = xr.open_dataset(dev, drop_variables=gcon.skip_these_vars)

    # ==================================================================
    # Update GCHP restart dataset (if any)
    # ==================================================================

    # Ref
    if any(v.startswith("SPC_") for v in refds.data_vars.keys()):
        refds = util.rename_and_flip_gchp_rst_vars(refds)

    # Dev
    if any(v.startswith("SPC_") for v in devds.data_vars.keys()):
        devds = util.rename_and_flip_gchp_rst_vars(devds)

    # ==================================================================
    # Make sure that all necessary meteorological variables are found
    # ==================================================================

    # Find the area variables in Ref and Dev
    try:
        ref_area = util.get_area_from_dataset(refds)
    except ValueError:
        if ref_met_extra != '':
            ref_met_extra = xr.open_dataset(ref_met_extra)
            ref_area = util.get_area_from_dataset(ref_met_extra)
        else:
            raise ValueError(
                'Must pass Met data if using a restart file without area')
    try:
        dev_area = util.get_area_from_dataset(devds)
    except ValueError:
        if dev_met_extra != '':
            dev_met_extra = xr.open_dataset(dev_met_extra)
            dev_area = util.get_area_from_dataset(dev_met_extra)
        else:
            raise ValueError(
                'Must pass Met data if using a restart file without area')
    # Find required meteorological variables in Ref
    # (or exit with an error if we can't find them)
    metvar_list = ["Met_DELPDRY", "Met_BXHEIGHT", "Met_TropLev"]
    refmet = util.get_variables_from_dataset(refds, metvar_list)
    devmet = util.get_variables_from_dataset(devds, metvar_list)

    # ==================================================================
    # Make sure that all necessary species are found
    # ==================================================================

    # Get lists of variables names in datasets
    vardict = util.compare_varnames(refds, devds, quiet=(not verbose))
    commonvars = vardict["commonvars3D"]
    refonly = vardict['refonly']
    devonly = vardict['devonly']

    # Narrow down the lists to only include species
    commonspc = [v for v in commonvars if "SpeciesRst_" in v]
    refonlyspc = [v for v in refonly if v.startswith('SpeciesRst_')]
    devonlyspc = [v for v in devonly if v.startswith('SpeciesRst_')]

    # Add ref only species to dev dataset with all nan values
    if refonlyspc:
        for v in refonlyspc:
            devds[v] = devds[commonspc[0]]
            devds[v].data = np.full(devds[v].shape, np.nan)
            devds[v].attrs['units'] = refds[v].units
            commonspc.append(v)

    # Add dev only species to ref dataset with all nan values
    if devonlyspc:
        for v in devonlyspc:
            refds[v] = refds[commonspc[0]]
            refds[v].data = np.full(refds[v].shape, np.nan)
            devds[v].attrs['units'] = refds[v].units
            commonspc.append(v)

    # Set list of variables to print in mass table. If this list was passed
    # as argument, check that all the vars are now in commonspc to ensure
    # in both datasets.
    if varlist:
        for v in varlist:
            if v not in commonspc:
                raise ValueError(
                    '{} folder error: Variable {} in varlist passed to make_benchmark_mass_tables ' + \
                    'is not present in ref and dev datasets'.format(dst, v))
    else:
        varlist = commonspc

    # Sort the list of species to be printed alphabetically
    varlist.sort()

    # ==================================================================
    # Create the mask arrays for the troposphere for Ref and Dev
    # ==================================================================
    ref_tropmask = get_troposphere_mask(refmet)
    dev_tropmask = get_troposphere_mask(devmet)

    # ==================================================================
    # Create a dictionary to hold all of the meterological
    # variables and mask variables that we need to pass down
    # ==================================================================
    met_and_masks = {
        "Ref_Area": ref_area,
        "Dev_Area": dev_area,
        "Ref_Delta_P": refmet["Met_DELPDRY"],
        "Dev_Delta_P": devmet["Met_DELPDRY"],
        "Ref_BxHeight": refmet["Met_BXHEIGHT"],
        "Dev_BxHeight": devmet["Met_BXHEIGHT"],
        "Ref_TropMask": ref_tropmask,
        "Dev_TropMask": dev_tropmask,
    }

    # ==================================================================
    # Create global mass table
    # ==================================================================
    if subdst is not None:
        mass_filename = "GlobalMass_TropStrat_{}.txt".format(subdst)
    else:
        mass_filename = "GlobalMass_TropStrat.txt"
    mass_file = os.path.join(dst, mass_filename)
    create_global_mass_table(
        refds,
        refstr,
        devds,
        devstr,
        varlist,
        met_and_masks,
        label,
        outfilename=mass_file,
        verbose=verbose,
        spcdb_dir=spcdb_dir
    )

    # ==================================================================
    # Create tropospheric mass table
    # ==================================================================
    if subdst is not None:
        mass_filename = 'GlobalMass_Trop_{}.txt'.format(subdst)
    else:
        mass_filename = 'GlobalMass_Trop.txt'
    mass_file = os.path.join(dst, mass_filename)
    create_global_mass_table(
        refds,
        refstr,
        devds,
        devstr,
        varlist,
        met_and_masks,
        label,
        outfilename=mass_file,
        trop_only=True,
        verbose=verbose,
        spcdb_dir=spcdb_dir
    )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    refds = xr.Dataset()
    devds = xr.Dataset()


def make_benchmark_oh_metrics(
        ref,
        refmet,
        refstr,
        dev,
        devmet,
        devstr,
        dst="./benchmark",
        overwrite=False,
):
    """
    Creates a text file containing metrics of global mean OH, MCF lifetime,
    and CH4 lifetime for benchmarking purposes.

    Args:
        ref: str
            Path name of "Ref" (aka "Reference") data set file.
        refmet: str
            Path name of ref meteorology data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name of "Dev" (aka "Development") data set file.
            The "Dev" data set will be compared against the "Ref" data set.
        devmet: list of str
            Path name of dev meteorology data set.
        devstr: str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: "./benchmark"
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
    """

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        os.makedirs(dst)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    refds = xr.open_dataset(ref, drop_variables=gcon.skip_these_vars)
    devds = xr.open_dataset(dev, drop_variables=gcon.skip_these_vars)
    refmetds = xr.open_dataset(refmet, drop_variables=gcon.skip_these_vars)
    devmetds = xr.open_dataset(devmet, drop_variables=gcon.skip_these_vars)

    # ==================================================================
    # Get tropopause mask
    # ==================================================================
    # Find the area variables in Ref and Dev
    #ref_area = util.get_area_from_dataset(refds)
    #dev_area = util.get_area_from_dataset(devds)

    # Find required meteorological variables in Ref
    # (or exit with an error if we can't find them)
    metvar_list = [
        "Met_AD",
        "Met_AIRDEN",
        "Met_BXHEIGHT",
        "Met_T",
        "Met_TropLev",
        "FracOfTimeInTrop",
    ]
    refmet = util.get_variables_from_dataset(refmetds, metvar_list)
    devmet = util.get_variables_from_dataset(devmetds, metvar_list)

    # Create the mask arrays for the troposphere for Ref and Dev
    ref_tropmask = get_troposphere_mask(refmetds)
    dev_tropmask = get_troposphere_mask(devmetds)

    # ==================================================================
    # Compute mass-weighted OH in the troposphere
    # ==================================================================

    # Ref
    ref_oh_trop = np.ma.masked_array(refds["OHconcAfterChem"].values,
                                     ref_tropmask)
    ref_airmass_trop = np.ma.masked_array(refmetds["Met_AD"].values,
                                          ref_tropmask)
    ref_oh_mass = ref_oh_trop * ref_airmass_trop
    ref_total_ohmass = np.sum(ref_oh_mass)
    ref_total_airmass = np.sum(ref_airmass_trop)
    ref_mean_oh = (ref_total_ohmass / ref_total_airmass) / 1e5

    # Dev
    dev_oh_trop = np.ma.masked_array(devds["OHconcAfterChem"].values,
                                     dev_tropmask)
    dev_airmass_trop = np.ma.masked_array(devmetds["Met_AD"].values,
                                          dev_tropmask)
    dev_oh_mass = dev_oh_trop * dev_airmass_trop
    dev_total_ohmass = np.sum(dev_oh_mass)
    dev_total_airmass = np.sum(dev_airmass_trop)
    dev_mean_oh = (dev_total_ohmass / dev_total_airmass) / 1e5

    # Diff
    oh_diff = dev_mean_oh - ref_mean_oh
    oh_pctdiff = ((dev_mean_oh - ref_mean_oh) / ref_mean_oh) * 100.0

    # ==================================================================
    # Compute MCF and CH4 lifetimes
    # ==================================================================

    # Select only boxes that are purely tropospheric
    # This excludes influence from the stratosphere
    ref_timetrop_mask = refmetds["FracOfTimeInTrop"].values != 1.0
    dev_timetrop_mask = devmetds["FracOfTimeInTrop"].values != 1.0

    # Get grid box volumes [cm3] (trop + strat)
    ref_vol = (refmetds["Met_BXHEIGHT"] * refmetds["AREA"]) * 1e6
    dev_vol = (devmetds["Met_BXHEIGHT"] * devmetds["AREA"]) * 1e6

    # Get grid box volumes [cm3] (trop only)
    ref_vol_trop = np.ma.masked_array(ref_vol.values, ref_timetrop_mask)
    dev_vol_trop = np.ma.masked_array(dev_vol.values, dev_timetrop_mask)

    # Get MCF and CH4 density [molec/cm3] (trop + strat)
    # Assume that species is evenly distributed in air, with
    # a mixing ratio of 1. Thus species density = air density.
    ref_dens = refmetds["Met_AIRDEN"] / 1e6
    dev_dens = devmetds["Met_AIRDEN"] / 1e6

    # Get MCF and CH4 density [molec/cm3] (trop only)
    ref_dens_trop = np.ma.masked_array(ref_dens.values, ref_timetrop_mask)
    dev_dens_trop = np.ma.masked_array(dev_dens.values, dev_timetrop_mask)

    # Get temperature [K] (trop only)
    ref_temp = np.ma.masked_array(refmetds["Met_T"].values, ref_timetrop_mask)
    dev_temp = np.ma.masked_array(devmetds["Met_T"].values, dev_timetrop_mask)

    # Compute Arrhenius parameter K [cm3/molec/s]
    ref_mcf_k = 1.64e-12 * np.exp(-1520e0 / ref_temp)
    dev_mcf_k = 1.64e-12 * np.exp(-1520e0 / dev_temp)
    ref_ch4_k = 2.45e-12 * np.exp(-1775e0 / ref_temp)
    dev_ch4_k = 2.45e-12 * np.exp(-1775e0 / dev_temp)

    # Numerator: Total atmospheric (trop+strat) burden
    ref_num = np.sum(ref_dens.values * ref_vol.values)
    dev_num = np.sum(dev_dens.values * dev_vol.values)

    # Denominator: Loss rate in troposphere
    ref_mcf_denom = np.sum(ref_mcf_k * ref_oh_trop *
                           ref_dens_trop * ref_vol_trop)
    dev_mcf_denom = np.sum(dev_mcf_k * dev_oh_trop *
                           dev_dens_trop * dev_vol_trop)
    ref_ch4_denom = np.sum(ref_ch4_k * ref_oh_trop *
                           ref_dens_trop * ref_vol_trop)
    dev_ch4_denom = np.sum(dev_ch4_k * dev_oh_trop *
                           dev_dens_trop * dev_vol_trop)

    # Compute lifetimes [years]
    sec_to_year = 365.25 * 86400.0
    ref_mcf_lifetime = (ref_num / ref_mcf_denom) / sec_to_year
    dev_mcf_lifetime = (dev_num / dev_mcf_denom) / sec_to_year
    ref_ch4_lifetime = (ref_num / ref_ch4_denom) / sec_to_year
    dev_ch4_lifetime = (dev_num / dev_ch4_denom) / sec_to_year

    # Compute differences
    mcf_diff = dev_mcf_lifetime - ref_mcf_lifetime
    ch4_diff = dev_ch4_lifetime - ref_ch4_lifetime

    mcf_pctdiff = ((dev_mcf_lifetime - ref_mcf_lifetime) /
                   ref_mcf_lifetime) * 100.0
    ch4_pctdiff = ((dev_ch4_lifetime - ref_ch4_lifetime) /
                   ref_ch4_lifetime) * 100.0

    # ==================================================================
    #  Define function for writing metrics to file
    # ==================================================================

    def print_metrics_to_file(f, title1, title2, ref, dev, diff, pctdiff):
        print("#" * 79, file=f)
        print("{}{}".format(title1.ljust(76), "###"), file=f)
        print("{}{}".format(title2.ljust(76), "###"), file=f)
        print("#" * 79, file=f)
        print("{}{}{}{}".format("  Ref".ljust(15),
                                "Dev".ljust(13), "Dev - Ref".ljust(13),
                                "% diff".ljust(11),), file=f)
        print("{:11.6f}  {:11.6f}  {:11.6f}  {:9.4f}".format(ref, dev, diff,
                                                             pctdiff), file=f,)

    # ==================================================================
    # Print metrics to file
    # ==================================================================

    # Create file
    outfilename = os.path.join(dst, "OH_metrics.txt")
    f = open(outfilename, "w")

    # Write mean OH
    title1 = "### Global mass-weighted OH concentration [1e5 molec/cm3]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)
    print_metrics_to_file(f, title1, title2, ref_mean_oh, dev_mean_oh,
                          oh_diff, oh_pctdiff)

    # Write MCF lifetime
    title1 = "### MCF lifetime w/r/t tropospheric OH [years]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)
    print_metrics_to_file(f, title1, title2, ref_mcf_lifetime,
                          dev_mcf_lifetime, mcf_diff, mcf_pctdiff)

    # Write CH4 lifetime
    title1 = "### CH4 lifetime w/r/t tropospheric OH [years]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)
    print_metrics_to_file(f, title1, title2, ref_ch4_lifetime,
                          dev_ch4_lifetime, ch4_diff, ch4_pctdiff)

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    refds = xr.Dataset()
    devds = xr.Dataset()
    refmetds = xr.Dataset()
    devmetds = xr.Dataset()


def make_benchmark_wetdep_plots(
        ref,
        refstr,
        dev,
        devstr,
        collection,
        dst="./benchmark",
        datestr=None,
        overwrite=False,
        verbose=False,
        benchmark_type="TransportTracersBenchmark",
        plots=["sfc", "500hpa", "zonalmean"],
        log_color_scale=False,
        normalize_by_area=False,
        areas=None,
        refmet=None,
        devmet=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates PDF files containing plots of species concentration
    for model benchmarking purposes.

    Args:
        ref: str
            Path name for the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.
        devstr: str
            A string to describe dev (e.g. version number)
        collection: str
            String name of collection to plot comparisons for.

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where a PDF
            file containing plots will be written.
            Default value: ./benchmark
        datestr: str
            A string with date information to be included in both the
            plot pdf filename and as a destination folder subdirectory
            for writing plots
            Default value: None
        benchmark_type: str
            A string denoting the type of benchmark output to plot,
            either FullChemBenchmark or TransportTracersBenchmark.
            Default value: "FullChemBenchmark"
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.
        verbose: bool
            Set this flag to True to print extra informational output.
            Default value: False.
        plots: list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']
        normalize_by_area: bool
            Set this flag to true to enable normalization of data
            by surfacea area (i.e. kg s-1 --> kg s-1 m-2).
            Default value: False
        areas: dict of xarray DataArray:
            Grid box surface areas in m2 on Ref and Dev grids.
            Default value: None
        refmet: str
            Path name for ref meteorology
            Default value: None
        devmet: str
            Path name for dev meteorology
            Default value: None
        n_job: int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the
            application to decide.
            Default value: -1
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
    """

    #  Make sure destination directory exists
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        os.mkdir(dst)

    # Make a collection subdirectory
    targetdst = os.path.join(dst, collection)
    if not os.path.isdir(targetdst):
        os.mkdir(targetdst)

    # If datestr is passed, create a further subdirectory
    if datestr is not None:
        targetdst = os.path.join(targetdst, datestr)
        if not os.path.isdir(targetdst):
            os.mkdir(targetdst)

    # Get the function that will read file(s) into a dataset
    reader = util.dataset_reader(time_mean)

    # Open datasets
    refds = reader(ref, drop_variables=gcon.skip_these_vars)
    devds = reader(dev, drop_variables=gcon.skip_these_vars)

    # Open met datasets if passed as arguments
    refmetds = None
    devmetds = None
    if refmet is not None:
        refmetds = reader(refmet, drop_variables=gcon.skip_these_vars)
    if devmet is not None:
        devmetds = reader(devmet, drop_variables=gcon.skip_these_vars)

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = util.dataset_mean(refds)
        devds = util.dataset_mean(devds)
        if refmet is not None:
            refmetds = util.dataset_mean(refmetds)
        if devmet is not None:
            devmetds = util.dataset_mean(devmetds)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    # Turn this off for now since add_missing_variables inserts GCC area into
    # GCHP files, which causes problems with area normalization (ewl)
    #[refds, devds] = add_missing_variables(refds, devds)

    # Get list of variables in collection
    vardict = util.compare_varnames(refds, devds, quiet=True)
    varlist = [v for v in vardict["commonvars3D"] if collection + "_" in v]
    varlist.sort()

    # Surface plots
    if "sfc" in plots:
        if datestr is not None:
            plotfilename = "{}_Surface_{}.pdf".format(collection, datestr)
        else:
            plotfilename = "{}_Surface.pdf".format(collection)
        pdfname = os.path.join(targetdst, plotfilename)
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            ilev=0,
            refmet=refmetds,
            devmet=devmetds,
            pdfname=pdfname,
            normalize_by_area=normalize_by_area,
            extra_title_txt=datestr,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        util.add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=collection + '_',
            verbose=verbose)

    # 500 hPa plots
    if "500hpa" in plots:
        if datestr is not None:
            plotfilename = "{}_500hPa_{}.pdf".format(collection, datestr)
        else:
            plotfilename = "{}_500hPa.pdf".format(collection)
        pdfname = os.path.join(targetdst, plotfilename)
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            ilev=22,
            refmet=refmetds,
            devmet=devmetds,
            pdfname=pdfname,
            normalize_by_area=normalize_by_area,
            extra_title_txt=datestr,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        util.add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=collection + '_',
            verbose=verbose
        )

    # Zonal mean plots
    if "zonalmean" in plots or "zm" in plots:

        # Full column
        if datestr is not None:
            plotfilename = "{}_FullColumn_ZonalMean_{}.pdf".format(
                collection,
                datestr
            )
        else:
            plotfilename = "{}_FullColumn_ZonalMean.pdf".format(collection)
        pdfname = os.path.join(targetdst, plotfilename)
        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            refmet=refmetds,
            devmet=devmetds,
            pdfname=pdfname,
            log_color_scale=log_color_scale,
            normalize_by_area=normalize_by_area,
            extra_title_txt=datestr,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        util.add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=collection + '_',
            verbose=verbose
        )

        # Stratosphere
        if datestr is not None:
            plotfilename = "{}_Strat_ZonalMean_{}.pdf".format(
                collection,
                datestr
            )
        else:
            plotfilename = "{}_Strat_ZonalMean.pdf".format(collection)
        pdfname = os.path.join(targetdst, plotfilename)
        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            refmet=refmetds,
            devmet=devmetds,
            pdfname=pdfname,
            pres_range=[1, 100],
            log_yaxis=True,
            extra_title_txt=datestr,
            normalize_by_area=normalize_by_area,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_dir=spcdb_dir
        )
        util.add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=collection + '_',
            verbose=verbose
        )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    refds = xr.Dataset()
    devds = xr.Dataset()
    refmetds = xr.Dataset()
    devmetds = xr.Dataset()


def make_benchmark_aerosol_tables(
        devdir,
        devlist_aero,
        devlist_spc,
        devlist_met,
        devstr,
        year,
        days_per_mon,
        dst='./benchmark',
        overwrite=False,
        is_gchp=False,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Compute FullChemBenchmark aerosol budgets & burdens

    Args:
        devdir: str
            Path to development ("Dev") data directory
        devlist_aero: list of str
            List of Aerosols collection files (different months)
        devlist_spc: list of str
            List of SpeciesConc collection files (different months)
        devlist_met: list of str
            List of meteorology collection files (different months)
        devstr: str
            Descriptive string for datasets (e.g. version number)
        year: str
            The year of the benchmark simulation (e.g. '2016').
        days_per_month: list of int
            List of number of days per month for all months

    Keyword Args (optional):
        dst: str
            Directory where budget tables will be created.
            Default value: './benchmark'
        overwrite: bool
            Overwrite burden & budget tables? (default=True)
            Default value: False
        is_gchp: bool
            Whether datasets are for GCHP
            Default value: False
        spcdb_dir: str
            Directory of species_datbase.yml file
            Default value: Directory of GCPy code repository

    """

    # Create the plot directory hierarchy if it doesn't already exist
    if os.path.isdir(dst) and not overwrite:
        err_str = "Pass overwrite=True to overwrite files in that directory"
        print("Directory {} exists. {}".format(dst, err_str))
        return
    if not os.path.isdir(dst):
        os.makedirs(dst)

    # List of species (and subsets for the trop & strat)
    species_list = ["BCPI", "OCPI", "SO4", "DST1", "SALA", "SALC"]

    # Read the species database
    path = os.path.join(spcdb_dir, "species_database.yml")
    spcdb = yaml.load(open(path), Loader=yaml.FullLoader)

    # Molecular weights [g mol-1], as taken from the species database
    mw = {}
    for v in species_list:
        mw[v] = spcdb[v]["MW_g"]
    mw["Air"] = gcon.MW_AIR_g

    # Get the list of relevant AOD diagnostics from a YAML file
    path = os.path.join(os.path.dirname(__file__), "aod_species.yml")
    aod = yaml.load(open(path), Loader=yaml.FullLoader)
    aod_list = [v for v in aod.keys() if "Dust" in v or "Hyg" in v]
    # different names for GCHP
    if is_gchp:
        aod_list = [v.replace('550nm', 'WL1') for v in aod_list]
    # Descriptive names
    spc2name = {"BCPI": "Black Carbon",
                "DST1": "Dust",
                "OCPI": "Organic Carbon",
                "SO4": "Sulfate",
                "SALA": "Sea Salt (accum)",
                "SALC": "Sea Salt (coarse)"
                }

    # Read data collections
    if LooseVersion(xr.__version__) < LooseVersion("0.15.0"):
        ds_aer = xr.open_mfdataset(
            devlist_aero,
            data_vars=aod_list,
            compat='override',
            coords='all')
        ds_spc = xr.open_mfdataset(
            devlist_spc, drop_variables=gcon.skip_these_vars)
        ds_met = xr.open_mfdataset(
            devlist_met, drop_variables=gcon.skip_these_vars)
    else:
        ds_aer = xr.open_mfdataset(
            devlist_aero,
            data_vars=aod_list,
            compat='override',
            coords='all')  # ,
        # combine="nested", concat_dim="time")
        ds_spc = xr.open_mfdataset(devlist_spc,
                                   drop_variables=gcon.skip_these_vars)  # ,
        # combine="nested", concat_dim="time")
        ds_met = xr.open_mfdataset(devlist_met,
                                   drop_variables=gcon.skip_these_vars)  # ,
        # combine="nested", concat_dim="time")

    # Get troposphere mask
    tropmask = get_troposphere_mask(ds_met)

    # Get number of months
    n_mon = len(days_per_mon)

    # --------------------------------
    # Surface area
    # (kludgey but it works - revisit this)
    # --------------------------------

    # Get number of vertical levels
    N_LEVS = ds_spc.dims["lev"]

    if is_gchp:
        area = ds_met["Met_AREAM2"].values
        a = area.shape
        area_m2 = np.zeros([a[0], N_LEVS, a[1], a[2], a[3]])
        for t in range(n_mon):
            for k in range(N_LEVS):
                area_m2[t, k, :, :, :] = area[t, :, :, :]
        total_area_m2 = np.sum(area_m2[0, 0, :, :, :])
    else:
        area = ds_met["AREA"].values
        a = area.shape
        area_m2 = np.zeros([a[0], N_LEVS, a[1], a[2]])
        for t in range(n_mon):
            for k in range(N_LEVS):
                area_m2[t, k, :, :] = area[t, :, :]
        total_area_m2 = np.sum(area_m2[0, 0, :, :])

    # ------------------------------
    # Conversion factors and time increments
    # ------------------------------
    # v/v dry --> Tg
    vv_to_Tg = {}
    for spc in species_list:
        vv_to_Tg[spc] = ds_met["Met_AD"].values * (mw[spc] / mw["Air"]) * 1e-9

    # Days in the benchmark duration
    days_per_yr = np.sum(days_per_mon)

    # ------------------------------
    # Define function to print tables
    # ------------------------------
    def print_aerosol_metrics(data, species_list, filename, title, label):

        with open(filename, "w+") as f:

            # Print top header
            print("%" * 79, file=f)
            print(" {} for {} in {}".format(title, year, devstr), file=f)
            print(" (weighted by the number of days per month)", file=f)
            print("%" * 79, file=f)
            line = "\n" + " " * 40 + "Strat         Trop         Strat+Trop\n"
            line += " " * 40 + "-----------   ----------   ----------"
            print(line, file=f)

            # Print data
            for spc in species_list:
                line = "{} ({}) {} :  {:11.9f}   {:10.8f}   {:10.8f}\n".format(
                    spc2name[spc].ljust(17),
                    spc.ljust(4),
                    label,
                    data[spc + "_s"],
                    data[spc + "_t"],
                    data[spc + "_f"])
                print(line, file=f)

    # --------------------------------------
    # Compute aerosol burdens [Tg] and print
    # --------------------------------------

    # Table info
    filename = "{}/Aerosol_Burdens.txt".format(dst)
    if n_mon == 12:
        title = "Annual average global aerosol burdens"
    else:
        title = "Average global aerosol burdens across {} months".format(n_mon)
    label = "burden [Tg]"

    # Initialize
    q = {}
    q_sum_f = np.zeros(n_mon)
    q_sum_t = np.zeros(n_mon)
    q_sum_s = np.zeros(n_mon)
    burdens = {}

    # Define the axes we need to sum over to make monthly sums
    if is_gchp:
        sum_axes = (1, 2, 3, 4)
    else:
        sum_axes = (1, 2, 3)

    # Loop over species
    for spc in species_list:

        # Whole-atmosphere and trop-only quantities [g]
        # NOTE: DryDep is by nature trop-only
        varname = "SpeciesConc_" + spc
        q[spc + "_f"] = ds_spc[varname].values * vv_to_Tg[spc]
        q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"], tropmask)

        # Compute monthly sums, weighted by the number of days per month
        q_sum_f = np.sum(q[spc + "_f"], axis=sum_axes) * days_per_mon
        q_sum_t = np.sum(q[spc + "_t"], axis=sum_axes) * days_per_mon
        q_sum_s = q_sum_f - q_sum_t

        # Compute annual averages
        burdens[spc + "_f"] = np.sum(q_sum_f) / days_per_yr
        burdens[spc + "_t"] = np.sum(q_sum_t) / days_per_yr
        burdens[spc + "_s"] = np.sum(q_sum_s) / days_per_yr

    print_aerosol_metrics(burdens, species_list, filename, title, label)

    # -------------------------------------------
    # Compute average AOD's [Tg] and print
    # -------------------------------------------

    # Table info
    filename = "{}/Global_Mean_AOD.txt".format(dst)
    if n_mon == 12:
        title = "Annual average global AODs"
    else:
        title = "Average global AODs across {} months".format(n_mon)
    label = "mean AOD [1]"

    # Initialize
    q = {}
    q_sum_f = np.zeros(n_mon)
    q_sum_t = np.zeros(n_mon)
    q_sum_s = np.zeros(n_mon)
    aods = {}

    # Define axes to sum over, and total surface area
    if is_gchp:
        sum_axes = (1, 2, 3, 4)
    else:
        sum_axes = (1, 2, 3)

    # Loop over AOD variables
    for varname in aod_list:

        # Get the corresponding species name
        if "Dust" in varname:
            spc = "DST1"
        else:
            spc = varname.split("_")[1]

        # Whole-atmosphere AOD [1]
        q[spc + "_f"] = ds_aer[varname].values

        # Tropospheric-only AOD [1]
        q[spc + "_t"] = np.ma.masked_array(q[spc + "_f"], tropmask)

        # Create monthly sums, weighted by the number of days per month
        q_sum_f = np.sum(q[spc + "_f"] * area_m2, axis=sum_axes) * days_per_mon
        q_sum_t = np.sum(q[spc + "_t"] * area_m2, axis=sum_axes) * days_per_mon
        q_sum_s = q_sum_f - q_sum_t

        # Take annual averages
        aods[spc + "_f"] = np.sum(q_sum_f) / total_area_m2 / days_per_yr
        aods[spc + "_t"] = np.sum(q_sum_t) / total_area_m2 / days_per_yr
        aods[spc + "_s"] = np.sum(q_sum_s) / total_area_m2 / days_per_yr

    print_aerosol_metrics(aods, species_list, filename, title, label)

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    ds_aer = xr.Dataset()
    ds_spc = xr.Dataset()
    ds_met = xr.Dataset()


def make_benchmark_operations_budget(
        refstr,
        reffiles,
        devstr,
        devfiles,
        ref_interval,
        dev_interval,
        benchmark_type=None,
        label=None,
        col_sections=["Full", "Trop", "PBL", "Strat"],
        operations=["Chemistry", "Convection", "EmisDryDep",
                    "Mixing", "Transport", "WetDep"],
        compute_accum=True,
        require_overlap=False,
        dst='.',
        species=None,
        overwrite=True
):
    """
    Prints the "operations budget" (i.e. change in mass after
    each operation) from a GEOS-Chem benchmark simulation.

    Args:
        refstr: str
            Labels denoting the "Ref" versions
        reffiles: list of str
            Lists of files to read from the "Ref" version.
        devstr: str
            Labels denoting the "Dev" versions
        devfiles: list of str
            Lists of files to read from "Dev" version.
        interval: float
            Number of seconds in the diagnostic interval.

    Keyword Args (optional):
        benchmark_type: str
            "TransportTracersBenchmark" or "FullChemBenchmark".
            Default value: None
        label: str
            Contains the date or date range for each dataframe title.
            Default value: None
        col_sections: list of str
            List of column sections to calculate global budgets for. May
            include Strat eventhough not calculated in GEOS-Chem, but Full
            and Trop must also be present to calculate Strat.
            Default value: ["Full", "Trop", "PBL", "Strat"]
        operations: list of str
            List of operations to calculate global budgets for. Accumulation
            should not be included. It will automatically be calculated if
            all GEOS-Chem budget operations are passed and optional arg
            compute_accum is True.
            Default value: ["Chemistry","Convection","EmisDryDep",
                            "Mixing","Transport","WetDep"]
        compute_accum: bool
            Optionally turn on/off accumulation calculation. If True, will
            only compute accumulation if all six GEOS-Chem operations budgets
            are computed. Otherwise a message will be printed warning that
            accumulation will not be calculated.
            Default value: True
        require_overlap: bool
            Whether to calculate budgets for only species that are present in
            both Ref or Dev.
            Default value: False
        dst: str
            Directory where plots & tables will be created.
            Default value: '.' (directory in which function is called)
        species: list of str
            List of species for which budgets will be created.
            Default value: None (all species)
        overwrite: bool
            Denotes whether to overwrite existing budget file.
            Default value: True
    """

    # ------------------------------------------
    # Column sections
    # ------------------------------------------

    # Print info. Only allow Strat if Trop and Full are present
    print("Column sections:")
    for col_section in col_sections:
        print("  {}".format(col_section))
    n_sections = len(col_sections)
    compute_strat = False
    if "Strat" in col_sections:
        compute_strat = True
        if "Full" not in col_sections or "Trop" not in col_sections:
            msg = "Strat budget table requires Full and Trop column sections"
            raise ValueError(msg)
        print("*** Will compute Strat budget as difference of Full"
              " and Trop ***")

    # Make a subset of column sections not including Strat
    gc_sections = [r for r in col_sections if "Strat" not in r]

    # ------------------------------------------
    # Operations
    # ------------------------------------------

    # Rename the optional argument to be clear it is GEOS-Chem budget
    # operations
    gc_operations = operations

    # Handle whether to compute accumulation.
    all_operations = gc_operations
    if compute_accum and len(gc_operations) == 6:
        all_operations = gc_operations + ["ACCUMULATION"]
    n_ops = len(all_operations)

    # Print info
    print("Operations:")
    for all_operation in all_operations:
        print("  {}".format(all_operation))
    if compute_accum:
        if "ACCUMULATION" in all_operations:
            print("*** Will compute ACCUMULATION operation as sum of all "
                  "GEOS-Chem operations ***")
        else:
            print("***Will not compute ACCUMULATION since not all GEOS-Chem"
                  " operation budgets will be computed.")

    # ------------------------------------------
    # Read data
    # ------------------------------------------

    # Assume this will be annual budget if interval greater than 3e7 sec
    annual = ref_interval > 3.0e7

    # Read data from disk (either one month or 12 months)
    print('Opening ref and dev data')
    skip_vars = gcon.skip_these_vars
    if annual:
        if LooseVersion(xr.__version__) < LooseVersion("0.15.0"):
            ref_ds = xr.open_mfdataset(reffiles, drop_variables=skip_vars)
            dev_ds = xr.open_mfdataset(devfiles, drop_variables=skip_vars)
        else:
            # , combine="nested", concat_dim="time")
            ref_ds = xr.open_mfdataset(reffiles, drop_variables=skip_vars)
            # , combine="nested", concat_dim="time")
            dev_ds = xr.open_mfdataset(devfiles, drop_variables=skip_vars)
    else:
        ref_ds = xr.open_dataset(reffiles, drop_variables=skip_vars)
        dev_ds = xr.open_dataset(devfiles, drop_variables=skip_vars)

    # ------------------------------------------
    # Species
    # ------------------------------------------

    # Get information about variables in data files
    vardict = util.compare_varnames(ref_ds, dev_ds, quiet=True)
    refonly = vardict["refonly"]
    devonly = vardict["devonly"]
    cmnvars = vardict["commonvars2D"]

    # Reduce each list to only variables containing "Budget" and not "Strat"
    refonly = [v for v in refonly if "Budget" in v and "Strat" not in v]
    devonly = [v for v in devonly if "Budget" in v and "Strat" not in v]
    cmnvars = [v for v in cmnvars if "Budget" in v and "Strat" not in v]

    # Get the species list, depending on if species was passed as argument.
    if species is not None:
        spclist = species
    else:
        # For each column section, get the union or intersection (depending
        # on optional argument require_overlap) of all budget diagnostic
        # species and sort alphabetically. A species is counted even if only
        # present for as few as one operation.
        spclists = {}
        for gc_section in gc_sections:
            refspc = [v.split("_")[1] for v in refonly if gc_section in v]
            devspc = [v.split("_")[1] for v in devonly if gc_section in v]
            cmnspc = [v.split("_")[1] for v in cmnvars if gc_section in v]
            if require_overlap:
                section_spc = list(set(cmnspc))
            else:
                section_spc = list(set(refspc + devspc + cmnspc))
            section_spc.sort()
            spclists[gc_section] = section_spc

        # If calculating Strat, use intersection of Trop and Full
        if "Strat" in col_sections:
            spclists["Strat"] = list(set(spclists["Trop"] + spclists["Full"]))

        # For now, define one species list as species combined across all
        # column sections. Compute budgets of these species for each section.
        # NOTE: eventually would want to be able to compute budget for only
        # species per column section, not all species for all section. If that
        # is implemented then handling of optional species list needs to be
        # revisited.
        spclist = []
        for s in spclists:
            spclist = spclist + spclists[s]

    # Make list items unique and put in alphabetical order
    spclist = list(set(spclist))
    spclist.sort()
    n_spc = len(spclist)

    # ------------------------------------------
    # Concentration units and if a wetdep species
    # ------------------------------------------

    # Load a YAML file containing species properties
    spc_properties = yaml.load(open(os.path.join(os.path.dirname(__file__),
                                                 "species_database.yml")),
                               Loader=yaml.FullLoader)

    # Determine what the converted units and conversion factor should be
    # based on benchmark type and species (tracer) name. Assume raw data [kg/s]
    ref_conv_fac = {}
    dev_conv_fac = {}
    units = {}
    is_wetdep = {}
    for spc in spclist:

        # Identify wetdep species
        is_wetdep[spc] = None
        properties = spc_properties.get(spc)
        if properties is not None:
            is_wetdep[spc] = properties.get("Is_WetDep")

        # Unit conversion factors and units
        ref_conv_fac[spc] = ref_interval * 1e-6
        dev_conv_fac[spc] = dev_interval * 1e-6
        units[spc] = '[Gg]'
        if benchmark_type is not None:
            if "TransportTracers" in benchmark_type and "Tracer" not in spc:
                ref_conv_fac[spc] = ref_interval
                dev_conv_fac[spc] = dev_interval
                if annual:
                    units[spc] = '[kg/yr]'
                else:
                    units[spc] = '[kg]'
            elif annual:
                ref_conv_fac[spc] = ref_interval * 1e-9
                dev_conv_fac[spc] = dev_interval * 1e-9
                units[spc] = '[Tg/yr]'

    # ------------------------------------------
    # Create dataframe
    # ------------------------------------------
    columns = ["Column_Section", "Species", "Operation", "Ref_raw", "Dev_raw",
               "Units_converted", "Ref", "Dev", "Diff", "Pct_diff"]

    # Make column data to initialize with
    col_section = list(itertools.chain.from_iterable(
        itertools.repeat(x, n_ops * n_spc) for x in col_sections))
    col_spc = list(itertools.chain.from_iterable(
        itertools.repeat(x, n_ops) for x in spclist)) * n_sections
    col_ops = all_operations * n_sections * n_spc

    # Put the column data together into a dictionary
    data = {
        'Species': col_spc,
        'Operation': col_ops,
        'Column_Section': col_section,
    }

    # Create the dataframe from the data dictionary and column names list
    df = pd.DataFrame(data, columns=columns)

    # ------------------------------------------
    # Populate dataframe for GEOS-Chem operations and column sections
    # ------------------------------------------
    print('Calculating budgets for all data operations and column sections...')

    # Loop over sections (only those with data in files)
    for gc_section in gc_sections:

        # Loop over species in that section
        for i, spc in enumerate(spclist):

            # Keep track of progress
            if (i + 1) % 50 == 0:
                print('  {}: species {} of {}'.format(gc_section, i + 1,
                                                      n_spc))

            # Loop over operations (only those with data in files)
            for gc_operation in gc_operations:

                # Get the dataframe row to fill. Skip if not found.
                dfrow = (df["Column_Section"] == gc_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == gc_operation)
                if not any(dfrow):
                    continue

                # Get the variable name in the datasets
                varname = "Budget" + gc_operation + gc_section + "_" + spc

                # Calculate Ref and dev raw as global sum
                if varname in ref_ds.data_vars.keys():
                    refraw = ref_ds[varname].values.sum()
                elif gc_operation == "WetDep" and is_wetdep[spc] is None:
                    refraw = 0.0
                else:
                    refraw = np.nan
                if varname in dev_ds.data_vars.keys():
                    devraw = dev_ds[varname].values.sum()
                elif gc_operation == "WetDep" and is_wetdep[spc] is None:
                    devraw = 0.0
                else:
                    devraw = np.nan

                # Convert units
                refconv = refraw * ref_conv_fac[spc]
                devconv = devraw * dev_conv_fac[spc]

                # Calculate diff and % diff from conc with converted units
                if not np.isnan(refconv) and not np.isnan(devconv):
                    diff = devconv - refconv
                    try:
                        pctdiff = diff / refconv * 100
                    except BaseException:
                        pctdiff = np.nan
                else:
                    diff = np.nan
                    pctdiff = np.nan

                # Fill dataframe
                df.loc[dfrow, "Ref_raw"] = refraw
                df.loc[dfrow, "Dev_raw"] = devraw
                df.loc[dfrow, "Units_converted"] = units[spc]
                df.loc[dfrow, "Ref"] = refconv
                df.loc[dfrow, "Dev"] = devconv
                df.loc[dfrow, "Diff"] = diff
                df.loc[dfrow, "Pct_diff"] = pctdiff

    # ------------------------------------------
    # Compute Strat for each data operation (if applicable)
    # ------------------------------------------
    if compute_strat:
        print('Computing Strat budgets from Trop and Full...')

        # Loop over species
        for i, spc in enumerate(spclist):

            # Keep track of progress
            if (i + 1) % 50 == 0:
                print('  Strat: species {} of {}'.format(i + 1, n_spc))

            # Loop over operations (only those with data in files)
            for gc_operation in gc_operations:

                # Get the strat dataframe row to fill. Skip if not found.
                dfrow = (df["Column_Section"] == "Strat") \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == gc_operation)
                if not any(dfrow):
                    continue

                # Get the "Full" row to use in the calculation
                dfrow_full = (df["Column_Section"] == "Full") \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == gc_operation)
                if not any(dfrow_full):
                    continue

                # Get the "Trop" row to use in the calculation
                dfrow_trop = (df["Column_Section"] == "Trop") \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == gc_operation)
                if not any(dfrow_trop):
                    continue

                # Calculate strat concentrations as Full minus Trop. Do
                # not bother with raw value computation.
                refstrat = df.loc[dfrow_full, "Ref"].values[0] \
                    - df.loc[dfrow_trop, "Ref"].values[0]
                devstrat = df.loc[dfrow_full, "Dev"].values[0] \
                    - df.loc[dfrow_trop, "Dev"].values[0]

                # Calculate diff and % diff
                if not np.isnan(refstrat) and not np.isnan(devstrat):
                    diff = devstrat - refstrat
                    try:
                        pctdiff = diff / refstrat * 100
                    except BaseException:
                        pctdiff = np.nan
                else:
                    diff = np.nan
                    pctdiff = np.nan

                # Fill dataframe
                df.loc[dfrow, "Units_converted"] = units[spc]
                df.loc[dfrow, "Ref"] = refstrat
                df.loc[dfrow, "Dev"] = devstrat
                df.loc[dfrow, "Diff"] = diff
                df.loc[dfrow, "Pct_diff"] = pctdiff

    # ------------------------------------------
    # Compute Accumulation for each column section (if applicable)
    # ------------------------------------------
    if compute_accum:
        print('Computing ACCUMULATION operation budgets...')

        # Loop over all column sections
        for col_section in col_sections:

            # Loop over species
            for i, spc in enumerate(spclist):

                # Keep track of progress
                if (i + 1) % 50 == 0:
                    print('  {}: species {} of {}'.
                          format(col_section, i + 1, n_spc))

                # Get the accumulation dataframe row to fill.Skip if not found.
                dfrow = (df["Column_Section"] == col_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == "ACCUMULATION")
                if not any(dfrow):
                    continue

                # Get the rows to sum
                dfrows_to_sum = (df["Column_Section"] == col_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"].isin(gc_operations))

                # Sum the concentrations. If there is one or more NaN values
                # then set to NaN since not enough data for accumulation.
                refsum = df.loc[dfrows_to_sum, "Ref"].sum()
                devsum = df.loc[dfrows_to_sum, "Dev"].sum()

                # Calculate diff and % diff
                if not np.isnan(refsum) and not np.isnan(devsum):
                    diff = devsum - refsum
                    try:
                        pctdiff = diff / refsum * 100
                    except BaseException:
                        pctdiff = np.nan
                else:
                    diff = np.nan
                    pctdiff = np.nan

                # Fill dataframe
                df.loc[dfrow, "Units_converted"] = units[spc]
                df.loc[dfrow, "Ref"] = refsum
                df.loc[dfrow, "Dev"] = devsum
                df.loc[dfrow, "Diff"] = diff
                df.loc[dfrow, "Pct_diff"] = pctdiff

    #  Sanity check write to csv (for testing. Keep commented out otherwise)
    #df.to_csv('df.csv', na_rep='NA')

    # ------------------------------------------
    # Make budget file
    # ------------------------------------------

    # Create the target output directory hierarchy if it doesn't already exist
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. ".format(dst)
        msg += "Pass overwrite=True to overwrite files in that directory"
        raise ValueError(msg)
    if not os.path.isdir(dst):
        os.makedirs(dst)

    # Print budgets to file
    if label is not None:
        filename = "{}/Budgets_After_Operations_{}.txt".format(dst, label)
    else:
        filename = "{}/Budgets_After_Operations.txt".format(dst)
    with open(filename, "w+") as f:
        print("#" * 78, file=f)
        if label is not None and benchmark_type is not None:
            print("{} budgets for {}".format(benchmark_type, label),
                  file=f)
        else:
            print("Budgets across {}/{} sec".format(ref_interval, dev_interval), file=f)
        print("\n", file=f)
        print("NOTES:", file=f)
        msg = " - When using the non-local mixing scheme (default), "\
              "'Mixing' includes\n   emissions and dry deposition "\
              "applied below the PBL. 'EmisDryDep'\n   therefore only "\
              "captures fluxes above the PBL.\n - When using full mixing, "\
              "'Mixing' and 'EmisDryDep' are fully separated.\n - Budgets "\
              "are calculated as the sum of [kg/s] tendencies\n - Strat "\
              "budget is calculated as Full minus Trop\n - ACCUMULATION "\
              "is calculated as sum of all other operations"
        print(msg, file=f)
        print("#" * 78, file=f)
        print(file=f)

        # Loop over species
        for i, spc in enumerate(spclist):
            print("{} budgets (Ref={}; Dev={})".format(
                spc, refstr, devstr), file=f)

            # Print a table for each column section
            for col_section in col_sections:

                # Get the dataframe rows. Skip if none found.
                dfrows = (df["Column_Section"] == col_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"].isin(all_operations))
                if not any(dfrows):
                    continue

                # Print dataframe subset to file
                print(
                    "{} {} : {}".format(
                        col_section,
                        units[spc],
                        spc),
                    file=f)
                print(tabulate(df.loc[dfrows, ["Operation",
                                               "Ref",
                                               "Dev",
                                               "Diff",
                                               "Pct_diff"]],
                               headers='keys',
                               tablefmt='psql',
                               showindex=False,
                               floatfmt=(".5f", ".5f", ".5f", ".5f", ".5f"),
                               ), file=f
                      )
            print("\n", file=f)

    # ------------------------------------------
    # Clean up
    # ------------------------------------------
    df = pd.DataFrame()
    ref_ds = xr.Dataset()
    dev_ds = xr.Dataset()


def make_benchmark_mass_conservation_table(
        datafiles,
        runstr,
        dst="./benchmark",
        overwrite=False,
        spcdb_dir=os.path.dirname(__file__)
):
    """
    Creates a text file containing global mass of the PassiveTracer
    from Transport Tracer simulations across a series of restart files.

    Args:
        datafiles: list of str
            Path names of restart files.
        runstr: str
            Name to put in the filename and header of the output file
        refstr: str
            A string to describe ref (e.g. version number)
        dev: str
            Path name of "Dev" (aka "Development") data set file.
            The "Dev" data set will be compared against the "Ref" data set.
        devmet: list of str
            Path name of dev meteorology data set.
        devstr: str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: "./benchmark"
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
    """

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    if not os.path.isdir(dst):
        os.makedirs(dst)

    # Load a YAML file containing species properties (such as
    # molecular weights), which we will need for unit conversions.
    properties_path = os.path.join(spcdb_dir, "species_database.yml")
    properties = yaml.load(open(properties_path), Loader=yaml.FullLoader)

    # Get the species name
    spc_name = 'PassiveTracer'

    # Get a list of properties for the given species
    species_properties = properties.get(spc_name)

    # Specify target units
    target_units = "Tg"

    dates = []
    masses = []

    # ==================================================================
    # Calculate global mass for the tracer at all restart dates
    # ==================================================================
    for f in datafiles:
        ds = xr.open_dataset(f, drop_variables=gcon.skip_these_vars)

        # Save date in desired format
        #datestr = str(pd.to_datetime(ds.time.values[0]))
        #dates.append(datestr[:4] + '-' + datestr[5:7] + '-' + datestr[8:10])

        # Assume typical restart file name format, but avoid using dates
        # from within files which may be incorrect for the initial restart
        datestr = f.split('/')[-1].split('.')[2][:9]
        dates.append(datestr[:4] + '-' + datestr[4:6] + '-' + datestr[6:8])

        area = util.get_area_from_dataset(ds)
        # Select for GCC or GCHP
        delta_p = ds['Met_DELPDRY'] if 'Met_DELPDRY' in list(ds.data_vars) else ds['DELP_DRY']

        # ==============================================================
        # Convert units of Ref and save to a DataArray
        # (or skip if Ref contains NaNs everywhere)
        # ==============================================================
        # Select for GCC or GCHP
        if 'SpeciesRst_PassiveTracer' in list(ds.data_vars):
            attrs = ds['SpeciesRst_PassiveTracer'].attrs
            da = ds['SpeciesRst_PassiveTracer'].astype(np.float64)
            da.attrs = attrs
        else:
            attrs = ds['SPC_PassiveTracer'].attrs
            da = ds['SPC_PassiveTracer'].astype(np.float64)
            da.attrs = attrs
        da = convert_units(
            da,
            spc_name,
            species_properties,
            target_units,
            area_m2=area,
            delta_p=delta_p
        )

        # Save total global mass
        masses.append(np.sum(da.values))
        # Clean up
        ds = xr.Dataset()
        da = xr.DataArray()

    # Calclate max and min mass, absolute diff, percent diff
    max_mass = np.max(masses)
    min_mass = np.min(masses)
    # Convert absdiff to grams
    absdiff = (max_mass-min_mass) * 10**12
    pctdiff = max_mass/min_mass

    # ==================================================================
    # Print masses to file
    # ==================================================================
    # Create file
    outfilename = os.path.join(dst, "Passive_mass.txt")

    with open(outfilename, 'w') as f:
        titlestr = '  Global Mass of Passive Tracer in ' + runstr + '  '
        #headers
        print('%' * (len(titlestr)+4), file=f)
        print(titlestr, file=f)
        print('%' * (len(titlestr)+4), file=f)
        print('', file=f)
        print(' Date' + ' ' * 8 + 'Mass [Tg]', file=f)
        print(' ' + '-' * 10 + '  ' + '-' * 16, file=f)
        #masses
        for i in range(len(masses)):
            print(' {}  {:11.13f}'.format(dates[i], masses[i]), file=f)
        print(' ', file=f)
        print(' Summary', file=f)
        print(' ' + '-' * 30, file=f)
        print(' Max mass =  {:2.13f} Tg'.format(max_mass), file=f)
        print(' Min mass =  {:2.13f} Tg'.format(min_mass), file=f)
        print(' Abs diff =  {:>16.3f} g'.format(absdiff), file=f)
        print(' Pct diff =  {:>16.5f} %'.format(pctdiff), file=f)
