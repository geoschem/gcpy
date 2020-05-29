"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.
"""

import os
import shutil
import yaml
import numpy as np
import xarray as xr
from cartopy.mpl.geoaxes import GeoAxes  # for assertion
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfFileWriter, PdfFileReader, PdfFileMerger
from .plot import WhGrYlRd, compare_single_level, compare_zonal_mean
from .regrid import make_regridder_C2L, make_regridder_L2L, create_regridders
from .grid import GEOS_72L_grid, GEOS_47L_grid, get_grid_extents, call_make_grid, get_vert_grid, \
    get_vert_grid, get_pressure_indices, pad_pressure_edges, convert_lev_to_pres, get_troposphere_mask
import gcpy.util as util
from .units import convert_units
import gcpy.constants as gcon
from joblib import Parallel, delayed, cpu_count, parallel_backend
from multiprocessing import current_process
import warnings

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
    interval=[2678400.0],
    template="Emis{}_",
    ref_area_varname="AREA",
    dev_area_varname="AREA",
):
    """
    Creates a table of emissions totals (by sector and by inventory)
    for a list of species in contained in two data sets.  The data sets,
    which typically represent output from two differnet model versions,
    are usually contained in netCDF data files.

    Args:
    -----
        refdata : xarray Dataset
            The first data set to be compared (aka "Reference" or "Ref").

        refstr : str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).

        devdata : xarray Dataset
            The second data set to be compared (aka "Development" or "Dev").

        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).

        species : dict
            Dictionary containing the name of each species and the target
            unit that emissions will be converted to. The format of
            species is as follows:

                { species_name : target_unit", etc. }

            where "species_name" and "target_unit" are strs.

        outfilename : str
            Name of the text file which will contain the table of
            emissions totals.

    Keyword Args (optional):
    ------------------------
        interval : float
            The length of the data interval in seconds. By default, interval
            is set to the number of seconds in a 31-day month (86400 * 31),
            which corresponds to typical benchmark simulation output.

        template : str
            Template for the diagnostic names that are contained both
            "Reference" and "Development" data sets.  If not specified,
            template will be set to "Emis{}", where {} will be replaced
            by the species name.

        ref_area_varname : str
            Name of the variable containing the grid box surface areas
            (in m2) in the ref dataset.
            Default value: 'AREA'

        dev_area_varname : str
            Name of the variable containing the grid box surface areas
            (in m2) in the dev dataset.
            Default value: 'AREA'

    Remarks:
    --------
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

    # Make sure that the area variable is present in both refdata and devdata
    if ref_area_varname not in refdata.data_vars.keys():
        msg = "Area variable {} is not in the ref Dataset!".format(
            ref_area_varname)
        raise ValueError(msg)
    if dev_area_varname not in devdata.data_vars.keys():
        msg = "Area variable {} is not in the dev Dataset!".format(
            dev_area_varname)
        raise ValueError(msg)

    # Load a YAML file containing species properties (such as
    # molecular weights), which we will need for unit conversions.
    # This is located in the "data" subfolder of this folder where
    # this benchmark.py file is found.
    properties_path = os.path.join(os.path.dirname(__file__), "species_database.yml")
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
            title1 = "### Emissions totals for inventory {}".format(species_name)
        else:
            print("Computing emissions totals for {}".format(species_name))
            title1 = "### Emissions totals for species {}".format(species_name)

        title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

        # Print header to file
        print("#" * 79, file=f)
        print("{}{}".format(title1.ljust(76), "###"), file=f)
        print("{}{}".format(title2.ljust(76), "###"), file=f)
        print("#" * 79, file=f)
        print(
            "{}{}{}{}".format(
                " ".ljust(22), "Ref".rjust(20), "Dev".rjust(20), "Dev - Ref".rjust(15)
            ),
            file=f,
        )

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
                    interval,
                    refdata[ref_area_varname],
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
                    interval,
                    devdata[dev_area_varname],
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
                    interval,
                    refdata[ref_area_varname],
                )
                devarray = convert_units(
                    devdata[v],
                    spc_name,
                    species_properties,
                    target_units,
                    interval,
                    devdata[dev_area_varname],
                )

            # ==========================================================
            # Print emission totals for Ref and Dev
            # ==========================================================
            util.print_totals(refarray, refstr, devarray, devstr, f)

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
):
    """
    Creates a table of global masses for a list of species in contained in
    two data sets.  The data sets,  which typically represent output from two
    differnet model versions, are usually contained in netCDF data files.

    Args:
    -----
        refdata : xarray Dataset
            The first data set to be compared (aka "Reference").

        refstr : str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).

        devdata : xarray Dataset
            The second data set to be compared (aka "Development").

        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).

        varlist : list of strings
            List of species concentation variable names to include
            in the list of global totals.

        met_and_masks : dict of xarray DataArray
            Dictionary containing the meterological variables and
            masks for the Ref and Dev datasets.

        label : str
            Label to go in the header string.  Can be used to
            pass the month & year.

    Keyword Args (optional):
    ------------------------
        trop_only : bool
            Set this switch to True if you wish to print totals
            only for the troposphere.
            Default value: False (i.e. print whole-atmosphere totals).

        outfilename : str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "GlobalMass_TropStrat.txt"

        verbose : bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
    --------
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
    properties_path = os.path.join(os.path.dirname(__file__), "species_database.yml")
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
    print("#" * 79, file=f)
    print("{}{}".format(title1.ljust(76), "###"), file=f)
    print("{}{}".format(title2.ljust(76), "###"), file=f)
    print("#" * 79, file=f)
    print(
        "{}{}{}{}{}".format(
            " ".ljust(13),
            "Ref".rjust(20),
            "Dev".rjust(20),
            "Dev - Ref".rjust(15),
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
                refstr,
                devarray,
                devstr,
                f,
                mass_tables=True,
                masks=met_and_masks,
            )
        else:
            util.print_totals(
                refarray,
                refstr,
                devarray,
                devstr,
                f,
                mass_tables=True
            )

    # ==================================================================
    # Close files
    # ==================================================================
    f.close()


def create_budget_table(
    devdata,
    devstr,
    region,
    species,
    varnames,
    outfilename,
    interval=[2678400.0],
    template="Budget_{}",
):
    """
    Creates a table of budgets by species and component for a data set.

    Args:
    -----
        devdata : xarray Dataset
            The second data set to be compared (aka "Development").

        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).

        region : str
            Name of region for which budget will be computed.

        species : List of strings
            List of species  to include in budget tables.

        varnames : List of strings
            List of variable names in the budget diagnostics.

        outfilename : str
            Name of the text file which will contain the table of
            emissions totals.

    Keyword Args (optional):
    ------------------------
        interval : list of float
            The length of the data interval in seconds. By default, interval
            is set to [2678400.0], which is the number of seconds in July
            (our 1-month benchmarking month).

        template : str
            Template for the diagnostic names that are contained in the
            data set. If not specified, template will be set to "Budget_{}",
            where {} will be replaced by the species name.

    Remarks:
    --------
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.
    """

    # ==================================================================
    # Initialization
    # ==================================================================

    # Error check arguments
    if not isinstance(devdata, xr.Dataset):
        raise TypeError("The devdata argument must be an xarray Dataset!")

    # Open file for output
    try:
        f = open(outfilename, "w")
    except FileNotFoundError:
        msg = "Could not open {} for writing!".format(outfilename)
        raise FileNotFoundError(msg)

    # ==================================================================
    # Loop over species
    # ==================================================================
    for spc_name in species:

        # Title string
        title = "### {} budget totals for species {}".format(devstr, spc_name)

        # Write header to file
        print("#" * 79, file=f)
        print("{}{}".format(title.ljust(76), "###"), file=f)
        print("#" * 79, file=f)

        # Get variable names for this species
        spc_vars = [v for v in varnames if v.endswith("_" + spc_name)]

        for v in spc_vars:

            # Component name
            comp_name = v.replace("Budget", "")
            comp_name = comp_name.replace("_" + spc_name, "")
            comp_name = comp_name.replace(region, "")

            # Convert from kg/s to Tg
            devarray = devdata[v] * interval * 1e-9
            units = "Tg"

            # Compute sum
            total_dev = np.sum(devarray.values)

            # Write output
            print(
                "{} : {:13.6e} {}".format(comp_name.ljust(12), total_dev, units), file=f
            )

        # Add new lines before going to the next species
        print(file=f)
        print(file=f)

    # Close file
    f.close()

def make_benchmark_plots(
    ref,
    refstr,
    dev,
    devstr,
    dst="./1mo_benchmark",
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
    areas=None,
    weightsdir='.',
    n_job=-1,
    secondref=None,
    secondrefstr='',
    seconddev=None,
    seconddevstr=''
):
    """
    Creates PDF files containing plots of species concentration
    for model benchmarking purposes.

    Args:
    -----
        ref: str
            Path name for the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where a PDF
            file containing plots will be written.
            Default value: ./1mo_benchmark

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False.

        plot_by_spc_cat: logical
            Set this flag to False to send plots to one file rather
            than separate file per category.
            Default value: True

        restrict_cats : list of strings
            List of benchmark categories in benchmark_categories.yml to make
            plots for. If empty, plots are made for all categories.
            Default value: empty

        plots : list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']

        log_color_scale: boolean
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False

        normalize_by_area: bool
            Set this flag to true to enable normalization of data
            by surfacea area (i.e. kg s-1 --> kg s-1 m-2).

        areas : dict of xarray DataArray:
            Grid box surface areas in m2 on Ref and Dev grids.
            Default value: None

        sigdiff_files : list of str
            Filenames that will contain the lists of species having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None

        weightsdir : str
            Directory in which to place (and possibly reuse) xESMF regridder netCDF files.
            Default value: '.'

        n_job : int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the application to decide.
            Default value: -1

        second_ref: str
            Path name for a second "Ref" (aka "Reference") data set for diff-of-diffs plotting. 
            This dataset should have the same model type and grid as ref.
            Default value: None

        second_refstr : str
            A string to describe second_ref (e.g. version number)

        second_dev: str
            Path name for a second "Ref" (aka "Reference") data set for diff-of-diffs plotting. 
            This dataset should have the same model type and grid as ref.
            Default value: None

        second_devstr : str
            A string to describe second_dev (e.g. version number)

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
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Define extra title text (usually a date string)
    # for the top-title of the plot
    if subdst is not None:
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Ref dataset
    try:
        refds = xr.open_dataset(ref, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        msg ="Could not find Ref file: {}".format(ref)
        raise FileNotFoundError(msg)

    # Dev dataset
    try:
        devds = xr.open_dataset(dev, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        msg = "Could not find Dev file: {}!".format(dev)
        raise FileNotFoundError(msg)
        
    secondrefds = None
    seconddevds = None
    if secondref is not None:
        secondrefds = xr.open_dataset(secondref, drop_variables=gcon.skip_these_vars)
    if seconddev is not None:
        seconddevds = xr.open_dataset(seconddev, drop_variables=gcon.skip_these_vars)

    # Create regridding files if necessary while not in parallel loop
    [ _ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

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
            if secondref is not None and 'AREA' not in secondrefds.data_vars:
                secondrefds = xr.merge([secondrefds, areas["Ref"]])
            if seconddev is not None and 'AREA' not in seconddevds.data_vars:
                seconddevds = xr.merge([seconddevds, areas["Dev"]])
        else:
            msg = "ERROR: normalize_by_area = True but " \
                + "the 'areas' argument was not passed!"
            raise ValueError(msg)

    # ==================================================================
    # If sending plots to one file then do all plots here and return
    # ==================================================================
    if not plot_by_spc_cat:
        [refds, devds] = util.add_missing_variables(refds, devds)
        var_prefix = 'SpeciesConc_'
        varlist = [k for k in refds.data_vars.keys() if var_prefix in k]
        varlist.sort()

        # Surface
        pdfname = os.path.join(dst,'SpeciesConc_Sfc.pdf')
        compare_single_level(refds, refstr, devds, devstr,
                             varlist=varlist,
                             pdfname=pdfname,
                             use_cmap_RdBu=use_cmap_RdBu,
                             log_color_scale=log_color_scale,
                             extra_title_txt=extra_title_txt,
                             normalize_by_area=normalize_by_area,
                             weightsdir=weightsdir,
                             second_ref=secondrefds,
                             
                             second_dev=seconddevds)

        util.add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                             verbose=verbose)
        # 500 hPa
        pdfname = os.path.join(dst,'SpeciesConc_500hPa.pdf')
        compare_single_level(refds, refstr, devds, devstr,
                             ilev=22,
                             varlist=varlist,
                             pdfname=pdfname,
                             use_cmap_RdBu=use_cmap_RdBu,
                             log_color_scale=log_color_scale,
                             normalize_by_area=normalize_by_area,
                             extra_title_txt=extra_title_txt,
                             weightsdir=weightsdir,
                             second_ref=secondrefds, 
                             second_dev=seconddevds)

        util.add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                             verbose=verbose)
        # Zonal mean
        pdfname = os.path.join(dst,'SpeciesConc_ZnlMn.pdf')
        compare_zonal_mean(refds, refstr, devds, devstr,
                           varlist=varlist,
                           pdfname=pdfname,
                           use_cmap_RdBu=use_cmap_RdBu,
                           log_color_scale=log_color_scale,
                           normalize_by_area=normalize_by_area,
                           extra_title_txt=extra_title_txt,
                           weightsdir=weightsdir,
                           second_ref=secondrefds, 
                           second_dev=seconddevds)

        util.add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=var_prefix,
                             verbose=verbose)
        return

    # ==================================================================
    # Otherwise plot by categories
    # ==================================================================

    # FullChemBenchmark has lumped species (TransportTracers does not)
    if "FullChem" in benchmark_type:
        print("\nAdding lumped species to ref dataset:")
        refds = util.add_lumped_species_to_dataset(refds)
        print("\nAdding lumped species to dev dataset:")
        devds = util.add_lumped_species_to_dataset(devds)
        if secondrefds is not None:
            secondrefds = util.add_lumped_species_to_dataset(secondrefds)
        if seconddevds is not None:
            seconddevds = util.add_lumped_species_to_dataset(seconddevds)
        util.archive_lumped_species_definitions(dst)

    # Get the list of species categories
    catdict = util.get_species_categories(benchmark_type)
    util.archive_species_categories(dst)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = util.add_missing_variables(refds, devds)

    if secondrefds is not None:
        [refds, secondrefds] = util.add_missing_variables(refds, secondrefds)
    if seconddevds is not None:
        [devds, seconddevds] = util.add_missing_variables(devds, seconddevds)

    # Collection prefix
    coll_prefix = collection.strip() + "_"

    # ==================================================================
    # Create the plots!
    # ==================================================================

    # Use dictionaries to maintain order of significant difference categories
    dict_sfc = {}
    dict_500 = {}
    dict_zm = {}
    def createplots(i, filecat):
        cat_diff_dict = {'sfc' : [], '500' : [], 'zm' : []}

        # Suppress harmless run-time warnings from all threads
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)

        # If restrict_cats list is passed,
        # skip all categories except those in the list
        if restrict_cats and filecat not in restrict_cats:
            return {filecat : cat_diff_dict}

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
                pdfname = os.path.join(catdir, "{}_Surface.pdf".format(filecat))

            diff_sfc = []
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                ilev=0,
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                normalize_by_area=normalize_by_area,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_sfc,
                weightsdir=weightsdir,
                second_ref=secondrefds,
                second_dev=seconddevds
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
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                normalize_by_area=normalize_by_area,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_500,
                weightsdir=weightsdir,
                second_ref=secondrefds,
                second_dev=seconddevds
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
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                log_color_scale=log_color_scale,
                normalize_by_area=normalize_by_area,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_zm,
                weightsdir=weightsdir,
                second_ref=secondrefds,
                second_dev=seconddevds
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
                pdfname=pdfname,
                use_cmap_RdBu=use_cmap_RdBu,
                pres_range=[1, 100],
                log_yaxis=True,
                extra_title_txt=extra_title_txt,
                log_color_scale=log_color_scale,
                normalize_by_area=normalize_by_area,
                weightsdir=weightsdir,
                second_ref=secondrefds,
                second_dev=seconddevds
            )
            util.add_nested_bookmarks_to_pdf(
                pdfname, filecat, catdict,
                warninglist, remove_prefix=coll_prefix
            )
            return {filecat : cat_diff_dict}
    # Create the plots in parallel
    results = Parallel(n_jobs=n_job)(
        delayed(createplots)(i, filecat) for i, filecat in enumerate(catdict)
    )

    dict_sfc = {list(result.keys())[0] : result[list(result.keys())[0]]['sfc'] for result in results}
    dict_500 = {list(result.keys())[0] : result[list(result.keys())[0]]['500'] for result in results}
    dict_zm  = {list(result.keys())[0] : result[list(result.keys())[0]]['zm']  for result in results}

    # ==============================================================
    # Write the list of species having significant differences,
    # which we need to fill out the benchmark approval forms.
    # ==============================================================
    if sigdiff_files != None:
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


def make_benchmark_emis_plots(
    ref,
    refstr,
    dev,
    devstr,
    dst="./1mo_benchmark",
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
):
    """
    Creates PDF files containing plots of emissions for model
    benchmarking purposes. This function is compatible with benchmark
    simulation output only. It is not compatible with transport tracers
    emissions diagnostics.

    Args:
    -----
        ref: str
            Path name for the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where
            PDF files containing plots will be written.
            Default value: './1mo_benchmark

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        plot_by_spc_cat : boolean
            Set this flag to True to separate plots into PDF files
            according to the benchmark species categories (e.g. Oxidants,
            Aerosols, Nitrogen, etc.)  These categories are specified
            in the YAML file benchmark_species.yml.
            Default value: False

        plot_by_hco_cat : boolean
            Set this flag to True to separate plots into PDF files
            according to HEMCO emissions categories (e.g. Anthro,
            Aircraft, Bioburn, etc.)
            Default value: False

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False

        flip_ref : boolean
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False

        flip_dev : boolean
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False

        log_color_scale: boolean
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False

         sigdiff_files : list of str
            Filenames that will contain the lists of species having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None

        weightsdir : str
            Directory in which to place (and possibly reuse) xESMF regridder netCDF files.
            Default value: '.'

        n_job : int
            Defines the number of simultaneous workers for parallel plotting.
            Set to 1 to disable parallel plotting. Value of -1 allows the application to decide.
            Default value: -1

    Remarks:
    --------
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

    # Ref dataset
    try:
        refds = xr.open_dataset(ref, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Dev dataset
    try:
        devds = xr.open_dataset(dev, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = util.add_missing_variables(refds, devds)

    # Create regridding files if necessary while not in parallel loop
    [ _ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

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
            weightsdir=weightsdir
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
                varnames = [k for k in emis_vars \
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
                pdfname = os.path.join(emisspcdir, "{}_Emissions.pdf".format(c))
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
                weightsdir=weightsdir
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

        results = Parallel(n_jobs=n_job)(delayed(createfile_hco_cat)(c) \
                                         for c in emis_cats)

        dict_emis = {list(result.keys())[0] : result[list(result.keys())[0]] \
                     for result in results}

        # =============================================================
        # Write the list of species having significant differences,
        # which we need to fill out the benchmark approval forms.
        # =============================================================
        if sigdiff_files != None:
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
        warninglist = (
            []
        )  # in case any emissions are skipped (for use in nested pdf bookmarks)
        allcatspc = (
            []
        )  # for checking if emissions species not defined in benchmark category file
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
                        emisvars = [v for v in emis_vars \
                                    if spc == v.split("_")[0][4:]]
                        for var in emisvars:
                            emisdict[filecat][spc].append(var.replace("Emis", ""))
                            varlist.append(var)
            if not varlist:
                print(
                    "\nWarning: no emissions species in benchmark species category {}".format(
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
            print(pdfname)
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
                weightsdir=weightsdir
            )
            util.add_nested_bookmarks_to_pdf(pdfname, filecat, emisdict, warninglist)
            return catspc
        results = Parallel(n_jobs=n_job)(
            delayed(createfile_bench_cat)(filecat) \
            for i, filecat in enumerate(catdict)
        )

        allcatspc = [spc for result in results for spc in result]
        # Give warning if emissions species is not assigned a benchmark category
        for spc in emis_spc:
            if spc not in allcatspc:
                print(
                    "Warning: species {} has emissions diagnostics but is not in benchmark_categories.yml".format(
                        spc
                    )
                )


def make_benchmark_emis_tables(
    reflist,
    refstr,
    devlist,
    devstr,
    dst="./1mo_benchmark",
    overwrite=False,
    interval=[2678400.0],
):
    """
    Creates a text file containing emission totals by species and
    category for benchmarking purposes.

    Args:
    -----
        reflist: list of str
             List with the path names of the emissions file, or emissions
             and met field files, that will constitute the "Ref" (aka
             "Reference") data set. If two files are passed in the list,
             the met field file must be second.

        refstr : str
            A string to describe ref (e.g. version number)

        devlist : list of str
             List with the path names of the emissions file, or emissions
             and met field files, that will constitute the "Dev" (aka
             "Development") data set. If two files are passed in the list,
             the met field file must be second. The "Dev" data set will
             be compared against the "Ref" data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./1mo_benchmark

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value : False

        interval : list of float
            The length of the data interval in seconds. By default, interval
            is set to [2678400.0], which is the number of seconds in July
            (our 1-month benchmarking month).
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
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Create the "Emissions" category folder if it does not exist
    emisdir = os.path.join(dst, "Emissions")
    if not os.path.isdir(emisdir):
        os.mkdir(emisdir)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Read the Ref dataset and make sure that the area variables are present
    if len(reflist) == 1:
        reflist = [reflist]
    refds = xr.open_mfdataset(reflist, drop_variables=gcon.skip_these_vars)
    refds = util.check_for_area(refds)

    # Read the Dev dataset and make sure that area variables are present
    if len(devlist) == 1:
        devlist = [devlist]
    devds = xr.open_mfdataset(devlist, drop_variables=gcon.skip_these_vars)
    devds = util.check_for_area(devds)

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
        interval,
        template="Emis{}_",
    )

    # Create table of emissions by inventory
    create_total_emissions_table(
        refds,
        refstr,
        devds,
        devstr,
        inventories,
        file_inv_totals,
        interval,
        template="Inv{}_",
    )


def make_benchmark_jvalue_plots(
    ref,
    refstr,
    dev,
    devstr,
    varlist=None,
    dst="./1mo_benchmark",
    subdst=None,
    local_noon_jvalues=False,
    plots=["sfc", "500hpa", "zonalmean"],
    overwrite=False,
    verbose=False,
    flip_ref=False,
    flip_dev=False,
    log_color_scale=False,
    sigdiff_files=None,
    weightsdir='.'
):
    """
    Creates PDF files containing plots of J-values for model
    benchmarking purposes.

    Args:
    -----
        ref: str
            Path name for the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        varlist : list of str
            List of J-value variables to plot.  If not passed,
            then all J-value variables common to both dev
            and ref will be plotted.  The varlist argument can be
            a useful way of restricting the number of variables
            plotted to the pdf file when debugging.
            Default value: None

        dst : str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./1mo_benchmark.

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        local_noon_jvalues : boolean
            Set this flag to plot local noon J-values.  This will
            divide all J-value variables by the JNoonFrac counter,
            which is the fraction of the time that it was local noon
            at each location.
            Default value : False

        plots : list of strings
            List of plot types to create.
            Default value: ['sfc', '500hpa', 'zonalmean']

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False

        flip_ref : boolean
            Set this flag to True to reverse the vertical level
            ordering in the "Ref" dataset (in case "Ref" starts
            from the top of atmosphere instead of the surface).
            Default value: False

        flip_dev : boolean
            Set this flag to True to reverse the vertical level
            ordering in the "Dev" dataset (in case "Dev" starts
            from the top of atmosphere instead of the surface).
            Default value: False

        log_color_scale: boolean
            Set this flag to True if you wish to enable plotting data
            (not diffs) on a log color scale.
            Default value: False

        sigdiff_files : list of str
            Filenames that will contain the lists of J-values having
            significant differences in the 'sfc', '500hpa', and
            'zonalmean' plots.  These lists are needed in order to
            fill out the benchmark approval forms.
            Default value: None

        weightsdir : str
            Directory in which to place (and possibly reuse) xESMF regridder netCDF files.
            Default value: '.'

    Remarks:
    --------
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
    elif not os.path.isdir(dst):
        os.mkdir(dst)

    # Ref dataset
    try:
        refds = xr.open_dataset(ref, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Dev dataset
    try:
        devds = xr.open_dataset(dev, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = util.add_missing_variables(refds, devds)

    # Create regridding files if necessary
    [ _ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # Get a list of the 3D variables in both datasets
    if varlist == None:
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
        if varlist == None:
            varlist = [v for v in cmn if prefix in v]

        # Make sure JNoonFrac (fraction of times it was local noon
        # in each column) is present in both Ref and Dev datasets
        if not "JNoonFrac" in cmn:
            msg = "JNoonFrac is not common to Ref and Dev datasets!"
            raise ValueError(msg)

        # JNoon_* are cumulative sums of local noon J-values; we need
        # to divide these by JNoonFrac to get the average value
        refds = util.divide_dataset_by_dataarray(refds, \
                                                 refds["JNoonFrac"], varlist)
        devds = util.divide_dataset_by_dataarray(devds, \
                                                 devds["JNoonFrac"], varlist)

        # Subfolder of dst where PDF files will be printed
        catdir = "JValuesLocalNoon"

    else:

        # Get a list of continuously averaged J-value variables
        # (or use the varlist passed via tha argument list)
        prefix = "Jval_"
        if varlist == None:
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
            pdfname = os.path.join(jvdir, "{}Surface_{}.pdf".format(
                prefix, subdst))
        else:
            pdfname = os.path.join(jvdir, "{}Surface.pdf".format(prefix))

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
            weightsdir=weightsdir
        )
        diff_sfc[:] = [v.replace(prefix, "") for v in diff_sfc]
        util.add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=prefix, verbose=verbose)

    # 500hPa plots
    if "500hpa" in plots:
        if subdst is not None:
            pdfname = os.path.join(jvdir, "{}500hPa_{}.pdf".format(
                prefix, subdst))
        else:
            pdfname = os.path.join(jvdir, "{}500hPa.pdf".format(prefix))

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
            weightsdir=weightsdir
        )
        diff_500[:] = [v.replace(prefix, "") for v in diff_500]
        util.add_bookmarks_to_pdf(pdfname, varlist,
                             remove_prefix=prefix, verbose=verbose)
    # Full-column zonal mean plots
    if "zonalmean" in plots:
        if subdst is not None:
            pdfname = os.path.join(
                jvdir, "{}FullColumn_ZonalMean_{}.pdf".format(prefix, subdst)
            )
        else:
            pdfname = os.path.join(jvdir, "{}FullColumn_ZonalMean.pdf".format(
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
            weightsdir=weightsdir
        )
        diff_zm[:] = [v.replace(prefix, "") for v in diff_zm]
        util.add_bookmarks_to_pdf(pdfname, varlist, remove_prefix=prefix, verbose=verbose)

        # Strat_ZonalMean plots will use a log-pressure Y-axis, with
        # a range of 1..100 hPa, as per GCSC request. (bmy, 8/13/19)
        if subdst is not None:
            pdfname = os.path.join(
                jvdir, "{}Strat_ZonalMean_{}.pdf".format(prefix, subdst)
            )
        else:
            pdfname = os.path.join(jvdir, "{}Strat_ZonalMean.pdf".format(prefix))

        compare_zonal_mean(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            pdfname=pdfname,
            pres_range=[0, 100],
            log_yaxis=True,
            flip_ref=flip_ref,
            flip_dev=flip_dev,
            extra_title_txt=extra_title_txt,
            log_color_scale=log_color_scale,
            weightsdir=weightsdir
        )
        util.add_bookmarks_to_pdf(pdfname, varlist,
                             remove_prefix=prefix, verbose=verbose)

        # ==============================================================
        # Write the lists of J-values that have significant differences,
        # which we need to fill out the benchmark approval forms.
        # ==============================================================
        if sigdiff_files != None:
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


def make_benchmark_aod_plots(
    ref,
    refstr,
    dev,
    devstr,
    varlist=None,
    dst="./1mo_benchmark",
    subdst=None,
    overwrite=False,
    verbose=False,
    log_color_scale=False,
    sigdiff_files=None,
    weightsdir='.'
):
    """
    Creates PDF files containing plots of column aerosol optical
    depths (AODs) for model benchmarking purposes.

    Args:
    -----
        ref: str
            Path name for the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : str
            Path name for the "Dev" (aka "Development") data set.
            This data set will be compared against the "Reference"
            data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        varlist : list of str
            List of AOD variables to plot.  If not passed, then all
            AOD variables common to both Dev and Ref will be plotted.
            Use the varlist argument to restrict the number of
            variables plotted to the pdf file when debugging.
            Default value: None

        dst : str
            A string denoting the destination folder where a
            PDF file  containing plots will be written.
            Default value: ./1mo_benchmark.

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False.

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False

        log_color_scale: boolean
            Set this flag to True to enable plotting data (not diffs)
            on a log color scale.
            Default value: False

        sigdiff_files : list of str
            Filenames that will contain the list of quantities having
            having significant differences in the column AOD plots.
            These lists are needed in order to fill out the benchmark
            approval forms.
            Default value: None

        weightsdir : str
            Directory in which to place (and possibly reuse) xESMF regridder netCDF files.
            Default value: '.'

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
    elif not os.path.isdir(dst):
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

    # Read the Ref dataset
    try:
        refds = xr.open_dataset(ref, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Ref file: {}".format(ref))

    # Read the Dev dataset
    try:
        devds = xr.open_dataset(dev, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Could not find Dev file: {}".format(dev))

    # Create regridding files if necessary
    [ _ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

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
    if varlist == None:
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
        weightsdir=weightsdir
    )
    diff_aod[:] = [v.replace("Column_AOD_", "") for v in diff_aod]
    util.add_bookmarks_to_pdf(
        pdfname, newvarlist, remove_prefix="Column_AOD_", verbose=verbose
    )

    # ==================================================================
    # Write the list of AOD quantities having significant differences,
    # which we will need to fill out the benchmark forms.
    # ==================================================================
    if sigdiff_files != None:
        for filename in sigdiff_files:
            if "sfc" in filename:
                with open(filename, "a+") as f:
                    print("* Column AOD: ", file=f, end="")
                    for v in diff_aod:
                        print("{} ".format(v), file=f, end="")
                    print(file=f)
                    f.close()


def make_benchmark_mass_tables(
    reflist,
    refstr,
    devlist,
    devstr,
    varlist=None,
    dst="./1mo_benchmark",
    subdst=None,
    overwrite=False,
    verbose=False,
    label="at end of simulation",
):
    """
    Creates a text file containing global mass totals by species and
    category for benchmarking purposes.

    Args:
    -----
        reflist : list of str
            List of files (i.e. pathnames) that will constitute
            the "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        dev : list of str
            List of files (i.e. pathnames) that will constitute
            the "Dev" (aka "Development") data set.  The "Dev"
            data set will be compared against the "Ref" data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        varlist : list of str
            List of variables to include in the list of totals.
            If omitted, then all variables that are found in either
            "Ref" or "Dev" will be included.  The varlist argument
            can be a useful way of reducing the number of
            variables during debugging and testing.
            Default value: None

        dst : str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./1mo_benchmark

        subdst : str
            A string denoting the sub-directory of dst where PDF
            files containing plots will be written.  In practice,
            subdst is only needed for the 1-year benchmark output,
            and denotes a date string (such as "Jan2016") that
            corresponds to the month that is being plotted.
            Default value: None

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value : False

        verbose : boolean
            Set this flag to True to print extra informational output.
            Default value: False.
    """

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    elif not os.path.isdir(dst):
        os.makedirs(dst)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Ref
    try:
        refds = xr.open_mfdataset(reflist, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Error opening Ref files: {}".format(reflist))

    # Dev dataset
    try:
        devds = xr.open_mfdataset(devlist, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError("Error opening Dev files: {}!".format(devlist))

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
    ref_area = util.get_area_from_dataset(refds)
    dev_area = util.get_area_from_dataset(devds)

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
                msg = msg.format(dst)
                raise ValueError('Variable {} in varlist passed to make_benchmark_mass_tables is not present in ref and dev datasets'.format(v))
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
    )


def make_benchmark_oh_metrics(
    reflist,
    refstr,
    devlist,
    devstr,
    dst="./1mo_benchmark",
    overwrite=False,
):
    """
    Creates a text file containing metrics of global mean OH, MCF lifetime,
    and CH4 lifetime for benchmarking purposes.

    Args:
    -----
        reflist: list of str
            List with the path names of files that will constitute the
            "Ref" (aka "Reference") data set.

        refstr : str
            A string to describe ref (e.g. version number)

        devlist : list of str
            List with the path names of files that will constitute the
            "Dev" (aka "Development") data set.  The "Dev" data set will be
            compared against the "Ref" data set.

        devstr : str
            A string to describe dev (e.g. version number)

    Keyword Args (optional):
    ------------------------
        dst : str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./1mo_benchmark

        overwrite : boolean
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value : False
    """

    # ==================================================================
    # Define destination directory
    # ==================================================================
    if os.path.isdir(dst) and not overwrite:
        msg = "Directory {} exists. Pass overwrite=True to overwrite " \
            + "files in that directory, if any."
        msg = msg.format(dst)
        raise ValueError(msg)
    elif not os.path.isdir(dst):
        os.makedirs(dst)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Ref
    try:
        refds = xr.open_mfdataset(reflist, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError(
            "Could not find one of the Ref files: {}".format(reflist)
        )

    # Dev
    try:
        devds = xr.open_mfdataset(devlist, drop_variables=gcon.skip_these_vars)
    except FileNotFoundError:
        raise FileNotFoundError(
            "Could not find one of the Dev files: {}".format(devlist)
        )

    # Make sure that required variables are found
    if "OHconcAfterChem" not in refds.data_vars.keys():
        raise ValueError('Could not find "OHconcAfterChem" in Ref!')
    if "OHconcAfterChem" not in devds.data_vars.keys():
        raise ValueError('Could not find "OHconcAfterChem" in Dev!')

    # ==================================================================
    # Make sure that all necessary variables are found
    # ==================================================================

    # Find the area variables in Ref and Dev
    ref_area = util.get_area_from_dataset(refds)
    dev_area = util.get_area_from_dataset(devds)

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
    refmet = util.get_variables_from_dataset(refds, metvar_list)
    devmet = util.get_variables_from_dataset(devds, metvar_list)

    # Create the mask arrays for the troposphere for Ref and Dev
    ref_tropmask = get_troposphere_mask(refmet)
    dev_tropmask = get_troposphere_mask(devmet)

    # Get the OH concentration
    ref_oh = refds["OHconcAfterChem"]
    dev_oh = devds["OHconcAfterChem"]

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Create file
    outfilename = os.path.join(dst, "{}_OH_metrics.txt".format(devstr))
    try:
        f = open(outfilename, "w")
    except FileNotFoundError:
        raise FileNotFoundError("Could not open {} for writing!".format(
            outfilename))

    # ==================================================================
    # Compute mass-weighted OH in the troposphere
    # ==================================================================

    # Physical constants
    Avo = gcon.AVOGADRO   # molec/mol
    mw_air = gcon.MW_AIR  # g/mole dry air
    g0 = gcon.G           # m/s2

    # Ref
    ref_oh_trop = np.ma.masked_array(ref_oh.values, ref_tropmask)
    ref_airmass_trop = np.ma.masked_array(refmet["Met_AD"].values, ref_tropmask)
    ref_oh_mass = ref_oh_trop * ref_airmass_trop
    ref_total_ohmass = np.sum(ref_oh_mass)
    ref_total_airmass = np.sum(ref_airmass_trop)
    ref_mean_oh = (ref_total_ohmass / ref_total_airmass) / 1e5

    # Dev
    dev_oh_trop = np.ma.masked_array(dev_oh.values, dev_tropmask)
    dev_airmass_trop = np.ma.masked_array(devmet["Met_AD"].values, dev_tropmask)
    dev_oh_mass = dev_oh_trop * dev_airmass_trop
    dev_total_ohmass = np.sum(dev_oh_mass)
    dev_total_airmass = np.sum(dev_airmass_trop)
    dev_mean_oh = (dev_total_ohmass / dev_total_airmass) / 1e5

    oh_diff = dev_mean_oh - ref_mean_oh
    oh_pctdiff = ((dev_mean_oh - ref_mean_oh) / ref_mean_oh) * 100.0

    # Title strings
    title1 = "### Global mass-weighted OH concentration [1e5 molec/cm3]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

    # Print header to file
    print("#" * 79, file=f)
    print("{}{}".format(title1.ljust(76), "###"), file=f)
    print("{}{}".format(title2.ljust(76), "###"), file=f)
    print("#" * 79, file=f)

    # Write results to file
    print(
        "{}{}{}{}".format(
            "  Ref".ljust(15),
            "Dev".ljust(13),
            "Dev - Ref".ljust(13),
            "% diff".ljust(11),
        ),
        file=f,
    )
    print(
        "{:11.6f}  {:11.6f}  {:11.6f}  {:9.4f}".format(
            ref_mean_oh, dev_mean_oh, oh_diff, oh_pctdiff
        ),
        file=f,
    )

    # ==================================================================
    # Compute MCF and CH4 lifetimes
    # ==================================================================

    # Select only boxes that are purely tropospheric
    # This excludes influence from the stratosphere
    ref_timetrop_mask = refmet["FracOfTimeInTrop"].values != 1.0
    dev_timetrop_mask = devmet["FracOfTimeInTrop"].values != 1.0

    # Get grid box volumes [cm3] (trop + strat)
    ref_vol = (refmet["Met_BXHEIGHT"] * ref_area) * 1e6
    dev_vol = (devmet["Met_BXHEIGHT"] * dev_area) * 1e6

    # Get grid box volumes [cm3] (trop only)
    ref_vol_trop = np.ma.masked_array(ref_vol.values, ref_timetrop_mask)
    dev_vol_trop = np.ma.masked_array(dev_vol.values, dev_timetrop_mask)

    # Get MCF and CH4 density [molec/cm3] (trop + strat)
    # Assume that species is evenly distributed in air, with
    # a mixing ratio of 1. Thus species density = air density.
    ref_dens = refmet["Met_AIRDEN"] / 1e6
    dev_dens = devmet["Met_AIRDEN"] / 1e6

    # Get MCF and CH4 density [molec/cm3] (trop only)
    ref_dens_trop = np.ma.masked_array(ref_dens.values, ref_timetrop_mask)
    dev_dens_trop = np.ma.masked_array(dev_dens.values, dev_timetrop_mask)

    # Get temperature [K] (trop only)
    ref_temp = np.ma.masked_array(refmet["Met_T"].values, ref_timetrop_mask)
    dev_temp = np.ma.masked_array(devmet["Met_T"].values, dev_timetrop_mask)

    # Compute Arrhenius parameter K [cm3/molec/s]
    ref_mcf_k = 1.64e-12 * np.exp(-1520e0 / ref_temp)
    dev_mcf_k = 1.64e-12 * np.exp(-1520e0 / dev_temp)
    ref_ch4_k = 2.45e-12 * np.exp(-1775e0 / ref_temp)
    dev_ch4_k = 2.45e-12 * np.exp(-1775e0 / dev_temp)

    # Numerator: Total atmospheric (trop+strat) burden
    ref_num = np.sum(ref_dens.values * ref_vol.values)
    dev_num = np.sum(dev_dens.values * dev_vol.values)

    # Denominator: Loss rate in troposphere
    ref_mcf_denom = np.sum(ref_mcf_k * ref_oh_trop * \
                           ref_dens_trop * ref_vol_trop)

    dev_mcf_denom = np.sum(dev_mcf_k * dev_oh_trop * dev_dens_trop * dev_vol_trop)
    ref_ch4_denom = np.sum(ref_ch4_k * ref_oh_trop * ref_dens_trop * ref_vol_trop)
    dev_ch4_denom = np.sum(dev_ch4_k * dev_oh_trop * dev_dens_trop * dev_vol_trop)

    # Compute lifetimes [years]
    sec_to_year = 365.25 * 86400.0
    ref_mcf_lifetime = (ref_num / ref_mcf_denom) / sec_to_year
    dev_mcf_lifetime = (dev_num / dev_mcf_denom) / sec_to_year
    ref_ch4_lifetime = (ref_num / ref_ch4_denom) / sec_to_year
    dev_ch4_lifetime = (dev_num / dev_ch4_denom) / sec_to_year

    # Compute differences
    mcf_diff = dev_mcf_lifetime - ref_mcf_lifetime
    ch4_diff = dev_ch4_lifetime - ref_ch4_lifetime

    mcf_pctdiff = ((dev_mcf_lifetime - ref_mcf_lifetime) / ref_mcf_lifetime) * 100.0
    ch4_pctdiff = ((dev_ch4_lifetime - ref_ch4_lifetime) / ref_ch4_lifetime) * 100.0

    # Title strings
    title1 = "### MCF lifetime w/r/t tropospheric OH [years]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

    # Print header to file
    print("", file=f)
    print("#" * 79, file=f)
    print("{}{}".format(title1.ljust(76), "###"), file=f)
    print("{}{}".format(title2.ljust(76), "###"), file=f)
    print("#" * 79, file=f)

    # Write results to file
    print(
        "{}{}{}{}".format(
            "  Ref".ljust(15),
            "Dev".ljust(13),
            "Dev - Ref".ljust(13),
            "% diff".ljust(11),
        ),
        file=f,
    )
    print(
        "{:11.6f}  {:11.6f}  {:11.6f}  {:9.4f}".format(
            ref_mcf_lifetime, dev_mcf_lifetime, mcf_diff, mcf_pctdiff
        ),
        file=f,
    )

    # Title strings
    title1 = "### CH4 lifetime w/r/t tropospheric OH [years]"
    title2 = "### Ref = {}; Dev = {}".format(refstr, devstr)

    # Print header to file
    print("", file=f)
    print("#" * 79, file=f)
    print("{}{}".format(title1.ljust(76), "###"), file=f)
    print("{}{}".format(title2.ljust(76), "###"), file=f)
    print("#" * 79, file=f)

    # Write results to file
    print(
        "{}{}{}{}".format(
            "  Ref".ljust(15),
            "Dev".ljust(13),
            "Dev - Ref".ljust(13),
            "% diff".ljust(11),
        ),
        file=f,
    )
    print(
        "{:11.6f}  {:11.6f}  {:11.6f}  {:9.4f}".format(
            ref_ch4_lifetime, dev_ch4_lifetime, ch4_diff, ch4_pctdiff
        ),
        file=f,
    )
