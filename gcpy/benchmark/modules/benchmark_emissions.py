#!/usr/bin/env python3
"""
Creates emissions plots & tables from GEOS-Chem benchmark simulation output.
"""
import os
import gc
import numpy as np
import xarray as xr
from joblib import Parallel, delayed
from gcpy.regrid import create_regridders
from gcpy.util import \
    add_bookmarks_to_pdf, add_missing_variables, add_nested_bookmarks_to_pdf, \
    compare_varnames, create_blank_dataarray, dataset_reader, \
    dataset_mean, get_emissions_varnames, insert_text_into_file, \
    make_directory, print_totals, read_config_file, read_species_metadata, \
    replace_whitespace, verify_variable_type
from gcpy.units import convert_units
from gcpy.constants import \
    COL_WIDTH, ENCODING, SKIP_THESE_VARS, TABLE_WIDTH
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.benchmark.modules.benchmark_utils import \
    EMISSION_SPC, EMISSION_INV, get_species_categories

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def create_total_emissions_table(
        refdata,
        refstr,
        devdata,
        devstr,
        species,
        spcdb_files,
        outfilename,
        ref_interval=[2678400.0],
        dev_interval=[2678400.0],
        template="Emis{}_",
        refmetdata=None,
        devmetdata=None,
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
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs
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

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species metadata (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Initialization
    # ==================================================================
    verify_variable_type(refdata, xr.Dataset)
    verify_variable_type(refstr, str)
    verify_variable_type(devdata, xr.Dataset)
    verify_variable_type(devstr, str)
    verify_variable_type(species, dict)
    verify_variable_type(spcdb_files, list)
    verify_variable_type(outfilename, str)

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

    # Read the species database files in the Ref & Dev rundirs,
    # and return a dict for each containing the species metadata.
    ref_metadata, dev_metadata = read_species_metadata(
        spcdb_files,
        quiet=True
    )

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Get the list of emission variables for which we will print totals
    # ==================================================================

    # Find all common variables between the two datasets
    # and get the lists of variables only in Ref and only in Dev,
    vardict = compare_varnames(refdata, devdata, quiet=True)
    cvars = vardict["commonvars"]
    refonly = vardict["refonly"]
    devonly = vardict["devonly"]

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refdata, devdata] = add_missing_variables(refdata, devdata)

    # =================================================================
    # Open the file for output
    # =================================================================
    # Create file
    try:
        f = open(outfilename, "w", encoding=ENCODING)
    except (IOError, OSError, FileNotFoundError) as e:
        raise e(f"Could not open {outfilename} for writing!") from e

    # Write a placeholder to the file that denotes where
    # the list of species with differences will be written
    placeholder = "@%% insert diff status here %%@"
    print(f"{placeholder}\n\n", file=f)

    # Define a list for differences
    diff_list = []

    # =================================================================
    # Loop through all of the species are in species_dict
    # =================================================================
    for species_name, target_units in species.items():

        # Get a list of emission variable names for each species
        diagnostic_template = template.replace("{}", species_name)
        varnames = get_emissions_varnames(cvars, diagnostic_template)

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
            msg = f"No emissions found for {species_name} ... skippping"
            print(msg)
            continue

        # Check if there is a total emissions variable in the list
        vartot = [v for v in varnames if "_TOTAL" in v.upper()]

        # Push the total variable to the last list element
        # so that it will be printed last of all
        if len(vartot) == 1:
            varnames.append(varnames.pop(varnames.index(vartot[0])))

        # Title strings
        if "Inv" in template:
            print(f"Computing inventory totals for {species_name}")
            title0 = f"for inventory {species_name}"
            title1 = f"### Emissions totals {title0} [Tg]"
        else:
            print(f"Computing emissions totals for {species_name}")
            title0 = f"for species {species_name}"
            title1 = f"### Emissions totals {title0} [Tg]"

        title2 = f"### Ref = {refstr}"
        title3 = f"### Dev = {devstr}"

        # Print header to file
        print("#" * TABLE_WIDTH, file=f)
        print(f"{title1 : <{TABLE_WIDTH-3}}{'###'}", file=f)
        print(f"{title2 : <{TABLE_WIDTH-3}}{'###'}", file=f)
        print(f"{title3 : <{TABLE_WIDTH-3}}{'###'}", file=f)
        print("#" * TABLE_WIDTH, file=f)
        print(f"{'' : <{COL_WIDTH-1}}{'Ref' : >{COL_WIDTH}}{'Dev' : >{COL_WIDTH}}{'Dev - Ref' : >{COL_WIDTH}}{'% diff' : >{COL_WIDTH}} {'diffs'}", file=f)

        # =============================================================
        # Loop over all emissions variables corresponding to this
        # species and print their totals in Ref and Dev to the file.
        # =============================================================
        for v in varnames:

            # KLUDGE, skip InvAFCID due to a file error in GCHP
            if "InvAFCID" in v:
                continue

            if "Inv" in template:
                spc_name = v.split("_")[1]
            else:
                spc_name = species_name

            # Get metadata for the given species
            ref_species_metadata = ref_metadata.get(spc_name)
            dev_species_metadata = dev_metadata.get(spc_name)
            if ref_species_metadata is None and dev_species_metadata is None:
                print(f"No metadata found for {spc_name} ... skipping")
                continue

            # Convert units of Ref and Dev and save to numpy ndarray objects
            # (or set to NaN if the variable is not found in Ref or Dev)
            if v in refonly and v not in devonly:

                # Convert units of Ref
                refarray = convert_units(
                    refdata[v],
                    spc_name,
                    ref_species_metadata,
                    target_units,
                    interval=ref_interval,
                    area_m2=refarea,
                )

                # Set Dev to NaN (missing values) everywhere
                devarray = create_blank_dataarray(
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
                    dev_species_metadata,
                    target_units,
                    interval=dev_interval,
                    area_m2=devarea,
                )

                # Set Ref to NaN (missing values) everywhere
                refarray = create_blank_dataarray(
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
                    ref_species_metadata,
                    target_units,
                    interval=ref_interval,
                    area_m2=refarea,
                )
                devarray = convert_units(
                    devdata[v],
                    spc_name,
                    dev_species_metadata,
                    target_units,
                    interval=dev_interval,
                    area_m2=devarea,
                )

            # ==========================================================
            # Print emission totals for Ref and Dev
            # ==========================================================
            print_totals(
                refarray,
                devarray,
                f,
                diff_list
            )

        # Add newlines before going to the next species
        print(file=f)
        print(file=f)

    # =================================================================
    # Cleanup and quit
    # =================================================================

    # Close file
    f.close()

    # Reopen file and replace placeholder with list of diffs
    insert_text_into_file(
        filename=outfilename,
        search_text=placeholder,
        replace_text=diff_list_to_text(
            refstr,
            devstr,
            diff_list),
        width=TABLE_WIDTH
    )


def make_benchmark_emis_plots(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files,
        dst="./benchmark",
        subdst=None,
        plot_by_spc_cat=False,
        plot_by_hco_cat=False,
        benchmark_type="FullChemBenchmark",
        cmpres=None,
        overwrite=False,
        verbose=False,
        flip_ref=False,
        flip_dev=False,
        log_color_scale=False,
        sigdiff_files=None,
        weightsdir='.',
        n_job=-1,
        time_mean=False,
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
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

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
            in the YAML file benchmark_categories.yml.
            Default value: False
        plot_by_hco_cat: bool
            Set this flag to True to separate plots into PDF files
            according to HEMCO emissions categories (e.g. Anthro,
            Aircraft, Bioburn, etc.)
            Default value: False
        benchmark_type: str
            A string denoting the type of benchmark output to plot, options are
            FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
            Default value: "FullChemBenchmark"
        cmpres: string
            Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
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
    verify_variable_type(ref, (str,list))
    verify_variable_type(refstr, str)
    verify_variable_type(dev, (str,list))
    verify_variable_type(devstr, str)
    verify_variable_type(spcdb_files, list)

    # Create the destination folder
    make_directory(dst, overwrite)

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

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Get the function that will read the dataset
    reader = dataset_reader(time_mean, verbose=verbose)

    # Ref dataset
    try:
        refds = reader(ref, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {ref}") from e

    # Dev dataset
    try:
        devds = reader(dev, drop_variables=SKIP_THESE_VARS)
    except (OSError, IOError, FileNotFoundError) as e:
        raise e(f"Could not find Ref file: {dev}") from e

    # Compute mean of data over the time dimension (if time_mean=True)
    if time_mean:
        refds = dataset_mean(refds)
        devds = dataset_mean(devds)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    # Create regridding files if necessary while not in parallel loop
    [_ for _ in create_regridders(refds, devds, weightsdir=weightsdir)]

    # Combine 2D and 3D variables into an overall list
    vardict = compare_varnames(refds, devds, quiet=not verbose)
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
            pdfname = os.path.join(
                emisdir,
                f"Emissions_{subdst}.pdf"
            )
        else:
            pdfname = os.path.join(
                emisdir,
                "Emissions.pdf"
            )

        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            cmpres=cmpres,
            pdfname=pdfname,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            weightsdir=weightsdir,
            n_job=n_job,
            spcdb_files=spcdb_files,
        )
        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix="Emis",
            verbose=verbose
        )
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
                    emisspcdir,
                    f"{c}_Emissions_{subdst}.pdf"
                )
            else:
                pdfname = os.path.join(
                    emisspcdir,
                    f"{c}_Emissions.pdf"
                )
            diff_dict = {}
            diff_emis = []
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varnames,
                cmpres=cmpres,
                ilev=0,
                pdfname=pdfname,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
                sigdiff_list=diff_emis,
                weightsdir=weightsdir,
                n_job=n_job,
                spcdb_files=spcdb_files
            )

            add_bookmarks_to_pdf(
                pdfname,
                varnames,
                remove_prefix="Emis",
                verbose=verbose
            )
            # Save the list of quantities with significant differences for
            # this category into the diff_dict dictionary for use below
            diff_emis[:] = [v.replace("Emis", "") for v in diff_emis]
            diff_emis[:] = [v.replace("_" + c, "") for v in diff_emis]
            diff_dict[c] = diff_emis
            return diff_dict

        # ---------------------------------------
        # Create plots in parallel
        # Turn off parallelization if n_job=1
        if n_job != 1:
            results = Parallel(n_jobs=n_job)(
                delayed(createfile_hco_cat)(c)
                for c in emis_cats
            )
        else:
            results = []
            for c in emis_cats:
                results.append(createfile_hco_cat(c))
        # ---------------------------------------

        dict_emis = {list(result.keys())[0]: result[list(result.keys())[0]]
                     for result in results}

        # =============================================================
        # Write the list of species having significant differences,
        # which we need to fill out the benchmark approval forms.
        # =============================================================
        if sigdiff_files is not None:
            for filename in sigdiff_files:
                if "emis" in filename:
                    with open(filename, "w+", encoding=ENCODING) as f:
                        for c, diff_list in dict_emis.items():
                            print(f"* {c}: ", file=f, end="")
                            for v in diff_list:
                                print(f"{v} ", file=f, end="")
                            print(file=f)
                        f.close()

    # ==================================================================
    # if plot_by_spc_cat is true, make a file for each benchmark
    # species category with emissions in the diagnostics file
    # ==================================================================
    if plot_by_spc_cat:

        catdict = get_species_categories(benchmark_type)
        # in case any emissions are skipped (for use in nested pdf bookmarks)
        warninglist = []
        # for checking if emissions species not defined in benchmark category
        # file
        allcatspc = []
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
                    f"category {filecat}"
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
                    catdir,
                    f"{filecat}_Emissions_{subdst}.pdf"
                )
            else:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_Emissions.pdf"
                )
            # Create the PDF
            compare_single_level(
                refds,
                refstr,
                devds,
                devstr,
                varlist=varlist,
                cmpres=cmpres,
                ilev=0,
                pdfname=pdfname,
                flip_ref=flip_ref,
                flip_dev=flip_dev,
                log_color_scale=log_color_scale,
                extra_title_txt=extra_title_txt,
                weightsdir=weightsdir,
                n_job=n_job,
                spcdb_files=spcdb_files,
            )
            add_nested_bookmarks_to_pdf(
                pdfname,
                filecat,
                emisdict,
                warninglist
            )
            return catspc

        #------------------------------------------------
        # Create plots in parallel
        # Turn of parallalization if n_job=1
        if n_job != 1:
            results = Parallel(n_jobs=n_job)(
                delayed(createfile_bench_cat)(filecat)
                for _, filecat in enumerate(catdict)
            )
        else:
            results = []
            for _, filecat in enumerate(catdict):
                results.append(createfile_bench_cat(filecat))
        #------------------------------------------------

        allcatspc = [spc for result in results for spc in result]
        # Give warning if emissions species is not assigned a benchmark
        # category
        for spc in emis_spc:
            if spc not in allcatspc:
                print(\
                    f"Warning: species {spc} has emissions diagnostics but is not"
                      " in benchmark_categories.yml"
                )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    gc.collect()


def make_benchmark_emis_tables(
        reflist,
        refstr,
        devlist,
        devstr,
        spcdb_files,
        dst="./benchmark",
        benchmark_type="FullChemBenchmark",
        refmet=None,
        devmet=None,
        overwrite=False,
        ref_interval=[2678400.0],
        dev_interval=[2678400.0],
):
    """
    Creates a text file containing emission totals by species and
    category for benchmarking purposes.

    Args:
        reflist: str | list of str
             Path name(s) of the emissions file(s) that constitute
             the "Ref" (aka "Reference") data set.
        refstr: str
            A string to describe ref (e.g. version number)
        devlist: str | list of str
             Path name(s) of the emissions file(s) that constitute
             the "Dev" (aka "Development") data set.
        devstr: str
            A string to describe dev (e.g. version number)
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: ./benchmark
        benchmark_type: str
            A string denoting the type of benchmark output to plot, options are
            FullChemBenchmark, TransportTracersBenchmark or CH4Benchmark.
            Default value: "FullChemBenchmark"
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

    """

    # ==================================================================
    # Initialization
    # ==================================================================
    verify_variable_type(reflist, (str,list))
    verify_variable_type(refstr, str)
    verify_variable_type(devlist, (str,list))
    verify_variable_type(devstr, str)
    verify_variable_type(spcdb_files, list)

    # Create the destination folder
    make_directory(dst, overwrite)

    # Create the "Tables" category folder if it does not exist
    emisdir = os.path.join(dst, "Tables")
    if not os.path.isdir(emisdir):
        os.mkdir(emisdir)

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ==================================================================
    # Read data from netCDF into Dataset objects
    # ==================================================================

    # Read the input datasets
    # Also read the meteorology datasets if passed. These are optional
    # since the refds and devds have variable AREA already (always true)
    # and unit conversions do not require any meteorology.
    if isinstance(reflist, str):
        reflist = [reflist]
    if isinstance(devlist, str):
        devlist = [devlist]
    refmetds = None
    devmetds = None

    refds = xr.open_mfdataset(
        reflist,
        drop_variables=SKIP_THESE_VARS
    )
    devds = xr.open_mfdataset(
        devlist,
        drop_variables=SKIP_THESE_VARS
    )
    if refmet is not None:
        refmetds = xr.open_mfdataset(
            refmet,
            drop_variables=SKIP_THESE_VARS
        )
    if devmet is not None:
        devmetds = xr.open_mfdataset(
            devmet,
            drop_variables=SKIP_THESE_VARS
        )

    # ==================================================================
    # Create table of emissions
    # ==================================================================

    # Emissions species dictionary
    spc_dict = read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            EMISSION_SPC
        ),
        quiet=True
    )
    species=spc_dict[benchmark_type]
    inv_dict = read_config_file(
        os.path.join(
            os.path.dirname(__file__),
            EMISSION_INV
        ),
        quiet=True
    )
    inventories=inv_dict[benchmark_type]

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
        spcdb_files,
        file_emis_totals,
        ref_interval,
        dev_interval,
        template="Emis{}_",
        refmetdata=refmetds,
        devmetdata=devmetds,
    )

    # Create table of emissions by inventory
    create_total_emissions_table(
        refds,
        refstr,
        devds,
        devstr,
        inventories,
        spcdb_files,
        file_inv_totals,
        ref_interval,
        dev_interval,
        template="Inv{}_",
        refmetdata=refmetds,
        devmetdata=devmetds,
    )

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    del refmetds
    del devmetds
    gc.collect()
