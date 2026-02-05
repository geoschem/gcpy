#!/usr/bin/env python3
"""
Specific utilities for creating plots from GEOS-Chem benchmark simulations.
"""
import os
import itertools
import gc
import numpy as np
import pandas as pd
import xarray as xr
from tabulate import tabulate
from gcpy.util import \
    add_missing_variables, \
    array_equals, compare_varnames, dataset_reader, \
    get_filepath, make_directory, \
    read_species_metadata, replace_whitespace, unique_values, \
    verify_variable_type, wrap_text
from gcpy.units import convert_units
from gcpy.constants import \
    ENCODING, SKIP_THESE_VARS, TABLE_WIDTH

# Suppress numpy divide by zero warnings to prevent output spam
np.seterr(divide="ignore", invalid="ignore")


def make_benchmark_operations_budget(
        refstr,
        reffiles,
        devstr,
        devfiles,
        ref_interval,
        dev_interval,
        spcdb_files,
        benchmark_type=None,
        label=None,
        col_sections=["Full", "Trop", "PBL", "FixedLevs", "Strat"],
        operations=["Chemistry", "Convection", "EmisDryDep",
                    "Mixing", "Transport", "WetDep"],
        compute_accum=True,
        compute_restart=False,
        require_overlap=False,
        dst='.',
        species=None,
        overwrite=True,
        verbose=False,
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
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

    Keyword Args (optional):
        benchmark_type: str
            A string denoting the type of benchmark output to plot, options are
            FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
            Default value: None
        label: str
            Contains the date or date range for each dataframe title.
            Default value: None
        col_sections: list of str
            List of column sections to calculate global budgets for. May
            include Strat eventhough not calculated in GEOS-Chem, but Full
            and Trop must also be present to calculate Strat.
            Default value: ["Full", "Trop", "PBL", "FixedLevs", "Strat"]
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
        compute_accum: bool
            Optionally turn on/off accumulation calculation. If True, will
            only compute accumulation if all six GEOS-Chem operations budgets
            are computed. Otherwise a message will be printed warning that
            accumulation will not be calculated.
            Default value: True
        compute_restart: bool
            Optionally turn on/off calculation of mass change based on restart
            file. Only functional for "Full" column section.
            Default value: False
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
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False
    """
    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # ------------------------------------------
    # Column sections
    # ------------------------------------------

    # Print info. Only allow Strat if Trop and Full are present
    print("Column sections:")
    for col_section in col_sections:
        print(f"  {col_section}")
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
    if compute_restart:
        all_operations = gc_operations + ["RESTART"]
    n_ops = len(all_operations)

    # Print info
    print("Operations:")
    for all_operation in all_operations:
        print(f"  {all_operation}")
    if compute_accum:
        if "ACCUMULATION" in all_operations:
            print("*** Will compute ACCUMULATION operation as sum of all "
                  "GEOS-Chem operations ***")
        else:
            print("***Will not compute ACCUMULATION since not all GEOS-Chem"
                  " operation budgets will be computed.")
    if compute_restart:
        print("*** Will compute RESTART operation as mass change "
              "based on simulation start and end restart files ***")

    # ------------------------------------------
    # Read data
    # ------------------------------------------

    # Assume this will be annual budget if interval greater than 3e7 sec
    annual = ref_interval > 3.0e7

    # Read data from disk (either one month or 12 months)
    print('Opening ref and dev data')
    skip_vars = SKIP_THESE_VARS
    if annual:
        ref_ds = xr.open_mfdataset(reffiles, drop_variables=skip_vars)
        dev_ds = xr.open_mfdataset(devfiles, drop_variables=skip_vars)
    else:
        ref_ds = xr.open_dataset(reffiles, drop_variables=skip_vars)
        dev_ds = xr.open_dataset(devfiles, drop_variables=skip_vars)

    # TODO: Add section for reading files for computing mass from restart file

    # ------------------------------------------
    # Species
    # ------------------------------------------

    # Get information about variables in data files
    vardict = compare_varnames(ref_ds, dev_ds, quiet=True)
    refonly = vardict["refonly"]
    devonly = vardict["devonly"]
    cmnvars = vardict["commonvars2D"]

    # Reduce each list to only variables containing "Budget" and not "Strat"
    refonly = [v for v in refonly if "Budget" in v and "Strat" not in v]
    devonly = [v for v in devonly if "Budget" in v and "Strat" not in v]
    cmnvars = [v for v in cmnvars if "Budget" in v and "Strat" not in v]

    # Special handling for fixed level budget diagnostic
    # Get variable name prefix, e.g. Levs1to35. Check that all fixed level
    # vars have the same prefix. Update section names used in table.
    fixedlevvars = [v for v in cmnvars if "Budget" in v and "Levs" in v]
    if fixedlevvars is not None:
        fixedlevnames = [v[v.index('Levs'):].split("_")[0] for v in fixedlevvars]
        if len(set(fixedlevnames)) > 1:
            msg = "Budget fixed level diagnostic name must be constant!"
            raise ValueError(msg)
        col_sections = [v.replace('FixedLevs',fixedlevnames[0]) for v in col_sections]
        gc_sections = [v.replace('FixedLevs',fixedlevnames[0]) for v in gc_sections]

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
    ref_metadata, dev_metadata = read_species_metadata(
        spcdb_files,
        quiet=True
    )

    # Determine what the converted units and conversion factor should be
    # based on benchmark type and species (tracer) name. Assume raw data [kg/s]
    ref_conv_fac = {}
    dev_conv_fac = {}
    units = {}
    is_wetdep = {}
    for spc in spclist:

        # Identify wetdep species
        is_wetdep[spc] = None
        ref_species_metadata = ref_metadata.get(spc)
        dev_species_metadata = dev_metadata.get(spc)
        if ref_species_metadata is not None and \
           dev_species_metadata is not None:
            is_wetdep[spc] = \
                ref_species_metadata.get("Is_WetDep") or \
                dev_species_metadata.get("Is_WetDep")

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

        # Keep track of progress in log
        print(f"  {gc_section}")

        # Loop over species in that section
        for i, spc in enumerate(spclist):

            # Keep track of progress (debugging print)
            #if (i + 1) % 50 == 0:
            #    print(f"  {gc_section}: species {i + 1} of {n_spc}")

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
        print('Computing Strat budgets from Trop and Full')

        # Loop over species
        for i, spc in enumerate(spclist):

            # Keep track of progress (debugging print)
            #if (i + 1) % 50 == 0:
            #    print(f"  Strat: species {i + 1} of {n_spc}")

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

            # Keep track of progress in log
            print(f"  {col_section}")

            # Loop over species
            for i, spc in enumerate(spclist):

                # Keep track of progress (debugging print)
                #if (i + 1) % 50 == 0:
                #    print(f"  {col_section}: species {i + 1} of {n_spc}")

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

    # ------------------------------------------
    # Compute mass change in restarts for each column section (if applicable)
    # ------------------------------------------
    if compute_restart:
        print('Computing RESTART operation budgets...')

        # Read the species database files in the Ref & Dev rundirs,
        # and return a dict containing metadata for each.
        ref_metadata, dev_metadata = read_species_metadata(
            spcdb_files,
            quiet=True
        )

        # Loop over all column sections
        for col_section in col_sections:

            # Loop over species
            for i, spc in enumerate(spclist):

                # Keep track of progress
                if (i + 1) % 50 == 0:
                    print(f"  {col_section}: species {i + 1} of {n_spc}")

                # Get the accumulation dataframe row to fill. Skip if not found
                # of if not Full column section.
                dfrow = (df["Column_Section"] == "Full") \
                    & (df["Species"] == spc) \
                    & (df["Operation"] == "RESTART")
                if not any(dfrow):
                    continue

                # Get ref and dev mass

                # Get species metadata for unit conversion. If none, skip.
                ref_species_metadata = ref_metadata.get(spc)
                dev_species_metadata = dev_metadata.get(spc)
                if ref_species_metadata is None and \
                   dev_species_metadata is None:
                    continue

                # Get molecular weights
                #ref_mol_wt_g = get_molwt_from_metadata(
                #    ref_species_metadata,
                #    spc
                #)
                #dev_mol_wt_g = get_molwt_from_metadata(
                #    ref_species_metadata,
                #    spc
                #)

                # Specify target units
                target_units = "Gg"

                # ==============================================================
                # Convert units of Ref and save to a DataArray
                # (or skip if Ref contains NaNs everywhere)
                # ==============================================================
                refarray = refdata[v]
                if not np.isnan(refdata[v].values).all():
                    refarray = convert_units(
                        refarray,
                        spc,
                        ref_species_metadata,
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
                        dev_species_metadata,
                        target_units,
                        area_m2=met_and_masks["Dev_Area"],
                        delta_p=met_and_masks["Dev_Delta_P"],
                        box_height=met_and_masks["Dev_BxHeight"],
                    )


                # Compute ref mass as end mass minus start mass in ref
                # TODO

                # Compute dev mass as end mass minus start mass in dev
                # TODO - copy above once compete. The rest should just work.

                # Calculate diff and % diff
                if not np.isnan(refmass) and not np.isnan(devmass):
                    diff = devmass - refmass
                    try:
                        pctdiff = diff / refmass * 100
                    except BaseException:
                        pctdiff = np.nan
                else:
                    diff = np.nan
                    pctdiff = np.nan

                # Fill dataframe
                df.loc[dfrow, "Units_converted"] = units[spc]
                df.loc[dfrow, "Ref"] = refmass
                df.loc[dfrow, "Dev"] = devmass
                df.loc[dfrow, "Diff"] = diff
                df.loc[dfrow, "Pct_diff"] = pctdiff

    #  Sanity check write to csv (for testing. Keep commented out otherwise)
    #df.to_csv('df.csv', na_rep='NA')

    # ------------------------------------------
    # Make budget file
    # ------------------------------------------

    # Create the target output directory hierarchy if it doesn't already exist
    make_directory(dst, overwrite)

    # Print budgets to file
    filename = f"{dst}/Budgets_After_Operations.txt"
    with open(filename, "w+", encoding=ENCODING) as f:
        print("#" * 78, file=f)
        if label is not None and benchmark_type is not None:
            print(f"{benchmark_type} budgets for {label}", file=f)
        else:
            print(f"Budgets across {ref_interval}/{dev_interval} sec", file=f)
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
            print(f"{spc} budgets (Ref={refstr}; Dev={devstr})", file=f)

            # Print a table for each column section
            for col_section in col_sections:

                # Get the dataframe rows. Skip if none found.
                dfrows = (df["Column_Section"] == col_section) \
                    & (df["Species"] == spc) \
                    & (df["Operation"].isin(all_operations))
                if not any(dfrows):
                    continue

                # Print dataframe subset to file
                print(f"{col_section} {units[spc]} : {spc}", file=f)
                print(tabulate(
                    df.loc[dfrows,
                           ["Operation",
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
    del df
    del ref_ds
    del dev_ds
    gc.collect()


def create_benchmark_summary_table(
        refpath,
        refstr,
        refdate,
        devpath,
        devstr,
        devdate,
        collections,
        dst="./benchmark",
        overwrite=False,
        outfilename="Summary.txt",
        verbose=False,
        ref_gchp=False,
        dev_gchp=False
):
    """
    Creates a benchmark summary table that shows which data collections
    have difference.  Useful for scanning the 1-hr and 1-month benchmark
    outputs.

    Args:
        refpath: str
            Path to the first data set to be compared (aka "Ref").
        refstr: str
            A string that can be used to identify refdata
            (e.g. a model version number or other identifier).
        refdate: np.datetime64
            Date/time stamp used by the "Ref" data files.
        ref_gchp: bool
            Set to True if the "Ref" data comes from a GCHP run.
            Default value: False
        devpath: str
            Path to the second data set to be compared (aka "Dev").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        dev_gchp: bool
            Set to True if the "Ref" data comes from a GCHP run.
            Default value: False
        devdate: np.datetime64
            Date/time stamp used by the "Dev" data files.
        collections: list of strings
            List of diagnostic collections to examine.

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: "./benchmark"
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "Summary.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.

        Species properties (such as molecular weights) are read from a
        YAML file called "species_database.yml".
    """

    # ==================================================================
    # Open file for output
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Create the directory for output
    make_directory(dst, overwrite)

    # Create file
    try:
        f = open(os.path.join(dst, outfilename), "w", encoding=ENCODING)
    except (IOError, OSError, FileNotFoundError) as e:
        msg = f"Could not open {outfilename} for writing!"
        raise e(msg) from e

    # Title strings
    title1 = "### Benchmark summary table"
    title2 = f"### Ref = {refstr}"
    title3 = f"### Dev = {devstr}"

    # Print header to file
    print("#" * 80, file=f)
    print(f"{title1 : <77}{'###'}", file=f)
    print(f"{'###'  : <77}{'###'}", file=f)
    print(f"{title2 : <77}{'###'}", file=f)
    print(f"{title3 : <77}{'###'}", file=f)
    print("#" * 80, file=f)
    print(file=f)

    # ==================================================================
    # Read data and look differences btw Ref & Dev versions
    # ==================================================================

    # Variables to skip
    skip_vars = SKIP_THESE_VARS
    skip_vars.append("AREA")

    # Pick the proper function to read the data
    reader = dataset_reader(
        multi_files=False,
        verbose=verbose
    )

    # Loop over diagnostic files
    for col in collections:

        # Read Ref data
        refdata = reader(
            get_filepath(
                refpath,
                col,
                refdate,
                is_gchp=ref_gchp
            ),
            drop_variables=skip_vars
        )

        # Get Dev data
        devdata = reader(
            get_filepath(
                devpath,
                col,
                devdate,
                is_gchp=dev_gchp
            ),
            drop_variables=skip_vars
        )

        # Make sure that Ref and Dev datasets have the same variables.
        # Variables that are in Ref but not in Dev will be added to Dev
        # with all missing values (NaNs). And vice-versa.
        [refdata, devdata] = add_missing_variables(
            refdata,
            devdata
        )

        # Find all common variables between the two datasets
        vardict = compare_varnames(
            refdata,
            devdata,
            quiet=True
        )

        # List of differences for this collection
        diff_list = []

        # Keep track of which variables are different
        # NOTE: Use 32-point float for comparisons since this is
        # the precision used for History diagnostics.
        for v in vardict["commonvarsData"]:
            if not array_equals(
                    refdata[v],
                    devdata[v],
                    dtype=np.float32
            ):
                diff_list.append(v)

        # Drop duplicate values from diff_list
        diff_list = unique_values(diff_list, drop=[None])

        if len(diff_list) == 0:
            print("-" *  79, file=f)
            print(f"{col}: {devstr} is identical to {refstr}", file=f)
            print(file=f)
        else:
            print("-" *  79, file=f)
            print(f"{col}: {devstr} differs from {refstr}", file=f)
            print("\n  Diagnostics that differ", file=f)
            for i, v in enumerate(diff_list):
                print(f"    {v}", file=f)
                if i > 10:
                    print(f"    ... and {len(diff_list) - 10} others", file=f)
                    break
            print(file=f)

    # ==================================================================
    # Close files
    # ==================================================================
    f.close()


def diff_list_to_text(
        refstr,
        devstr,
        diff_list,
        fancy_format=False
):
    """
    Converts a list of species/emissions/inventories/diagnostics that
    show differences between GEOS-Chem versions ot a printable text
    string.

    Args:
    -----
    diff_list : list
        List to be converted into text.  "None" values will be dropped.
    fancy_format: bool
        Set to True if you wish output text to be bookended with '###'.

    Returns:
    diff_text : str
        String with concatenated list values.
    """
    verify_variable_type(diff_list, list)

    # Use "Dev" and "Ref" for inserting into a header
    if fancy_format:
        refstr = "Ref"
        devstr = "Dev"

    # Strip out duplicates from diff_list
    # Prepare a message about species differences (or alternate msg)
    diff_list = unique_values(diff_list, drop=[None])

    # Print the text
    n_diff = len(diff_list)
    if n_diff > 0:
        diff_text = f"{devstr} and {refstr} show {n_diff} differences"
    else:
        diff_text = f"{devstr} and {refstr} are identical"

    # If we are placing the text in a header, trim the length of diff_text
    # to fit.  NOTE: TABLE_WIDTH-7 leaves room for the '### ' at the start
    # of the string and the '###' at the end of the string,
    if fancy_format:
        diff_text = f"### {diff_text : <{TABLE_WIDTH-7}}{'###'}"
        diff_text = wrap_text(
            diff_text,
            width=TABLE_WIDTH
        )

    return diff_text.strip()


def diff_of_diffs_toprow_title(config, model):
    """
    Creates the diff-of-diffs plot title for the top row of the
    six-plot output.  If the title string is too long (as empirically
    determined), then a newline will be inserted in order to prevent
    the title strings from overlapping.

    Args:
    -----
    config : dict
       Dictionary containing the benchmark options (as read from a
       YAML file such as 1mo_benchmark.yml, etc.)
    model: str
       The model to plot.  Accepted values are "gcc" or "gchp".

    Returns:
    --------
    title: str
        The plot title string for the diff-of-diff
    """
    verify_variable_type(config, dict)
    verify_variable_type(model, str)
    if not "gcc" in model and not "gchp" in model:
        msg = "The 'model' argument must be either 'gcc' or 'gchp'!"
        raise ValueError(msg)

    title = (
        config["data"]["dev"][model]["version"]
        + " - "
        + config["data"]["ref"][model]["version"]
    )

    if len(title) > 40:
        title = (
            config["data"]["dev"][model]["version"]
            + " -\n"
            + config["data"]["ref"][model]["version"]
        )

    return title


def create_benchmark_sanity_check_table(
        devpath,
        devstr,
        devdate,
        collections,
        dst="./benchmark",
        is_gchp=False,
        overwrite=False,
        outfilename="Diagnostic_Sanity_Check.txt",
        verbose=False,
):
    """
    Creates a diagnostic sanity check table that shows which diagnostic
    variables are zero or NaN everywhere.  This can help to identify
    bugs in diagnostic output.

    Args:
        devpath: str
            Path to the data set to be compared (aka "Dev").
        devstr: str
            A string that can be used to identify the data set specified
            by devfile (e.g. a model version number or other identifier).
        devdate: np.datetime64
            Date/time stamp used by the "Dev" data files.
        collections: list of strings
            List of diagnostic collections to examine.

    Keyword Args (optional):
        dst: str
            A string denoting the destination folder where the file
            containing emissions totals will be written.
            Default value: "./benchmark"
        is_gchp : bool
           Set this flag to true to denote if the data is from GCHP.
        overwrite: bool
            Set this flag to True to overwrite files in the
            destination folder (specified by the dst argument).
            Default value: False
        outfilename: str
            Name of the text file which will contain the table of
            emissions totals.
            Default value: "Summary.txt"
        verbose: bool
            Set this switch to True if you wish to print out extra
            informational messages.
            Default value: False

    Remarks:
        This method is mainly intended for model benchmarking purposes,
        rather than as a general-purpose tool.
    """

    # ==================================================================
    # Initial preparations
    # ==================================================================

    # Replace whitespace in the ref and dev labels
    devstr = replace_whitespace(devstr)

    # Create the directory for output (if necessary)
    make_directory(dst, overwrite)
    outfilename = os.path.join(dst, outfilename)

    # Pick the proper function to read the data
    reader = dataset_reader(
        multi_files=False,
        verbose=verbose
    )

    # Variables to skip
    skip_vars = SKIP_THESE_VARS
    skip_vars.append("AREA")

    # ==================================================================
    # Open output file and write header
    # ==================================================================
    with open(outfilename, "w", encoding=ENCODING) as ofile:

        # Title strings
        title1 = "### Benchmark diagnostic sanity check table"
        title2 = f"### Dev = {devstr}"

        # Print header to file
        print("#" * 80, file=ofile)
        print(f"{title1 : <77}{'###'}", file=ofile)
        print(f"{'###'  : <77}{'###'}", file=ofile)
        print(f"{title2 : <77}{'###'}", file=ofile)
        print("#" * 80, file=ofile)

        # ==============================================================
        # Loop over diagnostic collections and scan files
        # ==============================================================
        for col in collections:

            # Read data into an xr.DataSet object
            file_name = get_filepath(
                devpath,
                col,
                devdate,
                is_gchp=is_gchp,
            )
            dset = reader(
                file_name,
                drop_variables=skip_vars
            )

            # Determine which variables are all zeroes or NaN
            all_zeros_or_nans = []
            for var in dset.data_vars:
                data = dset[var].values
                if np.all(data == 0) or np.all(data == np.nan):
                    all_zeros_or_nans.append(var)

            # ===========================================================
            # Print results for each collection
            # ===========================================================
            print("", file=ofile)
            print("="*80, file=ofile)
            print(f"{os.path.basename(file_name)}", file=ofile)
            print("="*80, file=ofile)
            print("", file=ofile)

            if len(all_zeros_or_nans) == 0:
                print("No variables were all zero or all NaN", file=ofile)
            else:
                print("These variables were all zero or all NaN:", file=ofile)
                for var in all_zeros_or_nans:
                    print(f"   {var}", file=ofile)
