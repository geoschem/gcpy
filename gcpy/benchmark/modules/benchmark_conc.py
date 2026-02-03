#!/usr/bin/env python3
"""
Creates concentration plots from GEOS-Chem benchmark simulation output.
"""
import os
import warnings
import gc
import xarray as xr
import sparselt.esmf
import sparselt.xr
from joblib import \
    Parallel, delayed
from gcpy.regrid import create_regridders
from gcpy.util import \
    add_bookmarks_to_pdf, add_missing_variables, \
    add_nested_bookmarks_to_pdf, dataset_reader, dataset_mean, \
    make_directory, replace_whitespace, verify_variable_type
from gcpy.constants import \
    ENCODING, SKIP_THESE_VARS
from gcpy.plot.compare_single_level import compare_single_level
from gcpy.plot.compare_zonal_mean import compare_zonal_mean
from gcpy.benchmark.modules.benchmark_utils import \
    add_lumped_species_to_dataset, \
    archive_lumped_species_definitions, get_species_categories, \
    archive_species_categories, rename_speciesconc_to_speciesconcvv


def make_benchmark_conc_plots(
        ref,
        refstr,
        dev,
        devstr,
        spcdb_files,
        dst="./benchmark",
        subdst=None,
        overwrite=False,
        verbose=False,
        collection="SpeciesConc",
        benchmark_type="FullChemBenchmark",
        cmpres=None,
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
        spcdb_files : list
            Paths to species_database.yml files in Ref & Dev rundirs

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
            A string denoting the type of benchmark output to plot, options are
            FullChemBenchmark, TransportTracersBenchmark, or CH4Benchmark.
            Default value: "FullChemBenchmark"
        cmpres: string
            Grid resolution at which to compare ref and dev data, e.g. '1x1.25'
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
        time_mean : bool
            Determines if we should average the datasets over time
            Default value: False
    """
    # ==================================================================
    # Initialization and data read
    # ==================================================================
    verify_variable_type(ref, (str, list))
    verify_variable_type(refstr, str)
    verify_variable_type(dev, (str, list))
    verify_variable_type(devstr, str)
    verify_variable_type(spcdb_files, list)

    # Create the destination folder
    make_directory(dst, overwrite)

    # Define extra title text (usually a date string)
    # for the top-title of the plot
    if subdst is not None:
        extra_title_txt = subdst
    else:
        extra_title_txt = None

    # Replace whitespace in the ref and dev labels
    refstr = replace_whitespace(refstr)
    devstr = replace_whitespace(devstr)

    # Pick the proper function to read the data
    reader = dataset_reader(time_mean, verbose=verbose)

    # Open datasets
    refds = reader(ref, drop_variables=SKIP_THESE_VARS)
    devds = reader(dev, drop_variables=SKIP_THESE_VARS)

    # Rename SpeciesConc_ to SpeciesConcVV_ for consistency with new
    # naming introduced in GEOS-Chem 14.1.0
    refds = rename_speciesconc_to_speciesconcvv(refds)
    devds = rename_speciesconc_to_speciesconcvv(devds)

    # -----------------------------------------------------------------
    # Kludge, rename wrong variable name
    if "SpeciesConcVV_PFE" in refds.data_vars.keys():
        refds = refds.rename({"SpeciesConcVV_PFE": "SpeciesConcVV_pFe"})
    if "SpeciesConcVV_PFE" in devds.data_vars.keys():
        devds = devds.rename({"SpeciesConcVV_PFE": "SpeciesConcVV_pFe"})
    # -----------------------------------------------------------------

    if verbose:
        print('\nPrinting refds (comparison ref)\n')
        print(refds)
        print('\nPrinting devds (comparison dev)\n')
        print(devds)

    # Open met datasets if passed as arguments
    refmetds = None
    devmetds = None
    if refmet:
        refmetds = reader(refmet, drop_variables=SKIP_THESE_VARS)
    if devmet:
        devmetds = reader(devmet, drop_variables=SKIP_THESE_VARS)

    # Determine if doing diff-of-diffs
    diff_of_diffs = second_ref is not None and second_dev is not None

    if diff_of_diffs:

        # --------------------------------------------------------------
        # %%%%% We are plotting diff-of-diffs %%%%%
        #
        # Open the second Ref and Dev datasets, if they have been
        # passed as keyword arguments, as these are needed for the
        # diff-of-diffs plots.  Regrid to the same # horizontal grid
        # resolution if the two Refs or two Devs are not on the same
        # grid.
        #
        # Also, do not take the time mean (e.g. Annual Mean) of the
        # datasets, as this will compute the the difference of means.
        # Instead, we need to compute the mean of differences.  This
        # will be done in the plotting functions compare_single_level
        # and compare_zonal_mean.
        # --------------------------------------------------------------
        second_refds = reader(second_ref, drop_variables=SKIP_THESE_VARS)
        second_devds = reader(second_dev, drop_variables=SKIP_THESE_VARS)

        if verbose:
            print('\nPrinting second_refds (dev of ref for diff-of-diffs)\n')
            print(second_refds)
            print('\nPrinting second_devds (dev of dev for diff-of-diffs)\n')
            print(second_devds)

        # Only regrid if ref grid resolutions differ.
        # Assume only may differ if both cubed-sphere.
        # If different, assume the two resolutions are C24 and C48,
        # and do the comparison at C48.
        regrid_ref = False
        if "Xdim" in refds.dims and "Xdim" in second_refds.dims:
            if refds['Xdim'].size != second_refds['Xdim'].size:
                regrid_ref = True
        regrid_dev = False
        if "Xdim" in devds.dims and "Xdim" in second_devds.dims:
            if devds['Xdim'].size != second_devds['Xdim'].size:
                regrid_dev = True

        if regrid_ref or regrid_dev:
            # Assume regridding C24 to C48 to compute the difference at C48
            regridfile=os.path.join(weightsdir,'regrid_weights_c24_to_c48.nc')
            transform = sparselt.esmf.load_weights(
                regridfile,
                input_dims=[('nf', 'Ydim', 'Xdim'), (6, 24, 24)],
                output_dims=[('nf', 'Ydim', 'Xdim'), (6, 48, 48)],
            )
            if regrid_ref and refds['Xdim'].size == 24:
                print('\nRegridding ref dataset 1 from C24 to C48\n')
                refds = sparselt.xr.apply(transform, refds)
                print(refds)
                print('\nRegrid complete\n')
            if regrid_ref and second_refds['Xdim'].size == 24:
                print('\nRegridding ref dataset 2 from C24 to C48\n')
                second_refds = sparselt.xr.apply(transform, second_refds)
                print(second_refds)
                print('\nRegrid complete\n')
            if regrid_dev and devds['Xdim'].size == 24:
                print('\nRegridding dev dataset 1 from C24 to C48\n')
                devds = sparselt.xr.apply(transform, devds)
                print(devds)
                print('\nRegrid complete\n')
            if regrid_dev and second_devds['Xdim'].size == 24:
                print('\nRegridding dev dataset 2 from C24 to C48\n')
                second_devds = sparselt.xr.apply(transform, second_devds)
                print(second_devds)
                print('\nRegrid complete\n')

    else:

        # --------------------------------------------------------------
        # %%%%% We are not plotting diff-of-diffs %%%%%
        #
        # We do not need the second Ref and Dev data.
        # We can also compute the time mean (e.g. Annual Mean) here.
        # --------------------------------------------------------------
        second_refds = None
        second_devds = None
        if time_mean:
            refds = dataset_mean(refds)
            devds = dataset_mean(devds)
            refmetds = dataset_mean(refmetds)
            devmetds = dataset_mean(devmetds)

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
        [refds, devds] = add_missing_variables(refds, devds)
        var_prefix = 'SpeciesConcVV_'
        varlist = [k for k in refds.data_vars.keys() if var_prefix in k]
        varlist.sort()

        # Surface
        pdfname = os.path.join(dst, 'SpeciesConc_Sfc.pdf')
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            varlist=varlist,
            cmpres=cmpres,
            pdfname=pdfname,
            use_cmap_RdBu=use_cmap_RdBu,
            log_color_scale=log_color_scale,
            extra_title_txt=extra_title_txt,
            normalize_by_area=normalize_by_area,
            weightsdir=weightsdir,
            second_ref=second_refds,
            second_dev=second_devds,
            spcdb_files=spcdb_files
        )

        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=var_prefix,
            verbose=verbose
        )

        # 500 hPa
        pdfname = os.path.join(dst, 'SpeciesConc_500hPa.pdf')
        compare_single_level(
            refds,
            refstr,
            devds,
            devstr,
            ilev=22,
            varlist=varlist,
            cmpres=cmpres,
            pdfname=pdfname,
            use_cmap_RdBu=use_cmap_RdBu,
            log_color_scale=log_color_scale,
            normalize_by_area=normalize_by_area,
            extra_title_txt=extra_title_txt,
            weightsdir=weightsdir,
            second_ref=second_refds,
            second_dev=second_devds,
            spcdb_files=spcdb_files
        )

        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=var_prefix,
            verbose=verbose
        )

        # Zonal mean
        pdfname = os.path.join(dst, 'SpeciesConc_ZnlMn.pdf')
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
            weightsdir=weightsdir,
            second_ref=second_refds,
            second_dev=second_devds,
            spcdb_files=spcdb_files
        )

        add_bookmarks_to_pdf(
            pdfname,
            varlist,
            remove_prefix=var_prefix,
            verbose=verbose
        )
        return

    # ==================================================================
    # Otherwise plot by categories. Convert units to ug/m3 for
    # aerosol categories: Aerosols and Secondary Organic Aerosols.
    # ==================================================================

    # FullChemBenchmark has lumped species (TransportTracers, CH4 do not)
    if "FullChem" in benchmark_type:
        print("\nComputing lumped species for full chemistry benchmark")

        print("-->Adding lumped species to ref dataset")
        refds = add_lumped_species_to_dataset(refds)

        print("-->Adding lumped species to dev dataset")
        devds = add_lumped_species_to_dataset(devds)

        if diff_of_diffs:
            print("-->Adding lumped species to second ref and dev datasets")
            second_refds = add_lumped_species_to_dataset(second_refds)
            second_devds = add_lumped_species_to_dataset(second_devds)

        archive_lumped_species_definitions(dst)
        print("Lumped species computation complete.\n")

    # Get the list of species categories
    catdict = get_species_categories(benchmark_type)
    archive_species_categories(dst)

    # Make sure that Ref and Dev datasets have the same variables.
    # Variables that are in Ref but not in Dev will be added to Dev
    # with all missing values (NaNs). And vice-versa.
    [refds, devds] = add_missing_variables(refds, devds)

    if diff_of_diffs:
        [refds, second_refds] = add_missing_variables(refds, second_refds)
        [devds, second_devds] = add_missing_variables(devds, second_devds)

    # Collection prefix
    if "SpeciesConc" in collection:
        coll_prefix = "SpeciesConcVV_"
    else:
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

        # Get the list of variables in both Ref and Dev for each category
        # (this is computationally efficient)
        ref_vars = set(refds.data_vars)
        dev_vars = set(devds.data_vars)
        candidates = [
            coll_prefix + spc
            for subcat in catdict[filecat]
            for spc in catdict[filecat][subcat]
        ]
        varlist = \
            [var for var in candidates \
             if var in ref_vars and var in dev_vars
        ]
        warninglist = \
            [var for var in candidates \
             if var not in ref_vars or var not in dev_vars
        ]

        # Create new datasets containing only the variables for a
        # given category, as this will optimize performance.
        refds_cat = refds[varlist]
        devds_cat = devds[varlist]
        second_refds_cat = None
        if second_refds is not None:
            second_refds_cat = second_refds[varlist]
        second_devds_cat = None
        if second_devds is not None:
            second_devds_cat = second_devds[varlist]

        # -----------------------
        # Surface plots
        # -----------------------
        if "sfc" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_Surface_{subdst}.pdf"
                )
                print(f"creating {pdfname}")
            else:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_Surface.pdf"
                )

            diff_sfc = []
            compare_single_level(
                refds_cat,
                refstr,
                devds_cat,
                devstr,
                varlist=varlist,
                ilev=0,
                cmpres=cmpres,
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
                second_ref=second_refds_cat,
                second_dev=second_devds_cat,
                n_job=n_job,
                spcdb_files=spcdb_files,
            )
            diff_sfc[:] = [v.replace(coll_prefix, "") for v in diff_sfc]
            cat_diff_dict['sfc'] = diff_sfc
            add_nested_bookmarks_to_pdf(
                pdfname,
                filecat,
                catdict,
                warninglist,
                remove_prefix=coll_prefix
            )

        # -----------------------
        # 500 hPa plots
        # -----------------------
        if "500hpa" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_500hPa_{subdst}.pdf"
                )
            else:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_500hPa.pdf"
                )

            diff_500 = []
            compare_single_level(
                refds_cat,
                refstr,
                devds_cat,
                devstr,
                varlist=varlist,
                ilev=22,
                cmpres=cmpres,
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
                second_ref=second_refds_cat,
                second_dev=second_devds_cat,
                n_job=n_job,
                spcdb_files=spcdb_files
            )
            diff_500[:] = [v.replace(coll_prefix, "") for v in diff_500]
            #dict_500[filecat] = diff_500
            cat_diff_dict['500'] = diff_500
            add_nested_bookmarks_to_pdf(
                pdfname,
                filecat,
                catdict,
                warninglist,
                remove_prefix=coll_prefix
            )

        # -----------------------
        # Zonal mean plots
        # -----------------------
        if "zonalmean" in plots or "zm" in plots:

            if subdst is not None:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_FullColumn_ZonalMean_{subdst}.pdf"
                )
            else:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_FullColumn_ZonalMean.pdf"
                )

            diff_zm = []
            compare_zonal_mean(
                refds_cat,
                refstr,
                devds_cat,
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
                second_ref=second_refds_cat,
                second_dev=second_devds_cat,
                n_job=n_job,
                spcdb_files=spcdb_files
            )
            diff_zm[:] = [v.replace(coll_prefix, "") for v in diff_zm]
            #dict_zm = diff_zm
            cat_diff_dict['zm'] = diff_zm
            add_nested_bookmarks_to_pdf(
                pdfname,
                filecat,
                catdict,
                warninglist,
                remove_prefix=coll_prefix
            )

            # Strat_ZonalMean plots will use a log-pressure Y-axis, with
            # a range of 1..100 hPa, as per GCSC request. (bmy, 8/13/19)
            if subdst is not None:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_Strat_ZonalMean_{subdst}.pdf"
                )
            else:
                pdfname = os.path.join(
                    catdir,
                    f"{filecat}_Strat_ZonalMean.pdf"
                )

            compare_zonal_mean(
                refds_cat,
                refstr,
                devds_cat,
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
                second_ref=second_refds_cat,
                second_dev=second_devds_cat,
                n_job=n_job,
                spcdb_files=spcdb_files
            )
            add_nested_bookmarks_to_pdf(
                pdfname,
                filecat,
                catdict,
                warninglist,
                remove_prefix=coll_prefix
            )
        return {filecat: cat_diff_dict}

    # --------------------------------------------
    # Create the plots in parallel
    # Turn off parallelization if n_job=1
    if n_job != 1:
        results = Parallel(n_jobs=n_job)(
            delayed(createplots)(filecat)
            for _, filecat in enumerate(catdict)
        )
    else:
        results = []
        for _, filecat in enumerate(catdict):
            results.append(createplots(filecat))
    # --------------------------------------------

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
                    with open(filename, "a+", encoding=ENCODING) as f:
                        for c, diff_list in dict_sfc.items():
                            print(f"* {c}: ", file=f, end="")
                            for v in diff_list:
                                print(f"{v} ", file=f, end="")
                            print(file=f)
                        f.close()

            if "500hpa" in plots:
                if "500hpa" in filename:
                    with open(filename, "a+", encoding=ENCODING) as f:
                        for c, diff_list in dict_500.items():
                            print(f"* {c}: ", file=f, end="")
                            for v in diff_list:
                                print(f"{v} ", file=f, end="")
                            print(file=f)
                        f.close()

            if "zonalmean" in plots or "zm" in plots:
                if "zonalmean" in filename or "zm" in filename:
                    with open(filename, "a+", encoding=ENCODING) as f:
                        for c, diff_list in dict_zm.items():
                            print(f"* {c}: ", file=f, end="")
                            for v in diff_list:
                                print(f"{v} ", file=f, end="")
                            print(file=f)
                        f.close()

    # -------------------------------------------
    # Clean up
    # -------------------------------------------
    del refds
    del devds
    del refmetds
    del devmetds
    del second_ref
    del second_dev
    gc.collect()
