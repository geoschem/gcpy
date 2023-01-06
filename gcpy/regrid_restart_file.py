"""
This module takes a restart file, regridding weights generated by
ESMF_RegridWeightGen, a template file, and optional stretched-grid parameters,
to produce a regridded GCHP restart file.

Example:
    First create source and target grid specifications using `gridspec`
    (https://github.com/liambindle/gridspec), then create regridding weights
    from source to target grid using ESMF_RegridWeightGen from ESMF, then:

        $ python -m gcpy.regrid_restart_file \
                GEOSChem.Restart.fullchem.20190701_0000z.nc4 \
                4x5_to_c24_weights.nc \
                GCHP.Restart.fullchem.20190701_0000z.c24.nc4

    Or, for a stretched-grid:

        $ python -m gcpy.regrid_restart_file \
                --stretched-grid \
                --stretch-factor=2.0 \
                --target-latitude=32.0 \
                --target-longitude=-64.0 \
                GEOSChem.Restart.fullchem.20190701_0000z.nc4 \
                4x5_to_c24_weights.nc \
                GCHP.Restart.fullchem.20190701_0000z.c24.nc4

"""
import argparse
import logging
import os
from pathlib import Path
import re
import tempfile
import xarray as xr
import numpy as np
import sparselt.esmf
import sparselt.xr
import requests


TEMP_FILES = []


def file_path(path):
    """
    Checks whether or not a regular file exists at the passed path.

    Args:
        file_path (str): A path to a file.

    Returns:
        bool: True if a regular file exists at `file_path`.

    """
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError
    return path


def parse_command_line():
    """
    Parses command line arguments and options into a useful data structure.

    Returns:
        argparse.Namespace: A dict-like object containing command line
                            argument and option values.
    """
    description_text = (
        "regrid_restart_file - Regrid GCHP restart files"
        "\n\n"
        "regrid_restart_file is a tool for regridding  GCHP restart "
        "files. You can resize restart files from their original "
        "resolution to a new resolution, stretch an unstretched restart "
        "file, unstretch a stretched restart file, and re-stretch a "
        "stretched restart file."
        "\n\n"
        "To use this tool, you must first generate regridding weights for "
        "the regridding you would like to carry out, using "
        "ESMF_RegridWeightGen."
        "\n\n"
        "NOTE: GC-Classic regridding is not currently supported by this"
        "tool. If this is something you would like to be supported, please "
        "raise an issue via "
        "https://github.com/geoschem/gcpy/issues/new/choose"
    )

    epilog_text = (
        "Example usage (unstretched grid resizing): "
        "\n\n"
        "python -m gcpy.regrid_restart_file \\ "
        "\n\tGEOSChem.Restart.20190701_0000z.c90.nc4 \\ "
        "\n\tC90_to_C48_weights.nc \\ "
        "\n\tGEOSChem.Restart.20190701_0000z.c90.nc4"
        "\n\n"
        "Example usage (stretching a grid): "
        "\n\n"
        "python -m gcpy.regrid_restart_file \\ "
        "\n\t--stretched-grid \\ "
        "\n\t--stretch-factor 2.0 \\ "
        "\n\t--target-latitude 32.0 \\ "
        "\n\t--target-longitude -64.0 \\ "
        "\n\tGEOSChem.Restart.20190701_0000z.c90.nc4 \\ "
        "\n\tC90_to_C48_weights.nc \\ "
        "\n\tGEOSChem.Restart.20190701_0000z.c90.nc4"
    )

    parser = argparse.ArgumentParser(
        description=description_text,
        epilog=epilog_text,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "file_to_regrid",
        type=file_path,
        metavar="file_to_regrid",
        help="The GEOS-Chem restart file to be regridded",
    )
    parser.add_argument(
        "regridding_weights_file",
        type=file_path,
        metavar="regridding_weights_file",
        help=(
            "The regridding weights file for this regridding, generated "
            "by ESMF_RegridWeightGen"
        ),
    )
    parser.add_argument(
        "template_file",
        type=file_path,
        metavar="template_file",
        help=(
            "The GEOS-Chem restart file to use as a template for regridding - "
            "attributes, dimensions, and variables for the output file will "
            "be taken from this template"
        ),
    )

    parser.add_argument(
        "--stretched-grid",
        action="store_true",
        help=(
            "Create a stretched-grid restart file - you must also pass "
            "stretched-grid parameters!"
        ),
    )

    parser.add_argument(
        "--stretch-factor",
        type=np.float32,
        metavar="stretch_factor",
        help="The stretch factor, if creating a stretched-grid restart file",
        required=False,
    )
    parser.add_argument(
        "--target-latitude",
        type=np.float32,
        metavar="target_latitude",
        help="The target latitude, if creating a stretched-grid restart file",
        required=False,
    )
    parser.add_argument(
        "--target-longitude",
        type=np.float32,
        metavar="target_longitude",
        help="The target longitude, if creating a stretched-grid restart file",
        required=False,
    )

    return parser.parse_args()


def cleanup_tempfile():
    """
    Clean up temporary files created as part of the regridding.
    """
    global TEMP_FILES
    if len(TEMP_FILES) > 0:
        logging.debug("Deleting %d temp files", len(TEMP_FILES))
    for filepath in TEMP_FILES:
        Path(filepath).unlink(missing_ok=True)
    TEMP_FILES = []


def is_gchp_restart_file(dataset):
    """
    Checks whether or not an xarray dataset represents a GCHP restart file.

    Returns:
        bool: True if `dataset` represents a GCHP restart file.
    """
    is_gchp_restart = "SPC_O3" in dataset.data_vars
    is_gcclassic = "SpeciesRst_O3" in dataset.data_vars
    if not any((is_gchp_restart, is_gcclassic)):
        raise ValueError(
            "Couldn't determine if the provided file is a GC-Classic or GCHP "
            "restart file."
        )
    return is_gchp_restart


def open_dataset(file_or_url, chunk_size=8192):
    """
    Open a NetCDF-4 dataset from either file path or URL.

    Args:
        file_or_url (str): A file path on the local system or URL
        chunk_size  (int): Size of chunks to stream from remote dataset to
                           the local system.

    Returns:
        xarray.Dataset: An xarray dataset.
    """
    global TEMP_FILES
    is_url = bool(re.match(r"https?://", file_or_url))
    if is_url:
        logging.debug("Downloading %s", file_or_url)
        with requests.get(file_or_url, stream=True, timeout=30.0) as request:
            request.raise_for_status()  # raise HTTPError
            tempfile_fd, tempfile_path = tempfile.mkstemp()
            with open(tempfile_fd, "wb") as outfile:
                bytes_downloaded = 0
                for chunk in request.iter_content(chunk_size=chunk_size):
                    bytes_downloaded += len(chunk)
                    outfile.write(chunk)
            TEMP_FILES.append(tempfile_path)
    dataset_file_path = tempfile_path if is_url else file_or_url
    logging.debug("Opening %s", dataset_file_path)
    return xr.open_dataset(dataset_file_path)


def rename_variables(dataset, to_gchp=True):
    """
    Rename variables in passed dataset to match either GC-Classic or GCHP
    naming conventions.

    Args:
        datase  (xarray.Dataset): The dataset to have its variables renamed.
        to_gchp (bool)          : True if converting to GCHP naming convention,
                                  False if converting to GC-Classic.
    """
    to_gchp_re_sub = [
        (r"SpeciesRst_(.+)", r"SPC_\1"),
        (r"Met_(.+)", r"\1"),
        (r"Met_DELPDRY", r"DELP_DRY"),
        (
            r"Chem_(WetDepNitrogen|DryDepNitrogen|H2O2AfterChem|SO2AfterChem|KPPHvalue)",
            r"\1",
        ),
    ]
    to_gcclassic_re_sub = [
        (r"SPC_(.+)", r"SpeciesRst_\1"),
        (r"(TropLev|BXHEIGHT)", r"Met_\1"),
    ]
    re_sub_arg_list = to_gchp_re_sub if to_gchp else to_gcclassic_re_sub

    rename_dict = {}
    for re_sub_args in re_sub_arg_list:
        rename_dict.update(
            {
                name: re.sub(*re_sub_args, name)
                for name in dataset.data_vars
                if re.match(re_sub_args[0], name)
            }
        )
    logging.info("Renaming %d variables", len(rename_dict))
    return dataset.rename(rename_dict)


def reverse_lev(dataset):
    """
    Reverse the level index of the passed dataset.

    Args:
        dataset (xarray.Dataset): The dataset to have its level index reversed.

    Returns:
        xarray.Dataset: The input dataset with a reversed level index.
    """
    logging.info("Reversing coordinate 'lev'")
    dataset = dataset.reindex(lev=dataset.lev[::-1])
    dataset = dataset.assign_coords(lev=dataset.lev.values[::-1])
    return dataset


def drop_variables(dataset, output_template):
    """
    Drop variables in the passed dataset which aren't present in the regridding
    output template.

    Args:
        dataset         (xarray.Dataset): The dataset from which to drop
                                          variables.
        output_template (xarray.Dataset): The template from which to determine
                                          variables to drop.

    Returns:
        xarray.Dataset: The input dataset with variables dropped.
    """
    input_var_set = set(dataset.data_vars)
    output_var_set = set(output_template.data_vars)
    drop_vars = input_var_set - output_var_set
    missing_vars = output_var_set - input_var_set
    if len(drop_vars) > 0:
        info_message = (
            "Dropping %d variables from the input restart file "
            "that don't exist in the output template"
        )
        logging.info(info_message, len(drop_vars))

        debug_message = (
            "Variables being dropped from the input restart file: %s"
        )
        logging.debug(debug_message, drop_vars)

        dataset = dataset.drop(drop_vars)
    if len(missing_vars) > 0:
        warning_message = (
            "The input restart file is missing %d variables "
            "that exist in the output template"
        )
        logging.warning(warning_message, len(missing_vars))

        debug_message = "Variables missing in the input restart file: %s"
        logging.debug(debug_message, missing_vars)

        output_template = output_template.drop(missing_vars)
    return dataset, output_template


def regrid(dataset, output_template, weights_file):
    """
    Calculate and apply the regridding, based on passed regridding weights
    and input dataset attributes.

    Args:
        dataset         (xarray.Dataset): The dataset to be regridded.
        output_template (xarray.Dataset): The template file for the regridded
                                          output.
        weights_file    (xarray.Dataset): The precalculated regridding weights,
                                          generated by ESMF_RegridWeightGen.

    Returns:
        xarray.Dataset: The regridded dataset.
    """
    weights = open_dataset(weights_file)
    input_dims = [("lat", "lon"), (dataset.dims["lat"], dataset.dims["lon"])]

    output_template_shape = (
        output_template.dims["lat"],
        output_template.dims["lon"],
    )
    resize_output_template = (
        np.prod(output_template_shape) != weights.dst_grid_dims.item()
    )
    if resize_output_template:
        if is_gchp_restart_file(output_template):
            # This is useful for stretched-grid simulations because they usually
            # don't have a "normal" grid size
            cs_res = np.sqrt(weights.dst_grid_dims.item() / 6).astype(int)

            info_message = (
                "Reshaping the output restart file template to grid size C%f"
            )
            logging.info(info_message, cs_res)

            output_shape = (6 * cs_res, cs_res)
            func = lambda *args, **kwargs: np.ones(output_shape) * np.nan
            vfunc = np.vectorize(func, signature="(lat,lon)->(lat1,lon1)")
            new_output_template = xr.apply_ufunc(
                vfunc,
                output_template,
                keep_attrs=True,
                input_core_dims=[["lat", "lon"]],
                output_core_dims=[["lat1", "lon1"]],
            )
            new_output_template = new_output_template.rename(
                {"lat1": "lat", "lon1": "lon"}
            )
            new_output_template["lat"].attrs = output_template["lat"].attrs
            new_output_template["lon"].attrs = output_template["lat"].attrs
            new_output_template = new_output_template.assign_coords(
                lat=np.arange(
                    new_output_template.dims["lat"], dtype=np.float64
                ),
                lon=np.arange(
                    new_output_template.dims["lon"], dtype=np.float64
                ),
            )
            output_template = new_output_template
        else:
            error_message = (
                "GC-Classic restart resizing not implemented. "
                "Please provide a restart file template with "
                "the proper resolution."
            )
            raise ValueError(error_message)
    else:
        output_shape = output_template_shape

    output_dims = [("lat", "lon"), output_shape]
    logging.info("Regridding the input restart file")
    transform = sparselt.esmf.load_weights(weights, input_dims, output_dims)
    dataset = sparselt.xr.apply(transform, dataset, output_template)
    return dataset


def update_encoding(dataset):
    """
    Ensure dataset variables are encoded as float32.

    Args:
        dataset (xarray.Dataset): The dataset to have its encoding variable
                                  encoding checked and updated.

    Returns:
        xarray.Dataset: The input dataset with float32 variable encoding
                        applied.
    """
    logging.info("Updating encoding")
    for name in dataset.data_vars:
        dataset[name].encoding.update({"dtype": "float32"})
        if (
            "missing_value" in dataset[name].encoding
            and "_FillValue" in dataset[name].encoding
        ):
            del dataset[name].encoding["missing_value"]
    return dataset


def check_for_nans(dataset):
    """
    Check for the presence of NaN values in the passed dataset.

    Args:
        dataset (xarray.Dataset): The dataset to check for NaNs.
    """
    nan_vars = []
    for name in dataset.data_vars:
        if dataset[name].isnull().any().item():
            nan_vars.append(name)
    if len(nan_vars) > 0:
        warning_message = "Dataset has %f variables with NaN values"
        logging.warning(warning_message, len(nan_vars) / len(dataset.data_vars))

        logging.debug("Variables with NaN values: %s", nan_vars)


def regrid_restart_file(
    input_restart,
    regrid_weights,
    output_restart_template,
    stretch_factor=None,
    target_lat=None,
    target_lon=None,
):
    """
    Perform and end-to-end regridding from reading input gridded data and
    regridding weights to writing out the regridded data.

    Args:
        input_restart           (str)  : The path to the restart file that will
                                         be regridded.
        regrid_weights          (str)  : The path to the regridding weights,
                                         generated by ESMF_RegridWeightGen.
        output_restart_template (str)  : The path to the regridding output
                                         template file.
        stretch_factor          (float): An optional stretch factor, for use
                                         with stretched-regridding.
        target_lat              (float): An optional target latitude, for use
                                         with stretched-regridding.
        target_lon              (float): An optional target longitude, for use
                                         with stretched-regridding.
    """
    logging.info("Input restart file: %s", input_restart)
    logging.info("Regridding weights: %s", regrid_weights)
    logging.info("Output template restart file: %s", output_restart_template)

    dataset = open_dataset(input_restart)
    check_for_nans(dataset)
    output_template = open_dataset(output_restart_template)

    input_is_gchp = is_gchp_restart_file(dataset)
    output_is_gchp = is_gchp_restart_file(output_template)

    info_message = "Input restart file type is %s"
    logging.info(info_message, "GCHP" if input_is_gchp else "GC-Classic")

    info_message = "Output restart file type is %s"
    logging.info(info_message, "GCHP" if output_is_gchp else "GC-Classic")

    is_conversion = input_is_gchp != output_is_gchp
    if is_conversion:
        to_gchp = output_is_gchp
        dataset = rename_variables(dataset, to_gchp)
        dataset = reverse_lev(dataset)

    dataset, output_template = drop_variables(dataset, output_template)
    dataset = regrid(dataset, output_template, weights_file=regrid_weights)
    dataset = update_encoding(dataset)
    check_for_nans(dataset)

    if stretch_factor and target_lat and target_lon:
        try:
            dataset.attrs["STRETCH_FACTOR"] = np.float32(stretch_factor)
            dataset.attrs["TARGET_LAT"] = np.float32(target_lat)
            dataset.attrs["TARGET_LON"] = np.float32(target_lon)
        except Exception as exception:
            raise Exception(
                "Error when processing your stretched-grid parameters - are they correct?"
            ) from exception

    dataset.to_netcdf("new_restart_file.nc")

    info_message = "Wrote 'new_restart_file.nc' with %d variables"
    logging.info(info_message, len(dataset.data_vars))

    cleanup_tempfile()


if __name__ == "__main__":
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO").upper())
    COMMAND_LINE = parse_command_line()
    FILE_TO_REGRID = COMMAND_LINE.file_to_regrid
    REGRIDDING_WEIGHTS_FILE = COMMAND_LINE.regridding_weights_file
    TEMPLATE_FILE = COMMAND_LINE.template_file

    if COMMAND_LINE.stretched_grid:
        logging.info("Creating a stretched-grid restart file")

        if (
            (not COMMAND_LINE.stretch_factor)
            or (not COMMAND_LINE.target_latitude)
            or (not COMMAND_LINE.target_longitude)
        ):
            ERROR_MESSAGE = (
                "--stretched-grid was set but not all stretched-"
                "grid parameters were passed!"
            )
            raise RuntimeError(ERROR_MESSAGE)

        STRETCH_FACTOR = COMMAND_LINE.stretch_factor
        TARGET_LATITUDE = COMMAND_LINE.target_latitude
        TARGET_LONGITUDE = COMMAND_LINE.target_longitude

        regrid_restart_file(
            FILE_TO_REGRID,
            REGRIDDING_WEIGHTS_FILE,
            TEMPLATE_FILE,
            stretch_factor=STRETCH_FACTOR,
            target_lat=TARGET_LATITUDE,
            target_lon=TARGET_LONGITUDE,
        )
    else:
        regrid_restart_file(
            FILE_TO_REGRID, REGRIDDING_WEIGHTS_FILE, TEMPLATE_FILE
        )
