#!/usr/bin/env python3
"""
Creates the "HEMCO_sa_Spec.rc" file (needed for the HEMCO
standalone model) from a "geoschem_species_metadata.yml" file taken
from a GEOS-Chem simulation.

Calling sequence:
$ ./make_hemco_sa_spec.py geoschem_species_metadata.yml
"""
import sys
from gcpy.util import read_config_file, verify_variable_type


def write_to_file(metadata):
    """
    Writes species metadata to the "HEMCO_sa_Spec.rc" file for the
    HEMCO standaone mode.  Output includes species index, name, MW (g),
    and Henry's law K0, CR, PKA parameters.
    """
    verify_variable_type(metadata, dict)

    with open( "HEMCO_sa_Spec.rc", "w", encoding="utf-8") as ofile:

        # Write header
        header = """
# List species below. For each species, the following entries must be given:
# ID        : species ID
# NAME      : species name
# MW        : molecular weight (g/mol)
# K0        : Henry constant at 298 K (M/atm, liquid over gas)
# CR        : Henry temperature dependency (-d ln(KH)/d(1/T)
# pKA       : acid dissociation constant
#
# NOTE: These are the species corresponding to the standard simulation!
#
#ID NAME           MW      K0                CR        PKA
"""
        print(header, file=ofile)

        # Loop over species metadata
        for idx, var in enumerate(metadata):

            # Get Henry's law parameters (or set to 0 if not defined)
            henry_k0 = 0.0
            if "Henry_K0" in metadata[var]:
                henry_k0 = float(metadata[var]['Henry_K0'])
            henry_cr = 0.0
            if "Henry_CR" in metadata[var]:
                henry_cr = float(metadata[var]['Henry_CR'])
            henry_pka = 0.0
            if "Henry_pKa" in metadata[var]:
                henry_pka = float(metadata[var]['Henry_pKa'])

            # Write to file
            line=f"{idx+1:<3d} {var:11}  "
            line+=f"{metadata[var]['MW_g']:7.2f}  "
            line+=f"{henry_k0:13.6e} {henry_cr:9.2f} {henry_pka:9.2f}"
            print(line, file=ofile)


def make_hemco_sa_spec(argv):
    """
    Reads metadata from the "geoschem_species_metadata.yml" file
    and then calls the "write_to_file" routine to create the
    "HEMCO_sa_Spec.rc" configuration file.
    """
    if len(argv) != 2:
        msg = "The path to geoschem_species_metadata.yml was not passed!"
        raise FileNotFoundError(msg)

    try:
        metadata = read_config_file(sys.argv[1])
    except FileNotFoundError as exc:
        msg = "Could not find the 'geoschem_species_metadata.yml' file!"
        raise FileNotFoundError(msg) from exc

    write_to_file(metadata)


if __name__ == '__main__':
    make_hemco_sa_spec(sys.argv)
