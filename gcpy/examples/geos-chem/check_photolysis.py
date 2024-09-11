#!/usr/bin/env python3
"""
Convenience script to verify that photolysis reactions are specified
properly in GEOS-Chem.
"""

from gcpy.util import read_config_file, verify_variable_type

def get_photo_species_from_species_db(path):
    """
    Returns a list of species that have the "Is_Photolysis: true"
    setting in the species_database.yml file.
   
    Args
    path    : str : Absolute path to the species_database.yml file

    Returns
    species : dict : Names of species having "Is_Photolysis: true"
    """
    verify_variable_type(path, str)

    species_database = read_config_file(path)
    species = {}
    for spc in species_database:
        if "Is_Photolysis" in species_database[spc]:
            species[spc] = bool(species_database[spc]["Is_Photolysis"])

    return species


def get_photo_species_from_fjx_j2j(path):
    """
    Returns a list of photolysis species in the FJX_j2j.dat file.

    Args
    path    : str  : Absolute path to the FJX_j2j.dat file

    Returns
    species : list : Tuples containing species name and index
    """
    verify_variable_type(path, str)

    species = []
    indices = []
    with open(path, "r", encoding="utf-8") as ifile:
        for line in list(ifile):
            line = line.strip("\n")
            if "PHOTON" in line:
                sub_string = line.split()
                indices.append(int(sub_string[0]))
                species.append(sub_string[1])

    return list(zip(species, indices))


def get_photo_species_from_kpp_eqn(path):
    """
    Returns a list of species having photolysis reactions in
    the KPP fullchem.eqn dat file.

    Args
    path    : str  : Absolute path to the fullchem.eqn file

    Returns
    species : dict : Species name plus index of the PHOTOL array
    """
    verify_variable_type(path, str)
    
    species = []
    photol = []
    with open(path, "r", encoding="utf-8") as ifile:
        for line in list(ifile):
            line = line.strip("\n")
            name = ""
            if "+ hv" in line:
                species.append(line.split()[0])
            if "PHOTOL(" in line:
                photol.append(int(
                    line.split(":")[1].split(";")[0]\
                    .strip().strip("PHOTOL(").strip(")")
                ))
                
    return list(zip(species, photol))


def get_max_index(multi_list):
    """
    Computes the maximum index value in a "multi-index" list
    (i.e. a list of tuples).

    Args
    multi_list : list : Tuples containing species name and index

    Returns
    max_index  : int  : Maximum index value
    """
    max_index = 0
    for (species, index) in multi_list:
        if index > max_index:
            max_index = index
            
    return index


def print_results(spc_db, fjx_j2j, kpp_eqn):
    """
    Prints results 

    Args
    spc_db  : dict : Photolyzing species listed in species_database.yml
    fjx_j2j : dict : Photolyzing species listed in FJX_j2j.dat
    kpp_eqn : dict : Photolyzing species listed in fullchem.eqn
    """

    # Maximum number of entries in FJX_j2j.dat
    max_j2j = get_max_index(fjx_j2j)
    max_kpp = get_max_index(kpp_eqn)
    print(f"Number of entries in FJX_j2j.dat     : {max_j2j:<3}")
    print(f"Number of photo rxns in fullchem.eqn : {max_kpp:<3}")
    if max_kpp != max_j2j:
        msg = "ERROR: Mismatch between FJX_j2j.dat and fullchem.eqn!"
        raise ValueError(msg)
    
    # Print photolyzing species and the corresponding index
    # of the KPP PHOTOL array that it uses.  Ths should be the
    # same order as in the FJX_j2j.dat file.
    for (species, index) in kpp_eqn:
        if index > 0 and index <= max_j2j:
            msg = f"{species:<10} uses PHOTOL({index:<3}) "
            msg += f"and has Is_Photolysis={spc_db[species]}"
            print(msg)
        else:
            msg = "{species:<10} has an invalid PHOTOL index "
            msg += "({index}) > {max_j2j}!"
            print(msg)


def main():
    """
    """
    
    # Get all species with "Is_Photolysis: true" in the species database
    path = "/n/holyscratch01/jacob_lab/ryantosca/tests/14.5.0/gc2457/GCClassic/src/GEOS-Chem/run/shared/species_database.yml"
    spc_db = get_photo_species_from_species_db(path)

    # Get all species listed in FJX_j2j.dat
    path = "/n/holyscratch01/external_repos/GEOS-CHEM/gcgrid/data/ExtData/CHEM_INPUTS/CLOUD_J/v2024-09/FJX_j2j.dat"
    fjx_j2j = get_photo_species_from_fjx_j2j(path)

    # Get all species with a photo reaction in the f
    path = "/n/holyscratch01/jacob_lab/ryantosca/tests/14.5.0/gc2457/GCClassic/src/GEOS-Chem/KPP/fullchem/fullchem.eqn"
    kpp_eqn = get_photo_species_from_kpp_eqn(path)

    # Print results
    print_results(spc_db, fjx_j2j, kpp_eqn)
    

if __name__ == '__main__':
    main()
