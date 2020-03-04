'''Unit tests for functions in benchmark.py'''
#aerosols, jvalues, restart, budget, concafterchem, speciesconc, HEMCO

import gcpy.benchmark as bmk
import xarray as xr
import filecmp

def test_compare_single_level():
    refdata = xr.open_dataset('GCCData/SpeciesConcTest.nc')
    refstr = 'GCC Test Data'
    devdata = xr.open_dataset('GCHPData/SpeciesConcTest.nc')
    devstr = 'GCHP Test Data'
    pdfname = "TestOutput/NewCompareSingleLevel.pdf"
    try:
        bmk.compare_single_level(refdata, refstr, devdata, devstr, pdfname=pdfname)
        print("compare_single_level ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("compare_single_level failed")
        raise e

def test_compare_zonal_mean():
    refdata = xr.open_dataset('GCCData/SpeciesConcTest.nc')
    refstr = 'GCC Test Data'
    devdata = xr.open_dataset('GCHPData/SpeciesConcTest.nc')
    devstr = 'GCHP Test Data'
    pdfname = "TestOutput/NewCompareSingleLevel.pdf"
    try:
        bmk.compare_zonal_mean(refdata, refstr, devdata, devstr, pdfname=pdfname)
        print("compare_zonal_mean ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("compare_zonal_mean failed")
        raise e
    
def test_get_emissions_varnames():
    commonvars = ['EmisCO_Anthro', 'EmisNO_Anthro', 'AREA']
    correct_vars = ['EmisCO_Antrho', 'EmisNO_Anthro']
    assert bmk.get_emissions_varnames(commonvars) == correct_vars
    
def test_create_display_name():
    diag_name = 'EmisCO_Anthro'
    correct_name = "CO Anthro"
    assert bmk.create_display_name(diag_name) = correct_name
            
def test_create_total_emissions_table():
    refdata = xr.open_dataset('GCCData/HEMCODiagTest.nc')
    refstr = 'GCC Test Data'
    devdata = xr.open_dataset('GCHPData/HEMCODiagTest.nc')
    devstr = 'GCHP Test Data'
    species = {'CO' : 'Tg', 'ACET' : 'Tg C', 'ALK4' : 'Gg C'}
    outfilename = 'TestOutput/New_Emissions_totals.txt'
    correct_file = 'TestOutput/Original_Emissions_totals.txt'
    bmk.create_total_emissions_table(refdata, refstr, devdata, devstr, species, outfilename)
    assert filecmp.cmp(outfilename, correct_file) == True

def test_make_benchmark_plots():
    refdata = 'GCCData/SpeciesConcTest.nc'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/SpeciesConcTest.nc'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/'
    subdst = './TestOutput/make_benchmark_plots'
    try:
        bmk.make_benchmark_plots(refdata, refstr, devdata, devstr, subdst=subdst, overwrite=True)
        print("make_benchmark_plots ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("make_benchmark_plots failed")
        raise e
    
def test_make_benchmark_emis_plots():
    refdata = 'GCCData/HEMCODiagTest.nc'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/HEMCODiagTest.nc'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/'
    subdst = './TestOutput/make_benchmark_emis_plots'
    try:
        bmk.make_benchmark_emis_plots(refdata, refstr, devdata, devstr, dst=dst, subdst=subdst, overwrite=True)
        print("make_benchmark_emis_plots ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("make_benchmark_emis_plots failed")
        raise e

def test_make_benchmark_emis_tables():
    refdata = 'GCCData/HEMCODiagTest.nc'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/HEMCODiagTest.nc'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/EmisTables'
    correct_table = './TestOutput/EmisTables/OriginalEmissionsTable.txt'
    outfilename = ''
    bmk.make_benchmark_emis_tables(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
    assert filecmp.cmp(output_table, correct_table) == True

def test_make_benchmark_jvalue_plots():
    refdata = 'GCCData/JValuesTest.nc'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/JValuesTest.nc'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/'
    subdst = './TestOutput/make_benchmark_jvalue_plots'
    try:
        bmk.make_benchmark_jvalue_plots(refdata, refstr, devdata, devstr, subdst=subdst, overwrite=True)
        print("make_benchmark_jvalue_plots ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("make_benchmark_jvalue_plots failed")
        raise e

def test_make_benchmark_aod_plots():
    refdata = 'GCCData/AerosolsTest.nc'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/AerosolsTest.nc'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/'
    subdst = './TestOutput/make_benchmark_aod_plots'
    try:
        bmk.make_benchmark_plots(refdata, refstr, devdata, devstr, subdst=subdst, overwrite=True)
        print("make_benchmark_aod_plots ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("make_benchmark_aod_plots failed")
        raise e

def test_make_benchmark_mass_tables():
    refdata = 'GCCData/RestartTest.nc'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/RestartTest.nc'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/MassTables'
    correct_table = './TestOutput/OriginalMassTable.txt'
    output_table = ''
    bmk.make_benchmark_emis_tables(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
    assert filecmp.cmp(output_table, correct_table) == True
    
def test_make_benchmark_budget_tables():
    devdata = 'GCCData/BudgetTest.nc'
    devstr = 'GCC Test Data'
    dst = './TestOutput/Budget/'
    correct_table = './TestOutput/Budget/OriginalBudgetTable.txt'
    output_table =''
    bmk.make_benchmark_budget_tables(devdata, devstr, dst=dst, overwrite=True)
    assert filecmp.cmp(output_table, correct_table) == True
    
def test_make_benchmark_oh_metrics():
    refdata = ['./GCCData/ConcAfterChemTest.nc', './GCCData/MetTest.nc']
    refstr = 'GCC Test Data'
    devdata = ['./GCHPData/ConcAfterChemTest.nc', './GCHPData/MetTest.nc']
    devstr = 'GCHP Test Data'
    dst = './TestOutput/OHMetrics/'
    correct_file = './TestOutput/OHMetrics/OriginalOH.txt'
    output_file = ''
    bmk.make_benchmark_oh_metrics(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
    assert filecmp.cmp(correct_file, output_file) == True

def test_add_bookmarks_to_pdf():
    #use manual confirmation
    return
def test_add_nested_bookmarks_to_pdf():
    #use manual confirmation
    return
def test_add_missing_variables():
    refdata = './GCCData/SpeciesConcTest.nc'
    devdata = './GCHPData/HEMCODiagTest.nc'
    newref, newdev = bmk.add_missing_variables(refdata, devdata)
    assert refdata.equals(devdata)
def test_get_troposphere_mask():
    return
def test_check_units():
    refdata = './GCCData/SpeciesConcTest.nc'
    devdata = './GCHPData/SpeciesConcTest.nc'
    refunits, devunits = bmk.check_units(refdata, vardata, 'varnamehere')
    assert refunits == devunits
def test_reshape_MAPL_CS():
    return

def main():
    test_compare_single_level():        
    test_compare_zonal_mean():
    test_get_emissions_varnames():
    test_create_display_name():
    test_print_totals():
    test_create_total_emissions_table():
    test_create_global_mass_table():
    test_create_budget_table():
    test_get_species_categories():
    test_archive_species_categories():
    test_make_benchmark_plots():
    test_make_benchmark_emis_plots():
    test_make_benchmark_emis_tables():
    test_make_benchmark_jvalue_plots():
    test_make_benchmark_aod_plots():
    test_make_benchmark_mass_tables():
    test_make_benchmark_budget_tabels():
    test_make_benchmark_oh_metrics():
    test_add_bookmarks_to_pdf():
    test_add_nested_bookmarks_to_pdf():
    test_add_missing_variables():
    test_get_troposphere_mask():
    test_check_units():
    test_reshape_MAPL_CS():
    

