'''Unit tests for functions in benchmark.py'''
#aerosols, jvalues, restart, budget, concafterchem, speciesconc, HEMCO

import gcpy.benchmark as bmk
import xarray as xr
import filecmp
import os
from shutil import copyfile

'''
new_sigdiff_files = ['./TestOutput/New_sig_diffs_sfc.txt',
                     './TestOutput/New_sig_diffs_500hpa.txt',
                     './TestOutput/New_sig_diffs_zonalmean.txt',
                     './TestOutput/New_sig_diffs_emissions.txt']
old_sigdiff_files = ['./TestOutput/Correct_sig_diffs_sfc.txt',
                     './TestOutput/Correct_sig_diffs_500hpa.txt',
                     './TestOutput/Correct_sig_diffs_zonalmean.txt',
                     './TestOutput/Correct_sig_diffs_emissions.txt']
'''

def test_compare_single_level():
    refdata = xr.open_dataset('GCCData/SpeciesConcTest.nc4')
    refstr = 'GCC Test Data'
    devdata = xr.open_dataset('GCHPData/SpeciesConcTest.nc4')
    devstr = 'GCHP Test Data'
    pdfname = "TestOutput/NewCompareSingleLevel.pdf"
    try:
        bmk.compare_single_level(refdata, refstr, devdata, devstr, pdfname=pdfname)
        print("compare_single_level ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("compare_single_level failed")
        raise e
    if not os.path.exists('CorrectCompareSingleLevel.pdf'):
        copyfile(pdfname, 'TestOutput/CorrectCompareSingleLevel.pdf')

def test_compare_zonal_mean():
    refdata = xr.open_dataset('GCCData/SpeciesConcTest.nc4')
    refstr = 'GCC Test Data'
    devdata = xr.open_dataset('GCHPData/SpeciesConcTest.nc4')
    devstr = 'GCHP Test Data'
    #pdfname = "TestOutput/NewCompareZonalMean.pdf"
    try:
        bmk.compare_zonal_mean(refdata, refstr, devdata, devstr)#, pdfname=pdfname)
        print("compare_zonal_mean ran successfully, output should have appeared on screen")
    except Exception as e:
        print("compare_zonal_mean failed")
        raise e
    #if not os.path.exists('CorrectCompareZonalMean.pdf'):
    #    copyfile(pdfname, 'TestOutput/CorrectCompareZonalMean.pdf')
    
def test_get_emissions_varnames():
    commonvars = ['EmisCO_Anthro', 'EmisNO_Anthro', 'AREA']
    correct_vars = ['EmisCO_Anthro', 'EmisNO_Anthro']
    assert bmk.get_emissions_varnames(commonvars, 'Emis') == correct_vars
    
def test_create_display_name():
    diag_name = 'EmisCO_Anthro'
    correct_name = "CO Anthro"
    assert bmk.create_display_name(diag_name) == correct_name
            
def test_create_total_emissions_table():
    '''
    #Temporarily commented out because this functionality is included in test_make_benchmark_emis_tables
    refdata = xr.open_dataset('GCCData/EmissionsTest.nc4')
    refstr = 'GCC Test Data'
    devdata = xr.open_dataset('GCHPData/EmissionsTest.nc4')
    devstr = 'GCHP Test Data'
    species = {'CO' : 'Tg', 'ACET' : 'Tg C', 'ALK4' : 'Gg C'}
    outfilename = 'TestOutput/New_Emissions_totals.txt'
    correct_file = 'TestOutput/Original_Emissions_totals.txt'
    bmk.create_total_emissions_table(refdata, refstr, devdata, devstr, species, outfilename)
    if not os.path.exists(correct_file):
        copyfile(outfilename, correct_file)
    assert filecmp.cmp(outfilename, correct_file) == True
    '''
    return
def test_make_benchmark_plots():
    refdata = 'GCCData/SpeciesConcTest.nc4'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/SpeciesConcTest.nc4'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/'
    subdst = 'make_benchmark_plots/'
    try:
        bmk.make_benchmark_plots(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
        print("make_benchmark_plots ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("make_benchmark_plots failed")
        raise e
    
def test_make_benchmark_emis_plots():
    refdata = 'GCCData/EmissionsTest.nc4'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/EmissionsTest.nc4'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/'
    #subdst = 'make_benchmark_emis_plots'
    try:
        bmk.make_benchmark_emis_plots(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)#, subdst=subdst)
        print("make_benchmark_emis_plots ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("make_benchmark_emis_plots failed")
        raise e

def test_make_benchmark_emis_tables():
    refdata = ['GCCData/EmissionsTest.nc4']
    refstr = 'GCC Test Data'
    devdata = ['GCHPData/EmissionsTest.nc4', 'GCHPData/StateMetTest.nc4']
    devstr = 'GCHP Test Data'
    dst = './TestOutput/EmisTables'
    correct_tables = ['./TestOutput/EmisTables/Emissions/OriginalEmissionsTable.txt',
                      './TestOutput/EmisTables/Emissions/OriginalInventoryTable.txt']
    outfilenames = ['./TestOutput/EmisTables/Emissions/Emission_totals.txt',
                    './TestOutput/EmisTables/Emissions/Inventory_totals.txt']
    bmk.make_benchmark_emis_tables(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
    for i in [0, 1]:
        if not os.path.exists(correct_tables[i]):
            copyfile(outfilenames[i], correct_tables[i])
        assert filecmp.cmp(outfilenames[i], correct_tables[i]) == True

def test_make_benchmark_jvalue_plots():
    refdata = 'GCCData/JValuesTest.nc4'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/JValuesTest.nc4'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/'
    subdst = 'make_benchmark_jvalue_plots/'
    try:
        bmk.make_benchmark_jvalue_plots(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
        print("make_benchmark_jvalue_plots ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("make_benchmark_jvalue_plots failed")
        raise e

def test_make_benchmark_aod_plots():
    refdata = 'GCCData/AerosolsTest.nc4'
    refstr = 'GCC Test Data'
    devdata = 'GCHPData/AerosolsTest.nc4'
    devstr = 'GCHP Test Data'
    dst = './TestOutput/'
    #subdst = 'make_benchmark_aod_plots/'
    try:
        bmk.make_benchmark_aod_plots(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
        print("make_benchmark_aod_plots ran successfully, check output PDFs  to visually verify")
    except Exception as e:
        print("make_benchmark_aod_plots failed")
        raise e

def test_make_benchmark_mass_tables():
    refdata = 'GCCData/RestartTest.nc4'
    refstr = 'GCC_Test_Data'
    devdata = 'GCCData/RestartTest.nc4'
    devstr = 'GCC_Test_Data'
    dst = './TestOutput/MassTables'
    correct_tables = ['./TestOutput/MassTables/CorrectGlobalMass_TropStrat.txt',
                      './TestOutput/MassTables/CorrectGlobalMass_Trop.txt']
    outfilenames = ['./TestOutput/MassTables/GCC_Test_Data_GlobalMass_TropStrat.txt',
                    './TestOutput/MassTables/GCC_Test_Data_GlobalMass_Trop.txt']    
    bmk.make_benchmark_mass_tables(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
    for i in [0, 1]:
        if not os.path.exists(correct_tables[i]):
            copyfile(outfilenames[i], correct_tables[i])
        assert filecmp.cmp(outfilenames[i], correct_tables[i]) == True

    
def test_make_benchmark_budget_tables():
    devdata = 'GCCData/BudgetTest.nc4'
    devstr = 'GCC_Test_Data'
    dst = './TestOutput/Budget/'
    correct_table = './TestOutput/Budget/OriginalBudgetTable.txt'
    output_table = ''
    correct_tables = ['./TestOutput/Budget/Budget/CorrectBudget_Full.txt',
                      './TestOutput/Budget/Budget/CorrectBudget_gPBL.txt',
                      './TestOutput/Budget/Budget/CorrectBudget_nPBL.txt',
                      './TestOutput/Budget/Budget/CorrectBudget_pPBL.txt',
                      './TestOutput/Budget/Budget/CorrectBudget_tPBL.txt',
                      './TestOutput/Budget/Budget/CorrectBudget_Trop.txt',
                      './TestOutput/Budget/Budget/CorrectBudget_yPBL.txt']
    outfilenames = ['./TestOutput/Budget/Budget/Budget_Full.txt',
                    './TestOutput/Budget/Budget/Budget_gPBL.txt',
                    './TestOutput/Budget/Budget/Budget_nPBL.txt',
                    './TestOutput/Budget/Budget/Budget_pPBL.txt',
                    './TestOutput/Budget/Budget/Budget_tPBL.txt',
                    './TestOutput/Budget/Budget/Budget_Trop.txt',
                    './TestOutput/Budget/Budget/Budget_yPBL.txt']
    bmk.make_benchmark_budget_tables(devdata, devstr, dst=dst, overwrite=True)
    for i in [0, 1, 2, 3, 4, 5, 6]:
        if not os.path.exists(correct_tables[i]):
            copyfile(outfilenames[i], correct_tables[i])
        assert filecmp.cmp(outfilenames[i], correct_tables[i]) == True

    
def test_make_benchmark_oh_metrics():
    refdata = ['./GCCData/ConcAfterChemTest.nc4', './GCCData/StateMetTest.nc4']
    refstr = 'GCC_Test_Data'
    devdata = ['./GCCData/ConcAfterChemTest.nc4', './GCCData/StateMetTest.nc4']
    devstr = 'GCC_Test_Data'
    dst = './TestOutput/OHMetrics/'
    correct_file = './TestOutput/OHMetrics/Tables/OriginalOH.txt'
    output_file = './TestOutput/OHMetrics/Tables/GCC_Test_Data_OH_metrics.txt'
    bmk.make_benchmark_oh_metrics(refdata, refstr, devdata, devstr, dst=dst, overwrite=True)
    if not os.path.exists(correct_file):
        copyfile(output_file, correct_file)

    assert filecmp.cmp(correct_file, output_file) == True

def test_add_bookmarks_to_pdf():
    #use manual confirmation
    return
def test_add_nested_bookmarks_to_pdf():
    #use manual confirmation
    return
def test_add_missing_variables():
    refdata = './GCCData/SpeciesConcTest.nc4'
    refdata = xr.open_dataset(refdata)
    #Should change this to real, different datasets
    devdata = './GCCData/SpeciesConcTest.nc4'
    devdata = xr.open_dataset(devdata)
    newref, newdev = bmk.add_missing_variables(refdata, devdata)
    assert refdata.equals(devdata)
def test_get_troposphere_mask():
    return
def test_check_units():
    refdata = './GCCData/SpeciesConcTest.nc4'
    refdata = xr.open_dataset(refdata)
    devdata = './GCHPData/SpeciesConcTest.nc4'
    devdata = xr.open_dataset(devdata)
    refunits, devunits = bmk.check_units(refdata, devdata, 'SpeciesConc_CH4')
    assert refunits == devunits
def test_reshape_MAPL_CS():
    return

def main():
    test_get_emissions_varnames()
    test_create_display_name()
    test_create_total_emissions_table()
    test_make_benchmark_plots()
    test_make_benchmark_emis_plots()
    test_make_benchmark_emis_tables()
    test_make_benchmark_jvalue_plots()
    test_make_benchmark_aod_plots()
    test_make_benchmark_mass_tables()
    test_make_benchmark_budget_tables()
    test_make_benchmark_oh_metrics()
    test_add_bookmarks_to_pdf()
    test_add_nested_bookmarks_to_pdf()
    test_add_missing_variables()
    test_get_troposphere_mask()
    test_check_units()
    test_reshape_MAPL_CS()
    test_compare_single_level()        
    test_compare_zonal_mean()    
if __name__ == "__main__":
    main()
