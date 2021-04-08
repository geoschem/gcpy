# README for examples/diagnostics

## Author
Bob Yantosca (GCST), 08 Apr 2021

## Files:

  1. `compare_diags.py`: Script to compare diagnostic output between 2 GEOS-Chem model runs ("Ref" and "Dev").  Useful for debugging.

  2. `compare_diags.yml`: YAML file containing user-modifiable options for compare_diags.py.

  3. `compare_diagnostics.ipynb`: Jupyter Notebook that walks you through comparing diagnostic outputs.  (This is probably a little dated, we would recommend that you us the compare.diags.py script).

### Using compare_diags.py:

  1. We recommend that you make a copy of the `compare_diags.yml` file and then add your own modifications into that.  You can name it e.g. `my_compare_diags.yml`.  The YAML field names are more or less self-explanatory.

  2. Run `compare_diags.py` with this command:

    $ ./compare_diags.py my_compare_diags.yml
