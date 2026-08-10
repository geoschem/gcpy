# Environment files

This folder contains environment files that are used to install dependencies for GCPy and the ReadTheDocs documentation. 

## Installing GCPy dependencies

### With Mamba or Conda

Use one of these commands to build a Mamba/Conda environment with all of the GCPy dependencies.

```console
$ mamba env create -n gcpy_env --file=gcpy_environment_py312.yml   # If you wish to use Python 3.12
$ mamba env create -n gcpy_env --file=gcpy_environment_py313.yml   # If you wish to use Python 3.13
```
or
```console
$ conda env create -n gcpy_env --file=gcpy_environment_py312.yml   # If you wish to use Python 3.12
$ conda env create -n gcpy_env --file=gcpy_environment_py313.yml   # If you wish to use Python 3.13
```

## Installing ReadTheDocs dependencies

### With Mamba or Conda

The Sphinx/ReadTheDocs packages are now bundled into the GCPy environment
files, so there is no separate environment to create. Just build the
`gcpy_env` environment as described above
(`gcpy_environment_py312.yml`/`gcpy_environment_py313.yml`) and it will
already contain everything needed to build the documentation.

### With pip
Or use this command to install the dependencies from PyPi:
```console
$ pip install -r read_the_docs_requirements.txt
```
