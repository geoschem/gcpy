# Environment files

This folder contains environment files that are used to install dependencies for GCPy and the ReadTheDocs documentation.

NOTE: Most users will install GCPy with:

```console
$ pip install geoschem-gcpy
```
but GCPy developers may need to install the dependencies separately using the environment files in this folder.

## Installing GCPy dependencies

### With Mamba or Conda

Use one of these commands to build a Mamba/Conda environment with all of the GCPy dependencies.

```console
$ mamba env create -n gcpy_env --file=gcpy_environment.yml
```
or
```console
$ conda env create -n gcpy_env --file=gcpy_environment.yml
```

### With pip
Or use this command to install the dependencies from PyPI:
```console
$ pip install -r gcpy_requirements.txt
```

## Installing ReadTheDocs dependencies

### With Mamba or Conda

Use one of these commands to build a Mamba/Conda environment with the dependencies for building the GCPy ReadTheDocs documentation:

```console
$ mamba env create -n rtd_env --file=read_the_docs_environment.yml
```
or
```console
$ conda env create -n rtd_env --file=read_the_docs_environment.yml
```

### With pip
Or use this command to install the dependencies from PyPi:
```console
$ pip install -r read_the_docs_requirements.txt
```
