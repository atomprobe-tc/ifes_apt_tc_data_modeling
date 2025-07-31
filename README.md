# ifes_apt_tc_data_modeling

## Context:
`ifes_apt_tc_data_modeling` is a collection of Python modules for reading file formats
that are used in the research field of atom probe tomography and field ion microscopy.

The project combined all publicly available efforts of the atom probe research community
that individuals and AMETEK/Cameca as the leading technology partner in the field have
exchanged about how to work with different file formats that have been used in the field.

## Getting started:
To use this library create either a [Python virtual environment](https://docs.python.org/3/library/venv.html) or [conda environment](https://www.anaconda.com/docs/tools/working-with-conda/environments). We recommend using Python 3.12.

### Installation as a user
Decide or create a directory on your local computer or server where you would like to
install the software. Then follow these steps:

```
git clone https://www.github.com/atomprobe-tc/ifes_apt_tc_data_modeling.git
cd ifes_apt_tc_data_modeling
python -m pip install --upgrade pip
python -m pip install -e ".[ipynb]"
```
This will install the software and jupyterlab whereby the notebook with examples become executable.
The jupyterlab server can be started with

```
jupyter-lab
```

The notebook to run is the following `examples/ExamplesForUsersOrDevelopers.ipynb`

### Installation as a software developer
Developers should change the installation call line with the following:

```
python -m pip install -e ".[dev,docs,ipynb]"
```

This will install modules for linting, code styling and performing unit testing via the pytest framework.
Unit tests can then be started from the root directory of the installation with the following call:

```
pytest -sv tests
```

[Further documentation of the software is available here.]()

##Acknowledgements:
Contributions in this regard of individuals like Daniel Haley, Andrew London, Baptiste Gault,
David Reinhardt, Jim Payne, Andrew Breen, and Benjamin Caplins, Peter Felfer,
Mehrpad Monajem, and Martina Heller are noteworthy to mention here.

The work is supported by the International Field Emission Society (IFES). These efforts were
consolidated by Markus Kühbach who wrote and maintains this library, made possible thanks
to funding from the FAIRmat project that is a part of the German National Research Data Infrastructure.
