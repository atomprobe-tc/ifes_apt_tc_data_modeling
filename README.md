# ifes_apt_tc_data_modeling

## Context:
`ifes_apt_tc_data_modeling` is a collection of Python modules for reading file formats
that are used in the research field of atom probe tomography and field ion microscopy.

The project combined all publicly available efforts of the atom probe research community
that individuals and [AMETEK/Cameca](https://github.com/CamecaAPT) as the leading technology partner in the field
have exchanged about how to work with different file formats used for atom probe research.

## Getting started:
To use this library create either a [Python virtual environment](https://docs.python.org/3/library/venv.html) or [conda environment](https://www.anaconda.com/docs/tools/working-with-conda/environments).
### Installation as a user
The library can be used on Windows, Mac, or Unix, provided a Python version is installed.
Exemplified for using a Python virtual environment and Python 3.12, firstly you should create a
virtual environment — if you do not have one yet — in a directory on your local computer or server.

```shell
python3 -m venv .py3.12
source .py3.12/bin/activate
```

If you already have such environment or just created one, proceed with installing `ifes_apt_tc_data_modeling`:

```shell
pip install ifes_apt_tc_data_modeling[ipynb]
```

This will install the module and [jupyterlab](https://jupyterlab.readthedocs.io/en/latest/) whereby the notebook
with examples become executable. The jupyterlab server is started with

```shell
jupyter-lab
```

The notebook to run is the following `examples/ExamplesForUsersOrDevelopers.ipynb`

[Further documentation of the software is available here.]()

## Acknowledgements:
Contributions of individuals made the consolidation of open-source software for parsing atom probe data possible. Here, [Daniel Haley](https://orcid.org/0000-0001-9308-2620), [Andrew London](https://orcid.org/0000-0001-6959-9849), [Baptiste Gault](https://orcid.org/0000-0002-4934-0458), David Reinhard, Jim Payne, [Andrew Breen](https://orcid.org/0000-0002-3600-5108), and [Benjamin Caplins](https://orcid.org/0000-0002-4925-9537), [Peter Felfer](https://orcid.org/0000-0002-2338-1016), and [Mehrpad Monajem](https://orcid.org/0009-0002-6746-2835) are mentioned.

The work on the ifes_apt_tc_data_modeling library is supported by the International Field Emission Society (IFES).
The library was written and is maintained by [Markus Kühbach](0000-0002-7117-5196). The work is funded by the
Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - 460197019 (FAIRmat).
