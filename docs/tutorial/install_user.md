# Installation as a user

To use this library create either a [Python virtual environment](https://docs.python.org/3/library/venv.html) or [conda environment](https://www.anaconda.com/docs/tools/working-with-conda/environments).

The library can be used on Windows, Mac, or Unix, provided a Python version is installed.
Exemplified for using a Python virtual environment and Python 3.12, firstly you should create a
virtual environment — if you do not have one yet — in a directory on your local computer or server.

```shell
python3 -m venv .py3.12
source .py3.12/bin/activate
```

If you already have such environment or just created one, proceed with installing the library:

```shell
pip install ifes_apt_tc_data_modeling[ipynb]
```

This will install the module and [jupyterlab](https://jupyterlab.readthedocs.io/en/latest/) whereby the notebook
with examples become executable. The jupyterlab server is started with.

```shell
jupyter-lab
```

The notebook to run is the following `examples/ExamplesForUsersOrDevelopers.ipynb`
