# Installation as a user

It is recommended to use python 3.13 with a dedicated virtual environment for this package. Learn how to manage [python versions](https://github.com/pyenv/pyenv) and [virtual environments](https://realpython.com/python-virtual-environments-a-primer/).

There are many alternatives to managing virtual environments and package dependencies (requirements). We recommend using [`uv`](https://github.com/astral-sh/uv), an extremely fast manager Python package and project manager. In this tutorial, you will find paralleled descriptions, using either `uv` or a more classical approach using `venv` and `pip`. Using a [conda environment](https://www.anaconda.com/docs/tools/working-with-conda/environments) is yet another alternative.

The library can be used on Windows, Mac, or Unix, provided a Python version is installed.

## Setup

Start by creating a virtual environment, e.g., in a directory on your local computer:

=== "uv"
    `uv` is capable of creating a virtual environment and install the required Python version at the same time.

    ```bash
    uv venv --python 3.13
    ```

=== "venv"
    Note that you will need to install the Python version manually beforehand.

    ```bash
    python -m venv .venv
    source .venv/bin/activate
    ```
    
=== "conda"
    Note that you will need to install the Python version manually beforehand.
    
    ```bash
    conda create -n venv python=3.13
    conda activate venv
    conda install pip
    ```
That command creates a new virtual environment in a directory called `.venv`.

## Installation

Install the latest stable version of this package from PyPI with


=== "uv"

    ```bash
    uv pip install ifes_apt_tc_data_modeling[ipynb]
    ```

=== "pip"


    ```bash
    pip install ifes_apt_tc_data_modeling[ipynb]
    ```

=== "conda"

    ```bash
    python -m pip install ifes_apt_tc_data_modeling[ipynb]
    ```

This will install the module and [jupyterlab](https://jupyterlab.readthedocs.io/en/latest/) whereby the notebook
with examples become executable. 

## Start using `ifes_apt_tc_data_modeling`

The jupyterlab server is started with

```bash
jupyter-lab
```

The notebook to run is the following `examples/ExamplesForUsersOrDevelopers.ipynb`
That's it! You can now use the library that you have installed!



