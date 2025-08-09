# Installation as a developer

Developers should clone the repository and install in a Python virtual or conda environment:

```shell
git clone https://www.github.com/atomprobe-tc/ifes_apt_tc_data_modeling.git
cd ifes_apt_tc_data_modeling
python -m pip install --upgrade pip
python -m pip install -e ".[dev,docs,ipynb]"
```

This will install modules for linting, code styling, and unit testing via the [pytest](https://docs.pytest.org/en/stable) framework.
Unit tests can then be started from the root directory of the installation with the following call:

```shell
pytest -sv tests
```
