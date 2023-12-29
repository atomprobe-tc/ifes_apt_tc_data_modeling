#!/bin/bash

# execute within a local python virtual environment or within a conda environment
# in both cases the ifes_apt_tc_data_modeling module should have been installed in developer mode

python -m pycodestyle --ignore=E501 ifes_apt_tc_data_modeling
python -m pylint ifes_apt_tc_data_modeling --ignore build
python -m mypy --ignore-missing-imports ifes_apt_tc_data_modeling