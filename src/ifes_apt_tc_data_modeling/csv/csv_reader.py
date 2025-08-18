#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""CSV file format reader as sometimes used for reporting POS file content."""

# pylint: disable=duplicate-code

import os

import numpy as np
import pandas as pd

from ifes_apt_tc_data_modeling.nexus.nx_field import NxField


class ReadCsvFileFormat:
    """Read CSV file assuming (n_ions, 4) like in POS."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 4) or not file_path.lower().endswith(".csv"):
            raise ImportError(
                "WARNING::CSV file incorrect file_path ending or file type."
            )
        self.file_path = file_path

        self.file_size = os.path.getsize(self.file_path)
        self.number_of_events = None
        self.version = None

        shp = np.shape(pd.read_csv(self.file_path))
        if shp[0] > 0 and shp[1] == 4:
            self.number_of_events = shp[0]
        else:
            raise ImportError(
                "CSV file unsupported version because not formatted like POS."
            )

    def get_reconstructed_positions(self):
        """Read xyz columns."""

        xyz = NxField()
        xyz.values = np.zeros([self.number_of_events, 3], np.float32)
        xyz.unit = "nm"
        # there are too many assumption made here as to the content
        # in the csv file sure one could pass some configuration hints but
        # frankly there are nowadays much! better strategies to report
        # atom probe data than CSV, NeXus is one such, also csv files have
        # no magic number, de facto this works only because users know what
        # to expect in advance but how should a machine know this?
        for dim in [0, 1, 2]:
            xyz.values[:, dim] = pd.read_csv(self.file_path).iloc[:, dim]
        return xyz

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge-state-ratio column."""

        m_n = NxField()
        m_n.values = np.zeros([self.number_of_events, 1], np.float32)
        m_n.unit = "Da"
        # again such a strong assumption!
        # why reported in Da?
        # why in the third column
        # why at all a mass-to-charge-state-ratio value array?
        m_n.values[:, 0] = pd.read_csv(self.file_path).iloc[:, 3]
        return m_n
