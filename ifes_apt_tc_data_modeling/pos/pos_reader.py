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

"""POS file format reader used by atom probe microscopists."""

# pylint: disable=no-member,duplicate-code

import os
import numpy as np

from ifes_apt_tc_data_modeling.nexus.nx_field import NxField
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data


class ReadPosFileFormat():
    """Read *.pos file format."""

    def __init__(self, filename: str):
        """Initialize the reader."""
        if (len(filename) <= 4) or (filename.lower().endswith(".pos") is False):
            raise ImportError("WARNING::POS file incorrect filename ending or file type!")
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)
        assert self.filesize % 4 * 4 == 0, \
            "POS filesize not integer multiple of 4*4B!"
        assert np.uint32(self.filesize / (4 * 4)) < np.iinfo(np.uint32).max, \
            "POS file is too large, currently only 2*32 supported!"
        self.number_of_events = np.uint32(self.filesize / (4 * 4))
        # print("Initialized access to " + self.filename + " successfully")

        # https://doi.org/10.1007/978-1-4614-3436-8 for file format details
        # dtyp_names = ["Reconstructed position along the x-axis (nm)",
        #               "Reconstructed position along the y-axis (nm)",
        #               "Reconstructed position along the z-axis (nm)",
        #               "Reconstructed mass-to-charge-state ratio (Da)"]

    def get_reconstructed_positions(self):
        """Read xyz columns."""

        xyz = NxField()
        xyz.typed_value = np.zeros(
            [self.number_of_events, 3], np.float32)
        xyz.unit = "nm"

        xyz.typed_value[:, 0] = \
            get_memory_mapped_data(self.filename, ">f4",
                                   0 * 4, 4 * 4, self.number_of_events)  # x
        xyz.typed_value[:, 1] = \
            get_memory_mapped_data(self.filename, ">f4",
                                   1 * 4, 4 * 4, self.number_of_events)  # y
        xyz.typed_value[:, 2] = \
            get_memory_mapped_data(self.filename, ">f4",
                                   2 * 4, 4 * 4, self.number_of_events)  # z
        return xyz

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge-state-ratio column."""

        m_n = NxField()
        m_n.typed_value = np.zeros(
            [self.number_of_events, 1], np.float32)
        m_n.unit = "Da"

        m_n.typed_value[:, 0] = \
            get_memory_mapped_data(self.filename, ">f4",
                                   3 * 4, 4 * 4, self.number_of_events)
        return m_n
