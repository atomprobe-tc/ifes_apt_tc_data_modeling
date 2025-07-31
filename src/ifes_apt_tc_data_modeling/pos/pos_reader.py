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

# pylint: disable=duplicate-code

import os
import numpy as np

from ifes_apt_tc_data_modeling.nexus.nx_field import NxField
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data
from ifes_apt_tc_data_modeling.utils.custom_logging import logger


class ReadPosFileFormat:
    """Read *.pos file format."""

    def __init__(self, file_path: str):
        """Initialize the reader."""
        if (len(file_path) <= 4) or (not file_path.lower().endswith(".pos")):
            raise ImportError(
                "WARNING::POS file incorrect file_path ending or file type!"
            )
        self.file_path = file_path

        self.file_size = os.path.getsize(self.file_path)
        if self.file_size % (4 * 4) != 0:
            raise ValueError("POS file_size not integer multiple of 4*4B!")
        if np.uint32(self.file_size / (4 * 4)) >= np.iinfo(np.uint32).max:
            raise ValueError("POS file is too large, currently only 2*32 supported!")
        self.number_of_events = np.uint32(self.file_size / (4 * 4))
        logger.debug(f"Parsing {self.number_of_events} events from {self.file_path}")

        # https://doi.org/10.1007/978-1-4614-3436-8 for file format details
        # dtyp_names = ["Reconstructed position along the x-axis (nm)",
        #               "Reconstructed position along the y-axis (nm)",
        #               "Reconstructed position along the z-axis (nm)",
        #               "Reconstructed mass-to-charge-state ratio (Da)"]

    def get_reconstructed_positions(self):
        """Read xyz columns."""

        xyz = NxField()
        xyz.values = np.zeros([self.number_of_events, 3], np.float32)
        xyz.unit = "nm"

        xyz.values[:, 0] = get_memory_mapped_data(
            self.file_path, ">f4", 0 * 4, 4 * 4, self.number_of_events
        )  # x
        xyz.values[:, 1] = get_memory_mapped_data(
            self.file_path, ">f4", 1 * 4, 4 * 4, self.number_of_events
        )  # y
        xyz.values[:, 2] = get_memory_mapped_data(
            self.file_path, ">f4", 2 * 4, 4 * 4, self.number_of_events
        )  # z
        return xyz

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge-state-ratio column."""

        m_n = NxField()
        m_n.values = np.zeros([self.number_of_events, 1], np.float32)
        m_n.unit = "Da"

        m_n.values[:, 0] = get_memory_mapped_data(
            self.file_path, ">f4", 3 * 4, 4 * 4, self.number_of_events
        )
        return m_n
