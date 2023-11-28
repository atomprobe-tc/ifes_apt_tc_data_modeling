# POS file format reader used by atom probe microscopists.
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

# pylint: disable=no-member,duplicate-code

import os

import numpy as np

from ifes_apt_tc_data_modeling.nexus.nx_field import NxField
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data


class ReadAtoFileFormat():
    """Read Rouen group *.ato file format."""

    def __init__(self, filename: str):
        assert len(filename) > 4, "ATO file incorrect filename ending!"
        assert filename.lower().endswith(".ato"), \
            "ATO file incorrect file type!"
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)
        self.number_of_events = None
        self.version = None
        retval = self.get_ato_version()
        if retval in [3, 5]:
            # there also seems to exist a version 4 but I have never seen an example for it
            self.version = retval
            print(f"ATO file is in a supported version {self.version}")
            if self.version == 3:
                assert (self.filesize - 2 * 4) % 14 * 4 == 0, \
                    "ATO v3 filesize not integer multiple of 14*4B!"
                self.number_of_events = np.uint32((self.filesize - 2 * 4) / (14 * 4))
                print(f"ATO file contains {self.number_of_events} entries")
            if self.version == 5:
                assert (self.filesize - 5000) % 40 == 0, \
                    "ATO v5 filesize not integer multiple of 40B!"
                self.number_of_events = np.uint32((self.filesize - 5000) / 40)
                print(f"ATO file contains {self.number_of_events} entries")                
        else:
            raise ValueError("ATO file unsupported version!")
        # https://zenodo.org/records/8382828
        # details three versions of the Rouen/GPM ato format v3, v4, v5
        # Cameca/AMETEK's runrootl/FileConvert utility know two ATO flavours:
        # CamecaRoot v18.46.533g built Marc, 21, 2022 against ROOT 5.34/36
        # v3 LAWATOP and v5 current GPM
        # specifically an earlier parser
        # https://hg.sr.ht/~mycae/libatomprobe/browse/src/io/dataFiles.cpp?rev=tip
        # mentions that storage format may not be robust enough against overflow and
        # suggests that additional polishing of results is needed

    def get_ato_version(self):
        header = get_memory_mapped_data(self.filename, "<u4", 0, 4, 2)
        if header[1] in [3, 4, 5]:
            return header[1]
        return None

    def get_reconstructed_positions(self):
        """Read xyz columns."""

        xyz = NxField()
        xyz.typed_value = np.zeros(
            [self.number_of_events, 3], np.float32)
        xyz.unit = "nm"

        if self.version == 3:
            for dim in [0, 1, 2]:
                xyz.typed_value[:, dim] = \
                    get_memory_mapped_data(self.filename, "<f4",
                        2 * 4 + dim * 4, 14 * 4, self.number_of_events)
                # wpx -> x, wpy -> y, fpz -> z
        if self.version == 5:
            # publicly available sources are inconclusive whether coordinates are in angstroem or nm
            # based on the evidence of usa_denton_smith Si.epos converted to v5 ATO via CamecaRoot
            # the resulting x, y coordinates suggests that v5 ATO stores in angstroem, while fpz is stored in nm?
            # however https://zenodo.org/records/8382828 reports the reconstructed positions to be named
            # not at all wpx, wpy and fpz but x, y, z instead and here claims the nm
            xyz.typed_value[:, 0] = \
                np.float32(get_memory_mapped_data(self.filename, "<i2",
                           5000 + 0, 40, self.number_of_events) * 0.1)  # wpx -> x
            xyz.typed_value[:, 1] = \
                np.float32(get_memory_mapped_data(self.filename, "<i2",
                           5000 + 2, 40, self.number_of_events) * 0.1)  # wpy -> y
            xyz.typed_value[:, 2] = \
                get_memory_mapped_data(self.filename, "<f4",
                                       5000 + 4, 40, self.number_of_events)  # fpz -> z
        return xyz

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge-state-ratio column."""

        m_n = NxField()
        m_n.typed_value = np.zeros(
            [self.number_of_events, 1], np.float32)
        m_n.unit = "Da"

        if self.version == 3:
            m_n.typed_value[:, 0] = \
                get_memory_mapped_data(self.filename, "<f4",
                                    2 * 4 + 3 * 4, 14 * 4, self.number_of_events)
        if self.version == 5:
            m_n.typed_value[:, 0] = \
                get_memory_mapped_data(self.filename, "<f4",
                                    5000 + 8, 40, self.number_of_events)
        return m_n
