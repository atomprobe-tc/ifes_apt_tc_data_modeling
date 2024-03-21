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

"""ATO file format reader used by atom probe microscopists."""

import os
import numpy as np

from ifes_apt_tc_data_modeling.nexus.nx_field import NxField
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data


class ReadAtoFileFormat():
    """Read Rouen group *.ato file format."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 4) or (file_path.lower().endswith(".ato") is False):
            raise ImportError("WARNING::ATO file incorrect file_path ending or file type!")
        self.file_path = file_path

        self.file_size = os.path.getsize(self.file_path)
        self.number_of_events = None
        self.version = None
        retval = self.get_ato_version()
        if retval in [3, 5]:
            # there also seems to exist a version 4 but I have never seen an example for it
            self.version = retval
            print(f"ATO file is in a supported version {self.version}")
            if self.version == 3:
                assert (self.file_size - 2 * 4) % 14 * 4 == 0, \
                    "ATO v3 file_size not integer multiple of 14*4B!"
                self.number_of_events = np.uint32((self.file_size - 2 * 4) / (14 * 4))
                print(f"ATO file contains {self.number_of_events} entries")
            if self.version == 5:
                assert (self.file_size - 5000) % 40 == 0, \
                    "ATO v5 file_size not integer multiple of 40B!"
                self.number_of_events = np.uint32((self.file_size - 5000) / 40)
                print(f"ATO file contains {self.number_of_events} entries")
        else:
            raise ImportError("ATO file unsupported version!")
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
        """Identify if file_path matches a known ATO format version."""
        header = get_memory_mapped_data(self.file_path, "<u4", 0, 4, 2)
        # one can use little-endian <u4 as i8 and u8 for value 3 are degenerated
        # for little and big endian
        if header[1] in [3, 4, 5]:
            return header[1]
        return None

    def get_reconstructed_positions(self):
        """Read xyz columns."""

        xyz = NxField()
        xyz.values = np.zeros([self.number_of_events, 3], np.float32)
        xyz.unit = "nm"
        dtype_v3 = "<f4"  # we assume little-endian here, yes it is in contradiction
        # to the statement made in https://link.springer.com/content/pdf/bbm:978-1-4614-8721-0/1?pdf=chapter%20toc
        # analyzed examples are all consistent with all reported evidence that ATO v3 files
        # have a two times four byte header followed by records with 14 * 4 B each

        if self.version == 3:
            for dim in [0, 1, 2]:
                xyz.values[:, dim] = \
                    np.float32(get_memory_mapped_data(self.file_path, dtype_v3,
                                                      2 * 4 + dim * 4,
                                                      14 * 4, self.number_of_events) * 0.1)
                # wpx -> x, wpy -> y, fpz -> z
        if self.version == 5:
            # publicly available sources are inconclusive whether coordinates are in angstroem or nm
            # based on the evidence of usa_denton_smith Si.epos converted to v5 ATO via CamecaRoot
            # the resulting x, y coordinates suggests that v5 ATO stores in angstroem, while fpz is stored in nm?
            # however https://zenodo.org/records/8382828 reports the reconstructed positions to be named
            # not at all wpx, wpy and fpz but x, y, z instead and here claims the nm
            xyz.values[:, 0] = \
                np.float32(get_memory_mapped_data(self.file_path, "<i2",
                                                  5000 + 0,
                                                  40, self.number_of_events)) * 0.01  # wpx -> x
            xyz.values[:, 1] = \
                np.float32(get_memory_mapped_data(self.file_path, "<i2",
                                                  5000 + 2,
                                                  40, self.number_of_events)) * 0.01  # wpy -> y
            # angstroem to nm conversion for wpx and wpy was dropped to make results consistent with
            # APSuite based file format conversion tool, again a signature that the ATO format
            # demands better documentation by those who use it especially if claiming to perform
            # FAIR research, nothing about the documentation of this format is currently ticking the
            # FAIR principles but rather software development is prohibited because of contradictory/insufficient
            # documentation
            xyz.values[:, 2] = \
                np.float32(get_memory_mapped_data(self.file_path, "<f4",
                                                  5000 + 4, 40, self.number_of_events) * 0.1)  # fpz -> z
        return xyz

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge-state-ratio column."""

        m_n = NxField()
        m_n.values = np.zeros([self.number_of_events, 1], np.float32)
        m_n.unit = "Da"
        dtype_v3 = "<f4"  # see comment under respective function for xyz,
        # problem is m/q values typically are in the lower part of the byte
        # because of which some examples yield even physical reasonable
        # m/q values when reading with dtype_v3 = ">f4" !
        # what if the individual machines were different types of endianness?
        # this is another significant problem with just sharing files without
        # any self-documentation or context around it

        if self.version == 3:
            m_n.values[:, 0] = \
                get_memory_mapped_data(self.file_path, dtype_v3,
                                       2 * 4 + 3 * 4, 14 * 4, self.number_of_events)
        if self.version == 5:
            m_n.values[:, 0] = \
                get_memory_mapped_data(self.file_path, "<f4",
                                       5000 + 8, 40, self.number_of_events)
        return m_n
