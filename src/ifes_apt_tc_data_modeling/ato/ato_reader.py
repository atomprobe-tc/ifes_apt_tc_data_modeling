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

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data
from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg


class ReadAtoFileFormat:
    """Read Rouen group *.ato file format."""

    def __init__(self, file_path: str, verbose: bool = False):
        self.supported = False
        if not file_path.lower().endswith(".ato"):
            logger.warning(
                f"{file_path} is likely not an ATO file as those from atom probe"
            )
            return
        self.supported = True
        self.file_path = file_path
        self.verbose = verbose

        self.file_size = os.path.getsize(self.file_path)
        self.number_of_events = None
        self.version = None
        retval = self.get_ato_version()
        if retval in [3, 5]:
            # there also seems to exist a version 4 but I have never seen an example for it
            self.version = retval
            logger.info(f"ATO file is in a supported version {self.version}")
            if self.version == 3:
                if (self.file_size - 2 * 4) % (14 * 4) != 0:
                    raise ValueError("ATO v3 file_size not integer multiple of 14*4B.")

                self.number_of_events = np.uint32((self.file_size - 2 * 4) / (14 * 4))
                logger.info(f"ATO file contains {self.number_of_events} entries.")
            if self.version == 5:
                if (self.file_size - 5000) % 40 != 0:
                    raise ValueError("ATO v5 file_size not integer multiple of 40B.")
                self.number_of_events = np.uint32((self.file_size - 5000) / 40)
                logger.info(f"ATO file contains {self.number_of_events} entries.")
        else:
            raise ImportError("ATO file unsupported version.")
        # https://zenodo.org/records/8382828
        # details three versions of the Rouen/GPM ato format v3, v4, v5
        # Cameca/AMETEK's runrootl/FileConvert utility know two ATO flavours:
        # CamecaRoot v18.46.533g built Marc, 21, 2022 against ROOT 5.34/36
        # v3 LAWATAP and v5 current GPM
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

    def get_reconstructed_positions(self):
        """Read xyz columns."""

        values = np.zeros((self.number_of_events, 3), np.float32)
        # with "<f4" we assume little-endian here, yes it is in contradiction
        # to the statement made in https://link.springer.com/content/pdf/bbm:978-1-4614-8721-0/1?pdf=chapter%20toc
        # analyzed examples are all consistent with all reported evidence that ATO v3 files
        # have a two times four byte header followed by records with 14 * 4 B each
        if self.version == 3:
            # wpx -> x, wpy -> y, fpz -> z
            all_values = True
            for dim in [0, 1, 2]:
                data = get_memory_mapped_data(
                    self.file_path,
                    "<f4",
                    2 * 4 + dim * 4,
                    14 * 4,
                    self.number_of_events,
                )
                if data is not None:
                    np.multiply(data, 0.1, out=values[:, dim], casting="unsafe")
                else:
                    all_values = False
                    logger.warning("Unable to get_reconstructed positions dim {dim}")
            if all_values:
                return ureg.Quantity(values, ureg.nanometer)
        elif self.version == 5:
            # publicly available sources are inconclusive whether coordinates are in angstrom or nm
            # based on the evidence of usa_denton_smith Si.epos converted to v5 ATO via CamecaRoot
            # the resulting x, y coordinates suggests that v5 ATO stores in angstrom, while fpz is stored in nm?
            # however https://zenodo.org/records/8382828 reports the reconstructed positions to be named
            # not at all wpx, wpy and fpz but x, y, z instead and here claims the nm
            # wpx -> x, wpy -> y
            all_values = True
            for dim in [0, 1]:
                data = get_memory_mapped_data(
                    self.file_path,
                    "<i2",
                    5000 + (dim * 2),
                    40,
                    self.number_of_events,
                )
                if data is not None:
                    np.multiply(data, 0.01, out=values[:, dim], casting="unsafe")
                else:
                    all_values = False
                    logger.warning("Unable to get_reconstructed positions dim {dim}")
            # angstrom to nm conversion for wpx and wpy was dropped to make results consistent with
            # APSuite based file format conversion tool, again a signature that the ATO format
            # demands better documentation by those who use it especially if claiming to perform
            # FAIR research, nothing about the documentation of this format is currently ticking the
            # FAIR principles but rather software development is prohibited because of contradictory/insufficient
            # documentation
            # fpz -> z
            data = get_memory_mapped_data(
                self.file_path, "<f4", 5000 + 4, 40, self.number_of_events
            )
            if data is not None:
                np.multiply(data, 0.1, out=values[:, 2], casting="unsafe")
            else:
                all_values = False
                logger.warning("Unable to get_reconstructed positions dim 2")
        if all_values:
            return ureg.Quantity(values, ureg.nanometer)

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge-state-ratio column."""

        values = np.zeros((self.number_of_events,), np.float32)
        # see comment under respective function for xyz,
        # problem is m/q values typically are in the lower part of the byte
        # because of which some examples yield even physical reasonable
        # m/q values when reading with dtype_v3 = ">f4" !
        # what if the individual machines were different types of endianness?
        # this is another significant problem with just sharing files without
        # any self-documentation or context around it
        if self.version == 3:
            data = get_memory_mapped_data(
                self.file_path, "<f4", 2 * 4 + 3 * 4, 14 * 4, self.number_of_events
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.dalton)
        elif self.version == 5:
            data = get_memory_mapped_data(
                self.file_path, "<f4", 5000 + 8, 40, self.number_of_events
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.dalton)
        logger.warning("Unable to get_mass_to_charge_state_ratio")
