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
from ifes_apt_tc_data_modeling.utils.default_config import MAXIMUM_NUMBER_OF_IONS
from ifes_apt_tc_data_modeling.utils.memory_mapped_io import get_memory_mapped_data
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
        self.file_path = file_path
        self.verbose = verbose

        self.file_size = os.path.getsize(self.file_path)
        self.number_of_events = None
        self.version = None

        header = get_memory_mapped_data(self.file_path, "<u4", 0, (4,), (2,))
        # one can use little-endian <u4 as i8 and u8 for value 3 le and be degenerated
        if header[1] in [3, 5]:
            # there also seems to exist a version 4 but I have never seen an example for it
            # pyccapt source code suggests there is also a version 6, also here no example for it
            self.version = header[1]
            logger.info(
                f"{self.file_path} ATO file is a supported version {self.version}"
            )
            if self.version == 3:
                if (self.file_size - 2 * 4) % (14 * 4) != 0:
                    logger.warning(
                        f"{self.file_path} ATO v3 file_size is not an integer multiple of 14*4B."
                    )
                    return
                if (
                    np.uint32((self.file_size - 2 * 4) / (14 * 4))
                    > MAXIMUM_NUMBER_OF_IONS
                ):
                    logger.warning(
                        f"{self.file_path} with more than {MAXIMUM_NUMBER_OF_IONS} events is not supported."
                    )
                    return
                self.number_of_events = int((self.file_size - 2 * 4) / (14 * 4))
                logger.info(
                    f"{self.file_path} is supported, ATO v3, holds {self.number_of_events} events."
                )
                self.supported = True
            elif self.version == 5:
                if (self.file_size - 5000) % 40 != 0:
                    logger.warning(
                        f"{self.file_path} ATO v5 file_size is not an integer multiple of 40B."
                    )
                    return
                if np.uint32((self.file_size - 5000) / 40) > MAXIMUM_NUMBER_OF_IONS:
                    logger.warning(
                        f"{self.file_path} with more than {MAXIMUM_NUMBER_OF_IONS} events is not supported."
                    )
                    return
                self.number_of_events = int((self.file_size - 5000) / 40)
                logger.info(
                    f"{self.file_path} is supported, ATO v5, holds {self.number_of_events} events."
                )
                self.supported = True
        # https://zenodo.org/records/8382828
        # details three versions of the Rouen/GPM ato format v3, v4, v5
        # Cameca/AMETEK's runrootl/FileConvert utility know two ATO flavours:
        # CamecaRoot v18.46.533g built Marc, 21, 2022 against ROOT 5.34/36
        # v3 LAWATAP and v5 current GPM
        # specifically an earlier parser
        # https://hg.sr.ht/~mycae/libatomprobe/browse/src/io/dataFiles.cpp?rev=tip
        # mentions that storage format may not be robust enough against overflow and
        # suggests that additional polishing of results is needed

    def get_reconstructed_positions(self) -> ureg.Quantity | None:
        """Read xyz columns."""
        if self.supported:
            values = np.zeros((self.number_of_events, 3), np.float32)
            # with "<f4" we assume little-endian here, yes it is in contradiction
            # to the statement made in https://link.springer.com/content/pdf/bbm:978-1-4614-8721-0/1?pdf=chapter%20toc
            # analyzed examples are all consistent with all reported evidence that ATO v3 files
            # have a two times four byte header followed by records with 14 * 4 B each
            if self.version == 3:
                # wpx -> x, wpy -> y, fpz -> z
                type_literal = "<f4"
                item_size = np.dtype(type_literal).itemsize
                all_values = True
                for column_index in [0, 1, 2]:
                    data = get_memory_mapped_data(
                        self.file_path,
                        type_literal,
                        2 * 4 + column_index * item_size,
                        (14 * item_size,),
                        (self.number_of_events,),
                    )
                    if data is not None:
                        np.multiply(
                            data, 0.1, out=values[:, column_index], casting="unsafe"
                        )
                    else:
                        all_values = False
                        logger.warning(
                            "Unable to get_reconstructed positions dim {dim}"
                        )
                if all_values:
                    return ureg.Quantity(values, ureg.nanometer)
                else:
                    logger.warning("Unable to get_reconstructed_positions")
            elif self.version == 5:
                # publicly available sources are inconclusive whether coordinates are in angstrom or nm
                # based on the evidence of usa_denton_smith Si.epos converted to v5 ATO via CamecaRoot
                # the resulting x, y coordinates suggests that v5 ATO stores in angstrom, while fpz is stored in nm?
                # however https://zenodo.org/records/8382828 reports the reconstructed positions to be named
                # not at all wpx, wpy and fpz but x, y, z instead and here claims the nm
                # wpx -> x, wpy -> y
                item_size = np.dtype("<i2").itemsize
                all_values = True
                for column_index in [0, 1]:
                    data = get_memory_mapped_data(
                        self.file_path,
                        "<i2",
                        5000 + (column_index * 2),
                        (40,),
                        (self.number_of_events,),
                    )
                    if data is not None:
                        np.multiply(
                            data, 0.01, out=values[:, column_index], casting="unsafe"
                        )
                    else:
                        all_values = False
                        logger.warning(
                            "Unable to get_reconstructed positions column_index {column_index}"
                        )
                # angstrom to nm conversion for wpx and wpy was dropped to make results consistent with
                # APSuite based file format conversion tool, again a signature that the ATO format
                # demands better documentation by those who use and defined  it especially if claiming to perform
                # FAIR research, nothing about the documentation of this format is currently ticking the
                # FAIR principles but rather software development is prohibited because of contradictory/insufficient
                # documentation
                # fpz -> z
                data = get_memory_mapped_data(
                    self.file_path, "<f4", 5000 + 4, (40,), (self.number_of_events,)
                )
                if data is not None:
                    np.multiply(data, 0.1, out=values[:, 2], casting="unsafe")
                else:
                    all_values = False
                    logger.warning("Unable to get_reconstructed positions dim 2")
                if all_values:
                    return ureg.Quantity(values, ureg.nanometer)
                else:
                    logger.warning("Unable to get_reconstructed_positions")
        return None

    def get_mass_to_charge_state_ratio(self) -> ureg.Quantity | None:
        """Read mass-to-charge-state-ratio column."""
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = "<f4"
            item_size = np.dtype(type_literal).itemsize
            # see comment under respective function for xyz,
            # problem is m/q values typically are in the lower part of the byte
            # because of which some examples yield even physical reasonable
            # m/q values when reading with dtype_v3 = ">f4" !
            # what if the individual machines were different types of endianness?
            # this is another significant problem with just sharing files without
            # any self-documentation or context around it
            if self.version == 3:
                data = get_memory_mapped_data(
                    self.file_path,
                    type_literal,
                    2 * 4 + 3 * item_size,
                    (14 * item_size,),
                    (self.number_of_events,),
                )
                if data is not None:
                    np.copyto(values[:], data, casting="unsafe")
                    return ureg.Quantity(values, ureg.dalton)
                else:
                    logger.warning("Unable to get_mass_to_charge_state_ratio")
            elif self.version == 5:
                data = get_memory_mapped_data(
                    self.file_path,
                    type_literal,
                    5000 + 8,
                    (40,),
                    (self.number_of_events,),
                )
                if data is not None:
                    np.copyto(values[:], data, casting="unsafe")
                    return ureg.Quantity(values, ureg.dalton)
                else:
                    logger.warning("Unable to get_mass_to_charge_state_ratio")
        return None
