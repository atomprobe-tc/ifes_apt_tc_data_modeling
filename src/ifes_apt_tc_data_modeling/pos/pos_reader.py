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

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.default_config import MAXIMUM_NUMBER_OF_IONS
from ifes_apt_tc_data_modeling.utils.memory_mapped_io import get_memory_mapped_data
from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg


class ReadPosFileFormat:
    """Read *.pos file format."""

    def __init__(self, file_path: str, verbose: bool = False):
        """Initialize the reader."""
        self.supported = False
        if not file_path.lower().endswith(".pos"):
            logger.warning(f"{file_path} is likely not a POS file")
            return
        self.file_path = file_path
        self.verbose = verbose
        self.file_size = os.path.getsize(self.file_path)

        # https://doi.org/10.1007/978-1-4614-3436-8 for file format details
        # dtype_names = ["Reconstructed position along the x-axis (nm)",
        #               "Reconstructed position along the y-axis (nm)",
        #               "Reconstructed position along the z-axis (nm)",
        #               "Reconstructed mass-to-charge-state ratio (Da)"]
        # big-endian floats >f4
        item_size = np.dtype(">f4").itemsize
        if self.file_size % (item_size * item_size) != 0:
            logger.warning(
                f"{self.file_path} file_size does not match expectation for POS file."
            )
            return
        if np.uint32(self.file_size / (item_size * item_size)) > MAXIMUM_NUMBER_OF_IONS:
            logger.warning(
                f"{self.file_path} with more than {MAXIMUM_NUMBER_OF_IONS} events is not supported."
            )
            return
        self.number_of_events = int(self.file_size / (item_size * item_size))
        logger.info(
            f"{self.file_path} is supported, holds {self.number_of_events} events"
        )
        self.supported = True

    def get_reconstructed_positions(self) -> ureg.Quantity | None:
        """Read xyz columns."""
        if self.supported:
            values = np.zeros((self.number_of_events, 3), np.float32)
            type_literal = ">f4"
            item_size = np.dtype(type_literal).itemsize
            all_values = True
            for column_index in [0, 1, 2]:  # x, y, z
                data = get_memory_mapped_data(
                    self.file_path,
                    type_literal,
                    column_index * item_size,
                    (4 * item_size,),
                    (self.number_of_events,),
                )
                if data is not None:
                    np.copyto(values[:, column_index], data, casting="unsafe")
                else:
                    all_values = False
                    logger.warning(
                        f"Unable to get_reconstructed_positions for column_index {column_index}"
                    )
            if all_values:
                return ureg.Quantity(values, ureg.nanometer)
            else:
                logger.warning("Unable to get_reconstructed_positions")
        return None

    def get_mass_to_charge_state_ratio(self) -> ureg.Quantity | None:
        """Read mass-to-charge-state-ratio column."""
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = ">f4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                3 * item_size,
                (4 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.dalton)
            else:
                logger.warning("Unable to get_mass_to_charge_state_ratio")
        return None
