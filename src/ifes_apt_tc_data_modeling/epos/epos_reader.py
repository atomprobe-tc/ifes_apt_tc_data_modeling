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

"""ePOS file format reader used by atom probe microscopists."""

# pylint: disable=duplicate-code

import os

import numpy as np

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.default_config import MAXIMUM_NUMBER_OF_IONS
from ifes_apt_tc_data_modeling.utils.memory_mapped_io import get_memory_mapped_data
from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg


class ReadEposFileFormat:
    """Read *.epos file format."""

    def __init__(self, file_path: str, verbose: bool = False):
        self.supported = False
        if not file_path.lower().endswith(".epos"):
            logger.warning(f"{file_path} is likely not an ePOS file")
            return
        self.file_path = file_path
        self.verbose = verbose
        self.file_size = os.path.getsize(self.file_path)

        # https://doi.org/10.1007/978-1-4614-3436-8 for file format details
        # dtype_names = ["Reconstructed position along the x-axis (nm)",
        #               "Reconstructed position along the y-axis (nm)",
        #               "Reconstructed position along the z-axis (nm)",
        #               "Reconstructed mass-to-charge-state ratio (Da)",
        #               "Raw time-of-flight (ns)",
        #               "Standing voltage (V)",
        #               "Pulsed voltage (V)",
        #               "Ion impact x-coordinate at the detector (mm)",
        #               "Ion impact y-coordinate at the detector (mm)",
        #               "Number of pulses since the last detected ion (pulses)",
        #               "Hit multiplicity (ions)"]
        # raw = np.fromfile( fnm, dtype= {"names": dtype_names,
        # "formats": (, ">f4",">f4",">f4",">f4",">f4",">f4",">u4",">u4") } )

        if self.file_size % (11 * 4) != 0:
            logger.warning(
                f"{self.file_path} file_size does not match expectation for ePOS file."
            )
            return
        if np.uint32(self.file_size / (11 * 4)) > MAXIMUM_NUMBER_OF_IONS:
            logger.warning(
                f"{self.file_path} with more than {MAXIMUM_NUMBER_OF_IONS} events is not supported."
            )
            return
        self.number_of_events = int(self.file_size / (11 * 4))
        logger.info(
            f"{self.file_path} is supported, holds {self.number_of_events} events"
        )
        self.supported = True

    def get_reconstructed_positions(self) -> ureg.Quantity | None:
        """Read xyz columns."""
        if self.supported:
            values = np.zeros((self.number_of_events, 3), np.float32)
            type_literal = ">f4"
            item_size = np.dtype(
                type_literal
            ).itemsize  # >f4 and >u4 have the same itemsize
            all_values = True
            for column_index in [0, 1, 2]:  # x, y, z
                data = get_memory_mapped_data(
                    self.file_path,
                    type_literal,
                    column_index * item_size,
                    (11 * item_size,),
                    (self.number_of_events,),
                )
                if data is not None:
                    np.copyto(values[:, column_index], data, casting="unsafe")
                else:
                    all_values = False
                    logger.warning(
                        "Unable to get_reconstructed positions for column_index {column_index}"
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
            item_size = np.dtype(
                type_literal
            ).itemsize  # >f4 and >u4 have the same itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                3 * item_size,
                (11 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.dalton)
            else:
                logger.warning("Unable to get_mass_to_charge_state_ratio")
        return None

    def get_raw_time_of_flight(self) -> ureg.Quantity | None:
        """Read raw (uncorrected) time-of-flight."""
        # according to DOI: 10.1007/978-1-4899-7430-3 raw time-of-flight
        # i.e. this is an uncorrected time-of-flight
        # for which effects incorrect?
        # Only the proprietary IVAS/APSuite source code knows for sure
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = ">f4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                4 * item_size,
                (11 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.nanosecond)
            else:
                logger.warning("Unable to get_raw_time_of_flight")
        return None

    def get_standing_voltage(self) -> ureg.Quantity | None:
        """Read standing voltage."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # standing voltage on the specimen
        # according to DOI: 10.1007/978-1-4614-8721-0 also-known as DC voltage
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = ">f4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                5 * item_size,
                (11 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.kilovolt).to(ureg.volt)
            else:
                logger.warning("Unable to get_standing_voltage")
        return None

    def get_pulse_voltage(self) -> ureg.Quantity | None:
        """Read pulse voltage."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # additional voltage to trigger field evaporation in case
        # of high-voltage pulsing, 0 for laser pulsing
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = ">f4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                6 * item_size,
                (11 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.kilovolt).to(ureg.volt)
            else:
                logger.warning("Unable to get_pulse_voltage")
        return None

    def get_hit_positions(self) -> ureg.Quantity | None:
        """Read ion impact positions on detector."""
        if self.supported:
            values = np.zeros((self.number_of_events, 2), np.float32)
            type_literal = ">f4"
            item_size = np.dtype(type_literal).itemsize
            all_values = True
            for column_index in [0, 1]:  # x, y
                data = get_memory_mapped_data(
                    self.file_path,
                    type_literal,
                    (7 + column_index) * 4,
                    (11 * item_size,),
                    (self.number_of_events,),
                )
                if data is not None:
                    np.copyto(values[:, column_index], data, casting="unsafe")
                else:
                    all_values = False
                    logger.warning("Unable to get_hit_positions dim {dim}")
            if all_values:
                return ureg.Quantity(values, ureg.millimeter)
            else:
                logger.warning("Unable to get_hit_positions")
        return None

    def get_number_of_pulses(self) -> ureg.Quantity | None:
        """Read number of pulses."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # number of pulses since last event detected
        # 0 after the first ion per pulse
        # also known as $\Delta Pulse$
        if self.supported:
            values = np.zeros((self.number_of_events,), np.uint32)
            type_literal = ">u4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                ">u4",
                9 * item_size,
                (11 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values)
            else:
                logger.warning("Unable to get_number_of_pulses")
        return None

    def get_ions_per_pulse(self) -> ureg.Quantity | None:
        """Read ions per pulse."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # ions per pulse, 0 after the first ion
        if self.supported:
            values = np.zeros((self.number_of_events,), np.uint32)
            type_literal = ">u4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                10 * item_size,
                (11 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values)
            else:
                logger.warning("Unable to get_ions_per_pulse")
        return None
