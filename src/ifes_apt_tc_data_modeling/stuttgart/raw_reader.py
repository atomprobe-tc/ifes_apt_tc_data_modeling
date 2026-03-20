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

"""Stuttgart RAW file format reader used by atom probe microscopists."""

# pylint: disable=duplicate-code

import os

import numpy as np

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.default_config import MAXIMUM_NUMBER_OF_IONS
from ifes_apt_tc_data_modeling.utils.memory_mapped_io import get_memory_mapped_data
from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg


class ReadStuttgartApytRawFileFormat:
    """Read *.raw file format."""

    def __init__(self, file_path: str, verbose: bool = False):
        self.supported = False
        if not file_path.lower().endswith(".raw"):
            logger.warning(f"{file_path} is likely not a Stuttgart TAP-style raw file")
            return
        self.file_path = file_path
        self.verbose = verbose
        self.file_size = os.path.getsize(self.file_path)
        if self.file_size % (8 * 4) != 0:
            logger.warning(
                f"{self.file_path} APyT raw file_size not integer multiple of 8 * 4B."
            )
            return
        if np.uint32(self.file_size / (8 * 4)) > MAXIMUM_NUMBER_OF_IONS:
            logger.warning(
                f"{self.file_path} APyT raw file with more than {MAXIMUM_NUMBER_OF_IONS} events is not supported."
            )
            return
        self.number_of_events = int(self.file_size / (8 * 4))
        self.supported = True

        # https://github.com/sebi-85/apyt/blob/424b33ea5fee2d7e8b3fc0a9af68f382b8da341f/apyt/io/conv.py#L26-33 for file format details
        # Field         Data type  Description
        # ============  =========  =================================
        # U_base        float32    base voltage (V)
        # U_pulse       float32    pulse voltage (V)
        # U_reflectron  float32    reflectron voltage (V)
        # x_det         float32    `x` detector position (mm)
        # y_det         float32    `y` detector position (mm)
        # tof           float32    time of flight (ns)
        # epoch         int32      epoch of evaporation event
        # pulse_num     uint32     pulse number of evaporation event
        # dtype_names = []
        # raw = np.fromfile( fnm, dtype= {"names": dtype_names,
        # "formats": (, "<f4","<f4","<f4","<f4","<f4","<f4","<i4","<u4") } )

    def get_base_voltage(self) -> ureg.Quantity | None:
        """Read (uncorrected) base voltage column."""
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = "<f4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                0 * item_size,
                (8 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.volt)
            else:
                logger.warning("Unable to get_base_voltage")
        return None

    def get_pulse_voltage(self) -> ureg.Quantity | None:
        """Read (uncorrected) pulse voltage column."""
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = "<f4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                1 * item_size,
                (8 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.volt)
            else:
                logger.warning("Unable to get_pulse_voltage")
        return None

    def get_reflectron_voltage(self) -> ureg.Quantity | None:
        """Read (uncorrected) reflectron voltage column."""
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = "<f4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                2 * item_size,
                (8 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.volt)
            else:
                logger.warning("Unable to get_reflectron_voltage")
        return None

    def get_raw_detector_position(self) -> ureg.Quantity | None:
        """Read (uncorrected) detector position columns."""
        if self.supported:
            values = np.zeros((self.number_of_events, 2), np.float32)
            type_literal = "<f4"
            item_size = np.dtype(type_literal).itemsize
            all_values = True
            for column_index in [0, 1]:  # x, y
                data = get_memory_mapped_data(
                    self.file_path,
                    type_literal,
                    (3 + column_index) * item_size,
                    (8 * item_size,),
                    (self.number_of_events,),
                )
                if data is not None:
                    np.copyto(values[:, column_index], data, casting="unsafe")
                else:
                    all_values = False
                    logger.warning(
                        f"Unable to get_raw_detector_position dim {column_index}"
                    )
            if all_values:
                return ureg.Quantity(values, ureg.millimeter)
            else:
                logger.warning("Unable to get_raw_detector_position")
        return None

    def get_raw_time_of_flight(self) -> ureg.Quantity | None:
        """Read (uncorrected) time-of-flight column."""
        if self.supported:
            values = np.zeros((self.number_of_events,), np.float32)
            type_literal = "<f4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                5 * item_size,
                (8 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values, ureg.nanosecond)
            else:
                logger.warning("Unable to get_raw_time_of_flight")
        return None

    def get_epoch_evaporation_event(self) -> ureg.Quantity | None:
        """Read epoch of evaporation event column."""
        if self.supported:
            values = np.zeros((self.number_of_events,), np.int32)
            type_literal = "<i4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                (6 * item_size),
                (8 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values)
            else:
                logger.warning("Unable to get_epoch_evaporation_event")
        return None

    def get_pulse_number_evaporation_event(self) -> ureg.Quantity | None:
        """Read epoch of evaporation event column."""
        if self.supported:
            values = np.zeros((self.number_of_events,), np.uint32)
            type_literal = "<u4"
            item_size = np.dtype(type_literal).itemsize
            data = get_memory_mapped_data(
                self.file_path,
                type_literal,
                (7 * item_size),
                (8 * item_size,),
                (self.number_of_events,),
            )
            if data is not None:
                np.copyto(values[:], data, casting="unsafe")
                return ureg.Quantity(values)
            else:
                logger.warning("Unable to get_pulse_number_evaporation_event")
        return None
