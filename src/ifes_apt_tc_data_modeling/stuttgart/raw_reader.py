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
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data
from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg


class ReadStuttgartApytRawFileFormat:
    """Read *.raw file format."""

    def __init__(self, file_path: str, verbose: bool = False):
        self.supported = False
        if not file_path.lower().endswith(".raw"):
            logger.warning(f"{file_path} is likely not a Stuttgart TAP-style raw file")
            return
        self.supported = True
        self.verbose = verbose
        self.file_path = file_path
        self.file_size = os.path.getsize(self.file_path)
        if self.file_size % (8 * 4) != 0:
            raise ValueError("Stuttgart raw file_size not integer multiple of 8*4B.")
        if np.uint32(self.file_size / (8 * 4)) >= np.iinfo(np.uint32).max:
            raise ValueError(
                "Stuttgart raw file is too large, currently only 2*32 supported."
            )
        self.number_of_events = np.uint32(self.file_size / (8 * 4))

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

    def get_base_voltage(self):
        """Read (uncorrected) base voltage column."""
        values = np.zeros((self.number_of_events,), np.float32)
        data = get_memory_mapped_data(
            self.file_path, "<f4", 0 * 4, 8 * 4, self.number_of_events
        )
        if data is not None:
            np.copyto(values[:], data, casting="unsafe")
            return ureg.Quantity(values, ureg.volt)
        else:
            logger.warning("Unable to get_base_voltage")

    def get_pulse_voltage(self):
        """Read (uncorrected) pulse voltage column."""
        values = np.zeros((self.number_of_events,), np.float32)
        data = get_memory_mapped_data(
            self.file_path, "<f4", 1 * 4, 8 * 4, self.number_of_events
        )
        if data is not None:
            np.copyto(values[:], data, casting="unsafe")
            return ureg.Quantity(values, ureg.volt)
        else:
            logger.warning("Unable to get_pulse_voltage")

    def get_reflectron_voltage(self):
        """Read (uncorrected) reflectron voltage column."""
        values = np.zeros((self.number_of_events,), np.float32)
        data = get_memory_mapped_data(
            self.file_path, "<f4", 2 * 4, 8 * 4, self.number_of_events
        )
        if data is not None:
            np.copyto(values[:], data, casting="unsafe")
            return ureg.Quantity(values, ureg.volt)
        else:
            logger.warning("Unable to get_reflectron_voltage")

    def get_raw_detector_position(self):
        """Read (uncorrected) detector position columns."""
        values = np.zeros((self.number_of_events, 2), np.float32)
        all_values = True
        for dim in [0, 1]:  # x, y
            data = get_memory_mapped_data(
                self.file_path, "<f4", (3 + dim) * 4, 8 * 4, self.number_of_events
            )
            if data is not None:
                np.copyto(values[:, dim], data, casting="unsafe")
            else:
                all_values = False
                logger.warning(f"Unable to get_raw_detector_position dim {dim}")
        if all_values:
            return ureg.Quantity(values, ureg.millimeter)

    def get_raw_time_of_flight(self):
        """Read (uncorrected) time-of-flight column."""
        values = np.zeros((self.number_of_events,), np.float32)
        data = get_memory_mapped_data(
            self.file_path, "<f4", 5 * 4, 8 * 4, self.number_of_events
        )
        if data is not None:
            np.copyto(values[:], data, casting="unsafe")
            return ureg.Quantity(values, ureg.nanosecond)
        else:
            logger.warning("Unable to get_raw_time_of_flight")

    def get_epoch_evaporation_event(self):
        """Read epoch of evaporation event column."""
        values = np.zeros((self.number_of_events,), np.int32)
        data = get_memory_mapped_data(
            self.file_path, "<i4", 6 * 4, 8 * 4, self.number_of_events
        )
        if data is not None:
            np.copyto(values[:], data, casting="unsafe")
            return ureg.Quantity(values)
        else:
            logger.warning("Unable to get_epoch_evaporation_event")

    def get_pulse_number_evaporation_event(self):
        """Read epoch of evaporation event column."""
        values = np.zeros((self.number_of_events,), np.uint32)
        data = get_memory_mapped_data(
            self.file_path, "<u4", 7 * 4, 8 * 4, self.number_of_events
        )
        if data is not None:
            np.copyto(values[:], data, casting="unsafe")
            return ureg.Quantity(values)
        else:
            logger.warning("Unable to get_pulse_number_evaporation_event")
