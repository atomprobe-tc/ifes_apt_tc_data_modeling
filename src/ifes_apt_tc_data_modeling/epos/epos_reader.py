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

from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data


class ReadEposFileFormat:
    """Read *.epos file format."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 5) or not file_path.lower().endswith(".epos"):
            raise ImportError(
                "WARNING::ePOS file incorrect file_path ending or file type."
            )
        self.file_path = file_path
        self.file_size = os.path.getsize(self.file_path)
        if self.file_size % (11 * 4) != 0:
            raise ValueError("ePOS file_size not integer multiple of 11*4B.")
        if np.uint32(self.file_size / (11 * 4)) >= np.iinfo(np.uint32).max:
            raise ValueError("ePOS file is too large, currently only 2*32 supported.")
        self.number_of_events = np.uint32(self.file_size / (11 * 4))

        # https://doi.org/10.1007/978-1-4614-3436-8 for file format details
        # dtyp_names = ["Reconstructed position along the x-axis (nm)",
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
        # raw = np.fromfile( fnm, dtype= {"names": dtyp_names,
        # "formats": (, ">f4",">f4",">f4",">f4",">f4",">f4",">u4",">u4") } )

    def get_reconstructed_positions(self):
        """Read xyz columns."""
        values = np.zeros((self.number_of_events, 3), np.float32)
        for dim in [0, 1, 2]:  # x, y, z
            values[:, dim] = get_memory_mapped_data(
                self.file_path, ">f4", dim * 4, 11 * 4, self.number_of_events
            )
        return ureg.Quantity(values, ureg.nanometer)

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge-state-ratio column."""
        values = np.asarray((self.number_of_events,), np.float32)
        values[:] = get_memory_mapped_data(
            self.file_path, ">f4", 3 * 4, 11 * 4, self.number_of_events
        )
        return ureg.Quantity(values, ureg.dalton)

    def get_raw_time_of_flight(self):
        """Read raw (uncorrected) time-of-flight."""
        values = np.zeros((self.number_of_events,), np.float32)
        # according to DOI: 10.1007/978-1-4899-7430-3 raw time-of-flight
        # i.e. this is an uncorrected time-of-flight
        # for which effects uncorrect?
        # Only the proprietary IVAS/APSuite source code knows for sure
        values[:] = get_memory_mapped_data(
            self.file_path, ">f4", 4 * 4, 11 * 4, self.number_of_events
        )
        return ureg.Quantity(values, ureg.nanosecond)

    def get_standing_voltage(self):
        """Read standing voltage."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # standing voltage on the specimen
        # according to DOI: 10.1007/978-1-4614-8721-0 also-known as DC voltage
        values = np.zeros((self.number_of_events,), np.float32)
        values[:] = get_memory_mapped_data(
            self.file_path, ">f4", 5 * 4, 11 * 4, self.number_of_events
        )
        return ureg.Quantity(values, ureg.kilovolt)

    def get_pulse_voltage(self):
        """Read pulse voltage."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # additional voltage to trigger field evaporation in case
        # of high-voltage pulsing, 0 for laser pulsing
        values = np.zeros((self.number_of_events,), np.float32)
        values[:] = get_memory_mapped_data(
            self.file_path, ">f4", 6 * 4, 11 * 4, self.number_of_events
        )
        return ureg.Quantity(values, ureg.kilovolt)

    def get_hit_positions(self):
        """Read ion impact positions on detector."""
        values = np.zeros((self.number_of_events, 2), np.float32)
        for dim in [0, 1]:  # x, y
            values[:, dim] = get_memory_mapped_data(
                self.file_path, ">f4", (7 + dim) * 4, 11 * 4, self.number_of_events
            )
        return ureg.Quantity(values, ureg.millimeter)

    def get_number_of_pulses(self):
        """Read number of pulses."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # number of pulses since last event detected
        # 0 after the first ion per pulse
        # also known as $\Delta Pulse$
        values = np.zeros((self.number_of_events,), np.uint32)
        values[:] = get_memory_mapped_data(
            self.file_path, ">u4", 9 * 4, 11 * 4, self.number_of_events
        )
        return ureg.Quantity(values)

    def get_ions_per_pulse(self):
        """Read ions per pulse."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # ions per pulse, 0 after the first ion
        values = np.zeros((self.number_of_events,), np.uint32)
        values[:] = get_memory_mapped_data(
            self.file_path, ">u4", 10 * 4, 11 * 4, self.number_of_events
        )
        return ureg.Quantity(values)
