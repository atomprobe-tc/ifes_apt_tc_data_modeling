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

from ifes_apt_tc_data_modeling.nexus.nx_field import NxField
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data


class ReadEposFileFormat():
    """Read *.epos file format."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 5) or (file_path.lower().endswith(".epos") is False):
            raise ImportError("WARNING::ePOS file incorrect file_path ending or file type!")
        self.file_path = file_path
        self.file_size = os.path.getsize(self.file_path)
        assert self.file_size % 11 * 4 == 0, \
            "ePOS file_size not integer multiple of 11*4B!"
        assert np.uint32(self.file_size / (11 * 4)) < np.iinfo(np.uint32).max, \
            "ePOS file is too large, currently only 2*32 supported!"
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

        xyz = NxField()
        xyz.values = np.zeros([self.number_of_events, 3], np.float32)
        xyz.unit = "nm"

        xyz.values[:, 0] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   0 * 4, 11 * 4, self.number_of_events)  # x
        xyz.values[:, 1] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   1 * 4, 11 * 4, self.number_of_events)  # y
        xyz.values[:, 2] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   2 * 4, 11 * 4, self.number_of_events)  # z
        return xyz

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge-state-ratio column."""
        m_n = NxField()
        m_n.values = np.zeros([self.number_of_events, 1], np.float32)
        m_n.unit = "Da"

        m_n.values[:, 0] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   3 * 4, 11 * 4, self.number_of_events)
        return m_n

    def get_raw_time_of_flight(self):
        """Read raw (uncorrected) time-of-flight."""
        raw_tof = NxField()
        raw_tof.values = np.zeros([self.number_of_events, 1], np.float32)
        raw_tof.unit = "ns"

        # according to DOI: 10.1007/978-1-4899-7430-3 raw time-of-flight
        # i.e. this is an uncorrected time-of-flight
        # for which effects uncorrect?
        # Only the proprietary IVAS/APSuite source code knows for sure
        raw_tof.values[:, 0] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   4 * 4, 11 * 4, self.number_of_events)
        return raw_tof

    def get_standing_voltage(self):
        """Read standing voltage."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # standing voltage on the specimen
        # according to DOI: 10.1007/978-1-4614-8721-0 also-known as DC voltage
        dc_voltage = NxField()
        dc_voltage.values = np.zeros([self.number_of_events, 1], np.float32)
        dc_voltage.unit = "kV"
        # different to the above-mentioned references Gault et al. state
        # that standing and pulse_voltage are in V instead of kV

        dc_voltage.values[:, 0] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   5 * 4, 11 * 4, self.number_of_events)
        return dc_voltage

    def get_pulse_voltage(self):
        """Read pulse voltage."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # additional voltage to trigger field evaporation in case
        # of high-voltage pulsing, 0 for laser pulsing
        pu_voltage = NxField()
        pu_voltage.values = np.zeros([self.number_of_events, 1], np.float32)
        pu_voltage.unit = "kV"

        pu_voltage.values[:, 0] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   6 * 4, 11 * 4, self.number_of_events)
        return pu_voltage

    def get_hit_positions(self):
        """Read ion impact positions on detector."""
        hit_positions = NxField()
        hit_positions.values = np.zeros([self.number_of_events, 2], np.float32)
        hit_positions.unit = "mm"

        hit_positions.values[:, 0] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   7 * 4, 11 * 4, self.number_of_events)  # x
        hit_positions.values[:, 1] = \
            get_memory_mapped_data(self.file_path, ">f4",
                                   8 * 4, 11 * 4, self.number_of_events)  # y
        return hit_positions

    def get_number_of_pulses(self):
        """Read number of pulses."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # number of pulses since last event detected
        # 0 after the first ion per pulse
        # also known as $\Delta Pulse$
        npulses = NxField()
        npulses.values = np.zeros([self.number_of_events, 1], np.uint32)
        npulses.unit = ""

        npulses.values[:, 0] = \
            get_memory_mapped_data(self.file_path, ">u4",
                                   9 * 4, 11 * 4, self.number_of_events)
        return npulses

    def get_ions_per_pulse(self):
        """Read ions per pulse."""
        # according to DOI: 10.1007/978-1-4899-7430-3
        # ions per pulse, 0 after the first ion
        ions_per_pulse = NxField()
        ions_per_pulse.values = np.zeros([self.number_of_events, 1], np.uint32)
        ions_per_pulse.unit = ""

        ions_per_pulse.values[:, 0] = \
            get_memory_mapped_data(self.file_path, ">u4",
                                   10 * 4, 11 * 4, self.number_of_events)
        return ions_per_pulse
