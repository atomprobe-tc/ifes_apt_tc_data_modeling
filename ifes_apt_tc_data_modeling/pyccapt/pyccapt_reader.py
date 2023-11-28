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

import h5py

import numpy as np

import pandas as pd

from ase.data import atomic_numbers, chemical_symbols
from ifes_apt_tc_data_modeling.nexus.nx_ion import NxIon
from ifes_apt_tc_data_modeling.nexus.nx_field import NxField
from ifes_apt_tc_data_modeling.utils.utils import \
    isotope_to_hash, isotope_vector_to_nuclid_list, MAX_NUMBER_OF_ATOMS_PER_ION

# this implementation focuses on the following state of the pyccapt repository
# https://github.com/mmonajem/pyccapt/commit/e955beb4f2627befb8b4d26f2e74e4c52e00394e

# during the course of an atom probe measurement and analysis with FAU/Erlangen's Oxcart instrument
# several HDF5 files are generated with essentially two software tools. One is pyccapt which has a
# a control module, a calibration module (where the voltage/bowl calibration and reconstruction is performed),
# and a module/functionalities to document ranging i.e. ion type identification made
# The other software typically used by the FAU/Erlangen atom probe group is Atom Probe Toolbox;
# instructed as a set of Matlab live scripts this toolbox offers data analysis functionalities,
# results are stored via an HDF5 file

# specific comments 
# pyccapt/control
# an HDF5 file keeping relevant quantities 

# pyccapt/calibration
# unfortunately the generated HDF5 file has internally no provenance information
# with which pyccapt version it was generated, therefore developers of pyccapt should
# rather write the content of the HDF5 file explicitly dset by dset e.g. using h5py instead
# of the pandas HDF5 dump convenience functionality
# of course pandas stores its own version but that is not conclusive enough to infer with
# which pyccapt version and most importantly from which other context the file was generated
# this is an aspect of the FAIR RDM principles which the pyccapt approach currently ignores


class ReadPyccaptControlFileFormat():
    """Read FAU/Erlangen pyccapt (controle module) HDF5 file format."""

    def __init__(self, filename: str):
        assert len(filename) > 2, "H5 file incorrect filename ending!"
        assert filename.lower().endswith(".h5") or filename.lower().endswith(".hdf5"), \
            "HDF5 file incorrect file type!"
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)
        self.number_of_events = None
        self.version = "e955beb4f2627befb8b4d26f2e74e4c52e00394e"

        # check that the formatting matches that of an pyccapt control module output HDF5 file
        with h5py.File(self.filename, "r") as h5r:
            self.supported = 0  # voting-based
            required_groups = ["apt", "dld", "tdc"]
            for req_grpnm in required_groups:
                if req_grpnm in h5r.keys():
                    self.supported += 1
            if self.supported == 3:
                print(f"{self.filename} is a supported pyccapt/control HDF5 file!")
            else:
                print(f"{self.filename} is not a supported pyccapt/control HDF5 file!")
                return


class ReadPyccaptCalibrationFileFormat():
    """Read FAU/Erlangen pyccapt (calibration module) HDF5 file format."""

    def __init__(self, filename: str):
        assert len(filename) > 2, "H5 file incorrect filename ending!"
        assert filename.lower().endswith(".h5") or filename.lower().endswith(".hdf5"), \
            "HDF5 file incorrect file type!"
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)
        self.number_of_events = None
        self.version = "e955beb4f2627befb8b4d26f2e74e4c52e00394e"
        self.df = None

        with h5py.File(self.filename, "r") as h5r:
            self.supported = 0  # voting-based
            required_entries = ["df",
                                "df/axis0", "df/axis1",
                                "df/block0_items", "df/block0_values",
                                "df/block1_items", "df/block1_values"]
            for entry in required_entries:
                if entry in h5r.keys():
                    self.supported += 1
            if self.supported == 7:
                print(f"{self.filename} is a supported pyccapt/calibration HDF5 file!")
            else:
                print(f"{self.filename} is not a supported pyccapt/calibration HDF5 file!")
                return

        self.df = pd.read_hdf(self.filename)
        self.number_of_events = np.shape(self.df)[0]
    
    def get_named_quantities(self, term: str):
        if term in self.df.keys():
            return self.df[term]
        return None
    
    def get_reconstructed_positions(self):
        """Read xyz columns."""

        xyz = NxField()
        xyz.typed_value = np.zeros(
            [self.number_of_events, 3], np.float32)
        xyz.unit = "nm"

        dim = 0
        for quant in ["x (nm)", "y (nm)", "z (nm)"]:
            xyz.typed_value[:, dim] = np.asarray(self.get_named_quantities(quant), np.float32)
            dim += 1
        return xyz

    def get_mass_to_charge_state_ratio(self):
        """Read (calibrated) mass-to-charge-state-ratio column."""

        m_n = NxField()
        m_n.typed_value = np.zeros(
            [self.number_of_events, 1], np.float32)
        m_n.unit = "Da"

        m_n.typed_value[:, 0] = np.asarray(self.get_named_quantities("mc_c (Da)"), np.float32)
        return m_n


class ReadPyccaptRangingFileFormat():
    """Read FAU/Erlangen pyccapt (ranging module) HDF5 file format."""

    def __init__(self, filename: str):
        assert len(filename) > 2, "H5 file incorrect filename ending!"
        assert filename.lower().endswith(".h5") or filename.lower().endswith(".hdf5"), \
            "HDF5 file incorrect file type!"
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)
        self.number_of_events = None
        self.version = "e955beb4f2627befb8b4d26f2e74e4c52e00394e"
        self.df = None

        with h5py.File(self.filename, "r") as h5r:
            self.supported = 0  # voting-based
            required_entries = ["df",
                                "df/axis0", "df/axis1",
                                "df/block0_items", "df/block0_values",
                                "df/block1_items", "df/block1_values",
                                "df/block2_items", "df/block2_values"]
            for entry in required_entries:
                if entry in h5r.keys():
                    self.supported += 1
            if self.supported == 9:
                print(f"{self.filename} is a supported pyccapt/ranging HDF5 file!")
            else:
                print(f"{self.filename} is not a supported pyccapt/ranging HDF5 file!")
                return

        self.df = pd.read_hdf(self.filename)
        self.rng = {}
        self.rng["molecular_ions"] = []
        print(np.shape(self.df)[0])
        for idx in np.arange(0, np.shape(self.df)[0]):
            if isinstance(self.df.iloc[idx, 6], str) is True:
                if self.df.iloc[idx, 6] == "unranged":
                    continue

            elements = self.df.iloc[idx, 6]
            complexs = self.df.iloc[idx, 7]
            isotopes = self.df.iloc[idx, 8]
            # assertions
            ivec = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,), np.uint16)
            hashvector = []
            for idxj in np.arange(0, len(elements)):
                symbol = elements[idxj]
                if symbol in chemical_symbols and symbol != "X":
                    proton_number = atomic_numbers[symbol]
                    neutron_number = isotopes[idxj] - proton_number
                    for mult in np.arange(0, complexs[idxj]):
                        hashvector.append(isotope_to_hash(proton_number, neutron_number))
            ivec[0:len(hashvector)] = np.sort(np.asarray(hashvector, np.uint16), kind="stable")[::-1]
            
            m_ion = NxIon()
            m_ion.isotope_vector.typed_value = ivec
            m_ion.nuclid_list.typed_value = isotope_vector_to_nuclid_list(ivec)
            m_ion.charge_state.typed_value = np.int8(self.df.iloc[idx, 9])
            m_ion.add_range(self.df.iloc[idx, 3], self.df.iloc[idx, 4])
            m_ion.update_human_readable_name()
            # m_ion.report()
            self.rng["molecular_ions"].append(m_ion)
        print(f"{self.filename} parsed successfully")
