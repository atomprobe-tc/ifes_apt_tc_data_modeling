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

"""Reader for Cameca HDF5 files as shared here https://doi.org/10.18126/dqxb-9m77."""

# pylint: disable=fixme

import os
import h5py
import re
import numpy as np

from ifes_apt_tc_data_modeling.utils.nx_ion import NxIon
from ifes_apt_tc_data_modeling.cameca.cameca_utils import parse_elements
from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg
from ifes_apt_tc_data_modeling.utils.utils import (
    create_nuclide_hash,
    is_range_significant,
)
from ifes_apt_tc_data_modeling.utils.custom_logging import logger


class ReadCamecaHfiveFileFormat:
    """Read Cameca HDF5 files used in this study https://doi.org/10.18126/dqxb-9m77."""

    def __init__(self, file_path: str):
        if not file_path.lower().endswith((".hdf5")):
            raise ImportError(
                "WARNING::HDF5 file incorrect file_path ending or file type."
            )
        self.file_path = file_path
        self.file_size = os.path.getsize(self.file_path)
        self.number_of_events = None
        self.supported = 0  # voting-based
        with h5py.File(self.file_path, "r") as h5r:
            for entry in ("mass", "xyz", "ranges"):
                if entry in h5r.keys():
                    self.supported += 1
            self.number_of_events = np.shape(h5r["mass"][...])[0]
            if (
                self.number_of_events == 0
                or self.number_of_events != np.shape(h5r["xyz"][...])[0]
            ):
                logger.warning(
                    f"{self.file_path} length of mass and xyz arrays is zero or differs."
                )
                self.supported = 0
                return
        if self.supported == 3:
            logger.debug(f"{self.file_path} is a supported Cameca HDF5 file.")
        else:
            logger.warning(f"{self.file_path} is not a supported Cameca HDF5 file.")

    def get_reconstructed_positions(self):
        """Read xyz positions."""
        values = np.zeros((self.number_of_events,), np.float32)
        with h5py.File(self.file_path, "r") as h5r:
            values[:, :] = np.asarray(h5r["mass"][...], np.float32)
        return ureg.Quantity(values, ureg.nanometer)

    def get_mass_to_charge_state_ratio(self):
        """Read (calibrated?) mass-to-charge-state-ratio."""
        values = np.zeros((self.number_of_events,), np.float32)
        with h5py.File(self.file_path, "r") as h5r:
            values[:] = np.asarray(h5r["mass"][...], np.float32)
        return ureg.Quantity(values, ureg.dalton)

    def get_ranging_definitions(self):
        """Read ranging definitions."""
        self.rng: dict = {"ranges": {}, "ions": {}, "molecular_ions": []}

        with h5py.File(self.file_path, "r") as h5r:
            for entry in h5r["ranges"]:
                # clarify if reports following the creation order from Cameca or differently?
                logger.info(entry)
                match = re.compile(r"^range_\d+$", entry)
                if not match or not all(
                    attribute_name in h5r["ranges"][entry].attrs
                    for attribute_name in ["element", "min_da", "max_da"]
                ):  # , "volume"]:
                    continue
                # ion_id = int(entry.replace("range_", ""))
                atoms = parse_elements(f"""{h5r["ranges"][entry].attrs["element"]}""")
                min_mass_to_charge = np.float64(
                    h5r["ranges"][entry].attrs["min_da"][()]
                )  # in examples was already f64
                max_mass_to_charge = np.float64(
                    h5r["ranges"][entry].attrs["max_da"][()]
                )  # in examples was already f64
                if not is_range_significant(min_mass_to_charge, max_mass_to_charge):
                    logger.warning(f"Range for {entry} is not significant.")
                    continue
                m_ion = NxIon(nuclide_hash=create_nuclide_hash(atoms), charge_state=0)
                m_ion.comment = (
                    entry  # for debugging above-mentioned issue if order retained
                )
                m_ion.add_range(min_mass_to_charge, max_mass_to_charge)
                m_ion.apply_combinatorics()

                self.rng["molecular_ions"].append(m_ion)
        logger.info(f"{self.file_path} ranging definitions parsed successfully.")
