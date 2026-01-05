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

"""Reader for ranging defs extracted from FAU/Erlangen Atom Probe Toolbox Matlab figures FIG.TXT."""

# pylint: disable=too-many-locals

import re

import numpy as np
from ase.data import atomic_numbers

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.definitions import (
    MAX_NUMBER_OF_ATOMS_PER_ION,
    NEUTRON_NUMBER_FOR_ELEMENT,
)
from ifes_apt_tc_data_modeling.utils.molecular_ions import (
    get_chemical_symbols,
    isotope_to_hash,
)
from ifes_apt_tc_data_modeling.utils.nx_ion import NxIon


class ReadFigTxtFileFormat:
    """Read *.fig.txt file format."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 8) or not file_path.lower().endswith(".fig.txt"):
            raise ImportError(
                "WARNING::FIG.TXT file incorrect file_path ending or file type."
            )
        self.file_path = file_path
        self.fig: dict = {"ranges": {}, "ions": {}, "molecular_ions": []}
        self.read_fig_txt()

    def read_fig_txt(self):
        """Read FIG.TXT range file content."""
        with open(self.file_path, encoding="utf8") as fig_fp:
            txt = fig_fp.read()

        txt = txt.replace("\r\n", "\n")  # windows to unix EOL conversion
        txt = txt.replace(",", ".")  # use decimal dots instead of comma
        txt_stripped = [
            line
            for line in txt.split("\n")
            if line.strip() != "" and line.startswith("#") is False
        ]
        for molecular_ion in txt_stripped:
            tmp = molecular_ion.split(" ")
            mqmin = np.float64(tmp[len(tmp) - 2 : -1][0])
            mqmax = np.float64(tmp[len(tmp) - 1 :][0])
            ion_name = " ".join(tmp[:-2])
            # logger.debug(f"{ion_name} [{mqmin}, {mqmax}]")
            # ion_name = '16O 1H2 + + +  + '

            positive = ion_name.count("+")
            negative = ion_name.count("-")
            if (0 < positive <= 7) and (negative == 0):
                charge_state = positive
            elif (0 < negative <= 7) and (positive == 0):
                charge_state = -negative
            else:
                charge_state = 0

            tmp = ion_name.replace("+", "").replace("-", "").split(" ")
            ivec = []
            for isotope in tmp:
                if isotope != "":
                    prefix = re.findall("^[0-9]+", isotope)
                    mass_number = 0
                    if len(prefix) == 1:
                        if int(prefix[0]) > 0:
                            mass_number = int(prefix[0])
                    suffix = re.findall("[0-9]+$", isotope)
                    multiplier = 1
                    if len(suffix) == 1:
                        multiplier = int(suffix[0])
                    symbol = (
                        isotope.replace(f"{mass_number}", "")
                        .replace(f"{multiplier}", "")
                        .replace(" ", "")
                    )
                    if symbol in get_chemical_symbols():
                        proton_number = atomic_numbers[symbol]
                        neutron_number = NEUTRON_NUMBER_FOR_ELEMENT
                        if mass_number != 0:
                            neutron_number = mass_number - proton_number
                        ivec.extend(
                            [isotope_to_hash(proton_number, neutron_number)]
                            * multiplier
                        )
            ivec = np.sort(np.asarray(ivec, np.uint16))[::-1]
            ivector = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,), np.uint16)
            ivector[0 : len(ivec)] = ivec

            m_ion = NxIon(nuclide_hash=ivector, charge_state=charge_state)
            m_ion.add_range(mqmin, mqmax)
            m_ion.comment = ion_name
            m_ion.apply_combinatorics()
            # m_ion.report()

            self.fig["molecular_ions"].append(m_ion)
        logger.info(f"{self.file_path} parsed successfully.")
