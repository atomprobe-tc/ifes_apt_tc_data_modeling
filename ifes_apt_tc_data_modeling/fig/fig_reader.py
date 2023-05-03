# FIG.TXT range file reader used by atom probe microscopists.
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

import re

import numpy as np

from ifes_apt_tc_data_modeling.nexus.nx_ion import NxField, NxIon

from ifes_apt_tc_data_modeling.utils.utils import \
    create_isotope_vector, is_range_significant

from ifes_apt_tc_data_modeling.utils.definitions import \
    MQ_EPSILON, MAX_NUMBER_OF_ATOMS_PER_ION

from ifes_apt_tc_data_modeling.utils.molecular_ions import \
    isotope_to_hash, MolecularIonBuilder, \
    PRACTICAL_ABUNDANCE, PRACTICAL_ABUNDANCE_PRODUCT, \
    PRACTICAL_MIN_HALF_LIFE, VERBOSE, SACRIFICE_ISOTOPIC_UNIQUENESS

from ase.data import atomic_numbers, atomic_masses, chemical_symbols


class ReadFigTxtFileFormat():
    """Read *.fig.txt file format."""

    def __init__(self, filename: str):
        assert len(filename) > 7, "FIG.TXT file incorrect filename ending!"
        assert filename.lower().endswith(".fig.txt"), \
            "FIG.TXT file incorrect file type!"
        self.filename = filename

        self.fig: dict = {}
        self.fig["ranges"] = {}
        self.fig["ions"] = {}
        self.fig["molecular_ions"] = []
        self.read_fig_txt()

    def read_fig_txt(self):
        """Read FIG.TXT range file content."""
        with open(self.filename, mode="r", encoding="utf8") as figf:
            txt = figf.read()

        txt = txt.replace("\r\n", "\n")  # windows to unix EOL conversion
        txt = txt.replace(",", ".")  # use decimal dots instead of comma
        txt_stripped = [line for line in txt.split("\n")
                        if line.strip() != "" and line.startswith("#") is False]
        for molecular_ion in txt_stripped:
            tmp = molecular_ion.split(" ")
            mqmin = np.float64(tmp[len(tmp)-2:-1][0])
            mqmax = np.float64(tmp[len(tmp)-1:][0])
            ionname = " ".join(tmp[:-2])
            # print(f"{ionname} [{mqmin}, {mqmax}]")
            # ionname = '16O 1H2 + + +  + '

            positive = ionname.count("+")
            negative = ionname.count("-")
            if (0 < positive <= 7) and (negative == 0):
                charge_state = positive
            elif (0 < negative <= 7) and (positive == 0):
                charge_state = -negative
            else:
                charge_state = 0

            tmp = ionname.replace("+", "").replace("-", "").split(" ")
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
                    symbol = isotope.replace( \
                        str(mass_number), "").replace( \
                        str(multiplier), "").replace(" ", "")
                    if (symbol != "X") and (symbol in chemical_symbols):
                        for cnt in np.arange(0, multiplier):
                            proton_number = atomic_numbers[symbol]
                            neutron_number = mass_number - proton_number
                            ivec.append(isotope_to_hash(proton_number, neutron_number))
            ivec = np.sort(np.asarray(ivec, np.uint16))[::-1]
            ivector = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,), np.uint16)
            ivector[0:len(ivec)] = ivec
            # print(ivector)

            m_ion = NxIon(isotope_vector=ivector, charge_state=charge_state)
            m_ion.add_range(mqmin, mqmax)
            m_ion.comment = NxField(ionname, "")
            # m_ion.report()

            crawler = MolecularIonBuilder(
                min_abundance=PRACTICAL_ABUNDANCE,
                min_abundance_product=PRACTICAL_ABUNDANCE_PRODUCT,
                min_half_life=PRACTICAL_MIN_HALF_LIFE,
                sacrifice_uniqueness=SACRIFICE_ISOTOPIC_UNIQUENESS,
                verbose=VERBOSE)
            recovered_charge_state, m_ion_candidates = crawler.combinatorics(
                m_ion.isotope_vector.typed_value,
                m_ion.ranges.typed_value[0, 0],
                m_ion.ranges.typed_value[0, 1])
            # print(f"{recovered_charge_state}")
            m_ion.charge_state = NxField(np.int8(recovered_charge_state), "")
            m_ion.update_human_readable_name()
            m_ion.add_charge_state_model(
                {"min_abundance": PRACTICAL_ABUNDANCE,
                 "min_abundance_product": PRACTICAL_ABUNDANCE_PRODUCT,
                 "min_half_life": PRACTICAL_MIN_HALF_LIFE,
                 "sacrifice_isotopic_uniqueness": SACRIFICE_ISOTOPIC_UNIQUENESS},
                m_ion_candidates)

            self.fig["molecular_ions"].append(m_ion)
        print(self.filename + " parsed successfully")

if __name__ == "main":
    pass
    # testing
    # parsedFile = ReadFigTxtFileFormat("../../02763-v01.rng.fig.txt")
