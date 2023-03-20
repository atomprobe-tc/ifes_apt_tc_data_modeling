# RNG range file reader used by atom probe microscopists.
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

# pylint: disable=E1101

import re

import numpy as np

from ifes_apt_tc_data_modeling.nexus.nx_ion import NxField, NxIon
from ifes_apt_tc_data_modeling.utils.utils import \
    create_isotope_vector, is_range_significant
from ifes_apt_tc_data_modeling.utils.definitions import MQ_EPSILON
from ifes_apt_tc_data_modeling.utils.molecular_ions import MolecularIonBuilder
from ifes_apt_tc_data_modeling.utils.molecular_ions import \
    PRACTICAL_ABUNDANCE, PRACTICAL_ABUNDANCE_PRODUCT, \
    PRACTICAL_MIN_HALF_LIFE, VERBOSE, SACRIFICE_ISOTOPIC_UNIQUENESS

from ase.data import atomic_numbers, atomic_masses, chemical_symbols


def evaluate_rng_range_line(
        i: int, line: str, column_id_to_label: dict, n_columns: int) -> dict:
    """Represent information content of a single range line."""
    # example line: ". 107.7240 108.0960 1 0 0 0 0 0 0 0 0 0 3 0 0 0"
    info: dict = {}
    info["identifier"] = "Range" + str(i)
    info["range"] = np.asarray([0., MQ_EPSILON], np.float64)
    info["atoms"] = []
    info["volume"] = np.float64(0.)
    info["color"] = ""
    info["name"] = ""

    tmp = re.split(r"\s+", line)
    assert len(tmp) is n_columns, "Line " + line \
        + " inconsistent number columns!"
    assert tmp[0] == ".", "Line " + line \
        + " has inconsistent line prefix!"

    assert is_range_significant(np.float64(tmp[1]), np.float64(tmp[2])), \
        "Line " + line + " insignificant range!"
    info["range"] = np.asarray([tmp[1], tmp[2]], np.float64)

    # line encodes multiplicity of element via array of multiplicity counts
    element_multiplicity = np.asarray(tmp[3:len(tmp)], np.uint32)
    assert np.sum(element_multiplicity) > 0, \
        "Line " + line + " no element counts!"
    if np.sum(element_multiplicity) > 0:
        for j in np.arange(0, len(element_multiplicity)):
            assert element_multiplicity[j] >= 0, "Line " + line \
               + " no negative element counts!"
            if element_multiplicity[j] > 0:
                symbol = column_id_to_label[j + 1]
                if (symbol in chemical_symbols) and (symbol != "X"):
                    info["atoms"] = np.append(
                        info["atoms"], [column_id_to_label[j + 1]] \
                        * int(element_multiplicity[j]))
                else:
                    info["name"] = symbol
                    info["atoms"] = []  # will map to unknown type

    # color for RNG files can only be decoded by
    # loading the color of elements and polyatomic extensions
    # and then check to which category (element or polyatomic) a range
    # belongs and then take this color, there is almost no way
    # to make something so simple for disentangled
    return info


def evaluate_rng_ion_type_header(line: str) -> dict:
    """Represent information content in the key header line."""
    # line = "------------------- Fe Mg Al Mn Si V C Ga Ti Ca O Na Co H"
    # line = "---- a"
    # line = "----------------- Sc Fe O C Al Si Cr H unknown"
    info = {}
    info["column_id_to_label"] = {}
    tmp = re.split(r"\s+", line)
    assert len(tmp) > 1, "RNG file does not contain iontype labels!"
    for i in np.arange(1, len(tmp)):
        info["column_id_to_label"][i] = tmp[i]
    return info


class ReadRngFileFormat():
    """Read *.rng file format."""

    def __init__(self, filename: str):
        assert len(filename) > 4, "RNG file incorrect filename ending!"
        assert filename.lower().endswith(".rng"), \
            "RNG file incorrect file type!"
        self.filename = filename

        self.rng: dict = {}
        self.rng["ranges"] = {}
        self.rng["ions"] = {}
        self.rng["molecular_ions"] = []
        self.read_rng()

    def read_rng(self):
        """Read RNG range file content."""
        with open(self.filename, mode="r", encoding="utf8") as rngf:
            txt = rngf.read()

        txt = txt.replace("\r\n", "\n")  # pylint: disable=R0801 # windows to unix EOL conversion
        txt = txt.replace(",", ".")  # pylint: disable=R0801 # use decimal dots instead of comma
        txt_stripped = [line for line in txt.split("\n")  # pylint: disable=R0801
                        if line.strip() != "" and line.startswith("#") is False]  # pylint: disable=R0801
        del txt  # pylint: disable=R0801

        # see DOI: 10.1007/978-1-4899-7430-3 for further details to this
        # Oak Ridge National Lab / Oxford *.rng file format
        # only the first ------ line is relevant
        # it details all ion labels aka ions
        # AMETEK"s IVAS/APSuite-specific trailing
        # polyatomic extension is redundant info

        tmp = None
        current_line_id = int(0)  # search key header line
        for line in txt_stripped:
            tmp = re.search(r"----", line)
            if tmp is None:
                current_line_id += int(1)
            else:
                break
        assert tmp is not None, "RNG file does not contain key header line!"

        header = evaluate_rng_ion_type_header(txt_stripped[current_line_id])

        tmp = re.split(r"\s+", txt_stripped[0])
        assert tmp[0].isnumeric() is True, "Number of species corrupted!"
        n_element_symbols = int(tmp[0])
        assert n_element_symbols >= 0, "No species defined!"
        assert tmp[1].isnumeric() is True, "Number of ranges corrupted!"
        n_ranges = int(tmp[1])
        assert n_ranges >= 0, "No ranges defined!"

        for i in np.arange(current_line_id + 1,
                           current_line_id + 1 + n_ranges):
            dct = evaluate_rng_range_line(
                i - current_line_id, txt_stripped[i],
                header["column_id_to_label"],
                n_element_symbols + 3)
            assert dct is not None, \
                "Line " + txt_stripped[i] + " is corrupted!"

            m_ion = NxIon(isotope_vector=create_isotope_vector(
                dct["atoms"]), charge=0)
            m_ion.add_range(dct["range"][0], dct["range"][1])
            m_ion.comment = NxField(dct["name"], "")
            m_ion.color = NxField(dct["color"], "")
            m_ion.volume = NxField(dct["volume"], "")
            # m_ion.report()

            crawler = MolecularIonBuilder(
                min_abundance=PRACTICAL_ABUNDANCE,
                min_abundance_product=PRACTICAL_ABUNDANCE_PRODUCT,
                min_half_life=PRACTICAL_MIN_HALF_LIFE,
                sacrifice_uniqueness=SACRIFICE_ISOTOPIC_UNIQUENESS,
                verbose=VERBOSE)
            recovered_charge, m_ion_candidates = crawler.combinatorics(
                m_ion.isotope_vector.typed_value,
                m_ion.ranges.typed_value[0, 0],
                m_ion.ranges.typed_value[0, 1])
            # print(" " + str(recovered_charge))
            m_ion.charge = NxField(np.int8(recovered_charge), "")
            m_ion.update_human_readable_name()
            m_ion.add_charge_model(
                {"min_abundance": PRACTICAL_ABUNDANCE,
                 "min_abundance_product": PRACTICAL_ABUNDANCE_PRODUCT,
                 "min_half_life": PRACTICAL_MIN_HALF_LIFE,
                 "sacrifice_isotopic_uniqueness": SACRIFICE_ISOTOPIC_UNIQUENESS},
                m_ion_candidates)

            self.rng["molecular_ions"].append(m_ion)
        print(self.filename + " parsed successfully")

if __name__ == "main":
    pass
    # testing
    # parsedFile = ReadRngFileFormat("../../SeHoKim_R5076_44076_v02.rng")
