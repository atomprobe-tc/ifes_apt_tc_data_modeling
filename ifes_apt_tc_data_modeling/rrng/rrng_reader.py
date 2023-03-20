# RRNG range file reader used by atom probe microscopists.
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


def evaluate_rrng_range_line(i: int, line: str) -> dict:
    """Evaluate information content of a single range line."""
    # example line "Range7 = 42.8160 43.3110 vol:0.04543
    #     Al:1 O:1 Name: AlO Color:00FFFF"
    # according to DOI: 10.1007/978-1-4614-8721-0
    # mqmin, mqmax, vol, ion composition is required,
    #     name and color fields are optional
    info: dict = {}
    info["identifier"] = "Range" + str(i)
    info["range"] = np.asarray([0., MQ_EPSILON], np.float64)
    info["atoms"] = []
    info["volume"] = np.float64(0.)
    info["color"] = ""
    info["name"] = ""

    tmp = re.split(r"[\s=]+", line)
    assert len(tmp) >= 6, "Line " + line + \
        " does not contain all required fields!"
    assert tmp[0] == "Range" + str(i), "Line " + line + \
        " has inconsistent line prefix!"

    assert is_range_significant(np.float64(tmp[1]), np.float64(tmp[2])), \
        "Line " + line + " insignificant range!"
    info["range"] = np.asarray([tmp[1], tmp[2]], np.float64)

    if tmp[3].lower().startswith("vol:"):
        info["volume"] = np.float64(re.split(r":", tmp[3])[1])
    if (tmp[-1].lower().startswith("color:")) and \
        (len(re.split(r":", tmp[-1])[1]) == 6):
        info["color"] = "#" + re.split(r":", tmp[-1])[1]
    # HEX_COLOR_REGEX = r"^([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$"
    # replace r"^#( ...
    # regexp = re.compile(HEX_COLOR_REGEX)
    # if regexp.search(tmp[-1].split(r":")):

    for information in tmp[4:-1]:
        element_multiplicity = re.split(r":+", information)
        assert len(element_multiplicity) == 2, \
           "Element multiplicity is incorrectly formatted!"
        # skip vol, name, and color information
        if element_multiplicity[0].lower() == "name":
            info["name"] = str(element_multiplicity[1])
            # this captures properly formatted ranges with keyword
            # name whose name value is then however not a chemical
            # symbol but some user-defined string
            # this is a common habit to define custom names
        elif element_multiplicity[0].lower() not in ["vol", "color"]:
            # pick up what is an element name
            symbol = element_multiplicity[0]
            assert (symbol in chemical_symbols) and (symbol != "X"), \
                "Line " + line + " contains an invalid chemical symbol!"
            assert np.uint32(element_multiplicity[1]) > 0, \
                "Line " + line + " zero or negative multiplicity!"
            assert np.uint32(element_multiplicity[1]) < 256, \
                "Line " + line + " unsupport high multiplicity!"
            info["atoms"] = np.append(
                info["atoms"], [symbol] * int(element_multiplicity[1]))
    return info


class ReadRrngFileFormat():
    """Read *.rrng file format."""

    def __init__(self, filename: str):
        assert len(filename) > 5, "RRNG file incorrect filename ending!"
        assert filename.lower().endswith(".rrng"), \
            "RRNG file incorrect file type!"
        self.filename = filename

        self.rrng: dict = {}
        self.rrng["ionnames"] = []
        self.rrng["ranges"] = {}
        self.rrng["ions"] = {}
        self.rrng["molecular_ions"] = []
        self.read_rrng()

    def read_rrng(self):
        """Read content of an RRNG range file."""
        with open(self.filename, mode="r", encoding="utf8") as rrngf:
            txt = rrngf.read()

        txt = txt.replace("\r\n", "\n")  # pylint: disable=R0801 # windows to unix EOL conversion
        txt = txt.replace(",", ".")  # pylint: disable=R0801 # use decimal dots instead of comma
        txt_stripped = [line for line in txt.split("\n")  # pylint: disable=R0801
                        if line.strip() != "" and line.startswith("#") is False]  # pylint: disable=R0801
        del txt  # pylint: disable=R0801

        # see DOI: 10.1007/978-1-4899-7430-3 for further details to this  # pylint: disable=R0801
        # AMETEK/Cameca"s *.rrng file format  # pylint: disable=R0801

        # pylint: disable=R0801 # first, parse [Ions] section, which holds a list of element names
        # pylint: disable=R0801 # there are documented cases where experimentalists add custom strings
        # pylint: disable=R0801 # to specify ranges they consider special
        # pylint: disable=R0801 # these are loaded as user types
        # pylint: disable=R0801 # with isotope_vector np.iinfo(np.uint16).max
        where = [idx for idx, element in
                 enumerate(txt_stripped) if element == "[Ions]"]
        assert isinstance(where, list), "Section [Ions] not found!"
        assert len(where) == 1, "Section [Ions] not found or ambiguous!"
        current_line_id = where[0] + 1

        tmp = re.split(r"[\s=]+", txt_stripped[current_line_id])
        assert len(tmp) == 2, "[Ions]/Number line corrupted!"
        assert tmp[0] == "Number", "[Ions]/Number incorrectly formatted!"
        assert tmp[1].isnumeric(), "[Ions]/Number not a number!"
        number_of_ion_names = int(tmp[1])
        assert number_of_ion_names > 0, "No ion names defined!"
        current_line_id += 1
        for i in np.arange(0, number_of_ion_names):
            tmp = re.split(r"[\s=]+", txt_stripped[current_line_id + i])
            assert len(tmp) == 2, "[Ions]/Ion line corrupted!"
            assert tmp[0] == "Ion" + str(i + 1), \
                "[Ions]/Ion incorrectly formatted!"
            assert isinstance(tmp[1], str), "[Ions]/Name not a string!"
            self.rrng["ionnames"].append(tmp[1])  # [tmp[0]] = tmp[1]

        # second, parse [Ranges] section
        where = [idx for idx, element in
                 enumerate(txt_stripped) if element == "[Ranges]"]
        assert isinstance(where, list), "Section [Ranges] not found!"
        assert len(where) == 1, "Section [Ranges] not found or ambiguous!"
        current_line_id = where[0] + 1

        tmp = re.split(r"[\s=]+", txt_stripped[current_line_id])
        assert len(tmp) == 2, "[Ranges]/Number line corrupted!"
        assert tmp[0] == "Number", "[Ranges]/Number incorrectly formatted!"
        assert tmp[1].isnumeric(), "[Ranges]/Number not a number!"
        number_of_ranges = int(tmp[1])
        assert number_of_ranges > 0, "No ranges defined!"
        current_line_id += 1

        # print("Parsing range, progress via recovered charge state...")
        for i in np.arange(0, number_of_ranges):
            dct = evaluate_rrng_range_line(i + 1, txt_stripped[current_line_id + i])
            assert dct, \
                "Line " + txt_stripped[current_line_id + i] + " is corrupted!"

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

            self.rrng["molecular_ions"].append(m_ion)
        print(self.filename + " parsed successfully")

if __name__ == "main":
    pass
    # testing
    # parsedFile = ReadRrngFileFormat("../../R31_06365-v02.rrng")
