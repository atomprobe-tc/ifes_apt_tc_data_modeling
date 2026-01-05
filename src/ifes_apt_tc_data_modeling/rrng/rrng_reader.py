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

"""RRNG range file reader used by atom probe microscopists."""

# pylint: disable=too-many-branches,too-many-statements,duplicate-code

import re

import numpy as np

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.definitions import MQ_EPSILON
from ifes_apt_tc_data_modeling.utils.molecular_ions import get_chemical_symbols
from ifes_apt_tc_data_modeling.utils.nx_ion import (
    NxIon,
    try_to_reduce_to_unique_definitions,
)
from ifes_apt_tc_data_modeling.utils.utils import (
    create_nuclide_hash,
    is_range_significant,
)


def evaluate_rrng_range_line(i: int, line: str) -> dict:
    """Evaluate information content of a single range line."""
    # example line "Range7 = 42.8160 43.3110 vol:0.04543
    #     Al:1 O:1 Name: AlO Color:00FFFF"
    # according to DOI: 10.1007/978-1-4614-8721-0
    # mqmin, mqmax, vol, ion composition is required,
    #     name and color fields are optional
    info: dict = {
        "identifier": f"Range{i}",
        "range": np.asarray([0.0, MQ_EPSILON], np.float64),
        "atoms": [],
        "volume": np.float64(0.0),
        "color": "",
        "name": "",
    }

    tmp = re.split(r"[\s=]+", line)
    if len(tmp) < 6:
        # raise ValueError(f"Line {line} does not contain all required fields {len(tmp)}.")
        return None
    if not is_range_significant(np.float64(tmp[1]), np.float64(tmp[2])):
        # raise ValueError(f"Line {line} insignificant range.")
        return None
    info["range"] = np.asarray([tmp[1], tmp[2]], np.float64)

    if tmp[3].lower().startswith("vol:"):
        info["volume"] = np.float64(re.split(r":", tmp[3])[1])
    if (tmp[-1].lower().startswith("color:")) and (
        len(re.split(r":", tmp[-1])[1]) == 6
    ):
        info["color"] = "#" + re.split(r":", tmp[-1])[1]
    # HEX_COLOR_REGEX = r"^([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$"
    # replace r"^#( ...
    # regexp = re.compile(HEX_COLOR_REGEX)
    # if regexp.search(tmp[-1].split(r":")):

    for information in tmp[4:-1]:
        element_multiplicity = re.split(r":+", information)
        if len(element_multiplicity) != 2:
            raise ValueError(
                f"Line {line}, element multiplicity is not "
                f"correctly formatted {len(element_multiplicity)}."
            )
        # skip vol, name, and color information
        if element_multiplicity[0].lower() == "name":
            info["name"] = f"{element_multiplicity[1]}"
            # this captures properly formatted ranges with keyword
            # name whose name value is then however not a chemical
            # symbol but some user-defined string
            # this is a common habit to define custom names
        elif element_multiplicity[0].lower() not in ["vol", "color"]:
            # pick up what is an element name
            symbol = element_multiplicity[0]
            if symbol not in get_chemical_symbols():
                # raise ValueError(f"WARNING::Line {line} contains an invalid chemical symbol {symbol}.")
                return info
            # if np.uint32(element_multiplicity[1]) <= 0:
            # raise ValueError(f"Line {line} zero or negative multiplicity .")
            if np.uint32(element_multiplicity[1]) >= 256:
                # raise ValueError(f"Line {line} unsupported high multiplicity "
                #                  f"{np.uint32(element_multiplicity)}.")
                return info
            info["atoms"] = np.append(
                info["atoms"], [symbol] * int(element_multiplicity[1])
            )
    return info


class ReadRrngFileFormat:
    """Read *.rrng file format."""

    def __init__(self, file_path: str, unique=False, verbose=False):
        """Initialize the class."""
        if (len(file_path) <= 5) or not file_path.lower().endswith(".rrng"):
            raise ImportError(
                "WARNING::RRNG file incorrect file_path ending or file type."
            )
        self.file_path = file_path
        self.rrng: dict = {
            "ionnames": [],
            "ranges": {},
            "ions": {},
            "molecular_ions": [],
        }
        self.unique = unique
        self.verbose = verbose
        self.read_rrng()

    def read_rrng(self):
        """Read content of an RRNG range file."""
        with open(self.file_path, encoding="utf8") as rrngf:
            txt = rrngf.read()

        txt = txt.replace("\r\n", "\n")  # windows to unix EOL conversion
        txt = txt.replace(",", ".")  # use decimal dots instead of comma
        txt_stripped = [
            line
            for line in txt.split("\n")
            if line.strip() != "" and line.startswith("#") is False
        ]
        del txt

        # see DOI: 10.1007/978-1-4899-7430-3 for further details to this
        # AMETEK/Cameca"s *.rrng file format

        # first, parse [Ions] section, which holds a list of element names
        # there are documented cases where experimentalists add custom strings
        # to specify ranges they consider special
        # these are loaded as user types
        # with nuclide_hash np.iinfo(np.uint16).max
        where = [idx for idx, element in enumerate(txt_stripped) if element == "[Ions]"]
        if not isinstance(where, list):
            raise ValueError("Section [Ions] not found.")
        if len(where) != 1:
            raise ValueError("Section [Ions] not found or ambiguous.")
        current_line_id = where[0] + 1

        tmp = re.split(r"[\s=]+", txt_stripped[current_line_id])
        if len(tmp) != 2:
            raise ValueError(
                f"Line {txt_stripped[current_line_id]} [Ions]/Number line corrupted."
            )
        if tmp[0] != "Number":
            raise ValueError(
                f"Line {txt_stripped[current_line_id]} [Ions]/Number incorrectly formatted."
            )
        if not tmp[1].isnumeric():
            raise ValueError(
                f"Line {txt_stripped[current_line_id]} [Ions]/Number not a number."
            )
        number_of_ion_names = int(tmp[1])
        if number_of_ion_names <= 0:
            raise ValueError(
                f"Line {txt_stripped[current_line_id]} no ion names defined."
            )
        current_line_id += 1
        for i in np.arange(0, number_of_ion_names):
            tmp = re.split(r"[\s=]+", txt_stripped[current_line_id + i])
            if len(tmp) != 2:
                raise ValueError(
                    f"Line {txt_stripped[current_line_id + i]} [Ions]/Ion line corrupted."
                )
            if tmp[0] != f"Ion{i + 1}":
                raise ValueError(
                    f"Line {txt_stripped[current_line_id + i]} [Ions]/Ion incorrectly formatted."
                )
            if not isinstance(tmp[1], str):
                raise ValueError(
                    f"Line {txt_stripped[current_line_id + i]} [Ions]/Name not a string."
                )
            self.rrng["ionnames"].append(tmp[1])

        # second, parse [Ranges] section
        where = [
            idx for idx, element in enumerate(txt_stripped) if element == "[Ranges]"
        ]
        if not isinstance(where, list):
            raise ValueError("Section [Ranges] not found.")
        if len(where) != 1:
            raise ValueError("Section [Ranges] not found or ambiguous.")
        current_line_id = where[0] + 1

        tmp = re.split(r"[\s=]+", txt_stripped[current_line_id])
        if len(tmp) != 2:
            raise ValueError(
                f"Line {txt_stripped[current_line_id]} [Ranges]/Number line corrupted."
            )
        if tmp[0] != "Number":
            raise ValueError(
                f"Line {txt_stripped[current_line_id]} [Ranges]/Number incorrectly formatted."
            )
        if not tmp[1].isnumeric():
            raise ValueError(
                f"Line {txt_stripped[current_line_id]} [Ranges]/Number not a number."
            )
        number_of_ranges = int(tmp[1])
        if number_of_ranges <= 0:
            raise ValueError(
                f"Line {txt_stripped[current_line_id]}  No ranges defined."
            )
        current_line_id += 1

        m_ions = []
        for jdx in np.arange(0, number_of_ranges):
            if self.verbose:
                logger.debug(f"{txt_stripped[current_line_id + jdx]}")
            dct = evaluate_rrng_range_line(jdx + 1, txt_stripped[current_line_id + jdx])
            if dct is None:
                logger.warning(
                    f"RRNG line {txt_stripped[current_line_id + jdx]} is corrupted."
                )
                continue

            m_ion = NxIon(
                nuclide_hash=create_nuclide_hash(dct["atoms"]), charge_state=0
            )
            m_ion.add_range(dct["range"][0], dct["range"][1])
            m_ion.comment = dct["name"]
            m_ions.append(m_ion)
            # this set may contain duplicates or overlapping ranges if ranging definitions are ambiguous like here https://doi.org/10.5281/zenodo.7788883

        if self.unique:
            unique_m_ions = try_to_reduce_to_unique_definitions(m_ions)
            logger.info(
                f"Found {len(m_ions)} ranging definitions, performed reduction to {len(unique_m_ions)} unique ones."
            )
        else:
            unique_m_ions = m_ions.copy()
            logger.info(
                f"Found {len(m_ions)} ranging definitions, no reduction, {len(unique_m_ions)} remain."
            )
        del m_ions

        for m_ion in unique_m_ions:
            m_ion.apply_combinatorics()
            # m_ion.report()
            self.rrng["molecular_ions"].append(m_ion)
        logger.info(f"{self.file_path} parsed successfully.")
