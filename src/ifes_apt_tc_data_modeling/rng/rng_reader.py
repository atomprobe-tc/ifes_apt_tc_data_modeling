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

"""RNG range file reader used by atom probe microscopists."""

# pylint: disable=duplicate-code

import re
import numpy as np

from ifes_apt_tc_data_modeling.nexus.nx_ion import NxIon
from ifes_apt_tc_data_modeling.utils.utils import (
    create_nuclide_hash,
    is_range_significant,
)
from ifes_apt_tc_data_modeling.utils.definitions import MQ_EPSILON
from ifes_apt_tc_data_modeling.utils.molecular_ions import get_chemical_symbols
from ifes_apt_tc_data_modeling.utils.custom_logging import logger


# there are specific examples for unusual range files here:
# https://hg.sr.ht/~mycae/libatomprobe/browse/test/samples/ranges?rev=tip


def evaluate_rng_range_line(
    i: int, line: str, column_id_to_label: dict, n_columns: int
) -> dict:
    """Represent information content of a single range line."""
    # example line: ". 107.7240 108.0960 1 0 0 0 0 0 0 0 0 0 3 0 0 0"
    info: dict = {
        "identifier": f"Range{i}",
        "range": np.asarray([0.0, MQ_EPSILON], np.float64),
        "atoms": [],
        "volume": np.float64(0.0),
        "color": "",
        "name": "",
    }

    tmp = re.split(r"\s+", line)
    if len(tmp) != n_columns:
        raise ValueError(f"Line {line} inconsistent number columns {len(tmp)}.")
    if tmp[0] != ".":
        raise ValueError(f"Line {line} has inconsistent line prefix.")
    if not is_range_significant(np.float64(tmp[1]), np.float64(tmp[2])):
        # raise ValueError(f"Line {line} insignificant range.")
        return info
    info["range"] = np.asarray([tmp[1], tmp[2]], np.float64)

    # line encodes multiplicity of element via array of multiplicity counts
    element_multiplicity = np.asarray(tmp[3 : len(tmp)], np.uint32)
    if np.sum(element_multiplicity) < 0:
        # raise ValueError(f"Line {line} no element counts.")
        return info
    if np.sum(element_multiplicity) > 0:
        for jdx in np.arange(0, len(element_multiplicity)):
            if element_multiplicity[jdx] < 0:
                # raise ValueError(f"Line {line} no negative element counts.")
                raise ValueError(
                    f"element_multiplicity[jdx] {element_multiplicity[jdx]} needs to be positive."
                )
            if element_multiplicity[jdx] > 0:
                symbol = column_id_to_label[jdx + 1]
                if symbol in get_chemical_symbols():
                    info["atoms"] = np.append(
                        info["atoms"],
                        [column_id_to_label[jdx + 1]] * int(element_multiplicity[jdx]),
                    )
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
    info: dict = {"column_id_to_label": {}}
    tmp = re.split(r"\s+", line)
    if len(tmp) == 0:
        raise ValueError(f"Line {line} does not contain iontype labels {len(tmp)}.")
    for idx in np.arange(1, len(tmp)):
        info["column_id_to_label"][idx] = tmp[idx]
    return info


class ReadRngFileFormat:
    """Read *.rng file format."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 4) or not file_path.lower().endswith(".rng"):
            raise ImportError(
                "WARNING::RNG file incorrect file_path ending or file type."
            )
        self.file_path = file_path
        self.rng: dict = {"ranges": {}, "ions": {}, "molecular_ions": []}
        self.read_rng()

    def read_rng(self):
        """Read RNG range file content."""
        with open(self.file_path, mode="r", encoding="utf8") as rngf:
            txt = rngf.read()

        txt = txt.replace("\r\n", "\n")  # windows to unix EOL conversion
        txt = txt.replace(",", ".")  # use decimal dots instead of comma
        txt_stripped = [
            line
            for line in txt.split("\n")
            if line.strip() != "" and line.startswith("#") is False
        ]
        del txt

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
        if tmp is None:
            raise ValueError("RNG file does not contain key header line.")

        header = evaluate_rng_ion_type_header(txt_stripped[current_line_id])

        tmp = re.split(r"\s+", txt_stripped[0])
        if not tmp[0].isnumeric():
            raise ValueError(f"Line {txt_stripped[0]} number of species corrupted.")
        n_element_symbols = int(tmp[0])
        if n_element_symbols < 0:
            raise ValueError(f"Line {txt_stripped[0]} no species defined.")
        if not tmp[1].isnumeric():
            raise ValueError(f"Line {txt_stripped[0]} number of ranges corrupted.")
        n_ranges = int(tmp[1])
        if n_ranges < 0:
            raise ValueError(f"Line {txt_stripped[0]} no ranges defined.")

        for idx in np.arange(current_line_id + 1, current_line_id + 1 + n_ranges):
            dct = evaluate_rng_range_line(
                idx - current_line_id,
                txt_stripped[idx],
                header["column_id_to_label"],
                n_element_symbols + 3,
            )
            if dct is None:
                logger.warning(f"RNG line {txt_stripped[idx]} is corrupted.")
                continue

            m_ion = NxIon(
                nuclide_hash=create_nuclide_hash(dct["atoms"]), charge_state=0
            )
            m_ion.add_range(dct["range"][0], dct["range"][1])
            m_ion.comment = dct["name"]
            m_ion.apply_combinatorics()
            # m_ion.report()

            self.rng["molecular_ions"].append(m_ion)
        logger.info(f"{self.file_path} parsed successfully.")
