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

"""ENV file format reader for GPM/Rouen ENV system configuration and range files"""

# pylint: disable=too-many-nested-blocks,duplicate-code

import re
import numpy as np
from ifes_apt_tc_data_modeling.nexus.nx_ion import NxField, NxIon
from ifes_apt_tc_data_modeling.utils.utils import \
    create_nuclide_hash, is_range_significant, get_smart_chemical_symbols
from ifes_apt_tc_data_modeling.utils.definitions import MQ_EPSILON


def evaluate_env_range_line(line: str):
    """Represent information content of a single range line."""
    # example line: ". 107.7240 108.0960 1 0 0 0 0 0 0 0 0 0 3 0 0 0"
    info: dict = {"identifier": None,
                  "range": np.asarray([0., MQ_EPSILON], np.float64),
                  "atoms": [],
                  "volume": np.float64(0.),
                  "color": "",
                  "name": ""}

    tmp = line.split()
    # interpret zeroth token into a list of chemical symbols
    # interpret first token as inclusive left of m/q interval
    # interpret second token as inclusive right bound of m/q interval
    if len(tmp) < 3:
        print(f"WARNING::ENV file ranging definition {line} has insufficient information!")
        return None
    if is_range_significant(np.float64(tmp[1]), np.float64(tmp[2])) is False:
        print(f"WARNING::ENV file ranging definition {line} has insignificant range!")
        return None

    info["range"] = np.asarray([tmp[1], tmp[2]], np.float64)
    lst: list = []
    if tmp[0] == "Hyd":
        lst = []
    elif tmp[0] in get_smart_chemical_symbols():
        lst.append(tmp[0])
    else:
        tokens = re.split(r'(\d+)', tmp[0])
        for jdx in np.arange(0, len(tokens)):
            kdx = 0
            for sym in get_smart_chemical_symbols():
                if tokens[jdx][kdx:].startswith(sym) is True:
                    mult = 1
                    if jdx < len(tokens) - 1:
                        if (tokens[jdx][kdx:] == sym) and (tokens[jdx + 1].isdigit() is True):
                            mult = int(tokens[jdx + 1])
                            kdx += len(tokens[jdx + 1])
                    lst.extend([sym] * mult)
                    kdx += len(sym)
    info["atoms"] = lst
    return info


class ReadEnvFileFormat():
    """Read GPM/Rouen *.env file format."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 4) or (file_path.lower().endswith(".env") is False):
            raise ImportError("WARNING::ENV file incorrect file_path ending or file type!")
        self.file_path = file_path
        self.env: dict = {"ranges": {},
                          "ions": {},
                          "molecular_ions": []}
        self.read_env()

    def read_env(self):
        """Read ENV system configuration and ranging definitions."""
        # GPM/Rouen ENV file format is neither standardized nor uses magic number
        with open(self.file_path, mode="r", encoding="utf-8") as envf:
            txt = envf.read()
            txt = txt.replace("\r\n", "\n")  # windows to unix EOL conversion
            txt = txt.replace(",", ".")  # use decimal dots instead of comma
            txt_stripped = [line for line in txt.split("\n") if line.strip() != ""]
            # search for ranging definitions "# Definition of"
            rng_s = None
            rng_e = None
            for idx in np.arange(0, len(txt_stripped)):
                if txt_stripped[idx].startswith("# Definition of") is False:
                    continue
                rng_s = idx
                break
            for idx in np.arange(rng_s + 1, len(txt_stripped)):
                if txt_stripped[idx].startswith("# Atom probe definition") is False:
                    continue
                rng_e = idx
                break
            if rng_s is None or rng_e is None:
                print("WARNING:: No ranging definitions were found!")
                return

            for idx in np.arange(rng_s + 1, rng_e):
                dct = evaluate_env_range_line(txt_stripped[idx])
                if dct is None:
                    continue

                m_ion = NxIon(nuclide_hash=create_nuclide_hash(dct["atoms"]),
                              charge_state=0)
                m_ion.add_range(dct["range"][0], dct["range"][1])
                m_ion.comment = NxField(dct["name"], "")
                m_ion.apply_combinatorics()
                # m_ion.report()

                self.env["molecular_ions"].append(m_ion)
            print(f"{self.file_path} parsed successfully")
