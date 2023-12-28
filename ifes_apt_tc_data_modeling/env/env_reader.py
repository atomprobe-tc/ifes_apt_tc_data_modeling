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


import re
import numpy as np
from ase.data import chemical_symbols
from ifes_apt_tc_data_modeling.nexus.nx_ion import NxField, NxIon
from ifes_apt_tc_data_modeling.utils.utils import \
    create_isotope_vector, is_range_significant
from ifes_apt_tc_data_modeling.utils.definitions import MQ_EPSILON
from ifes_apt_tc_data_modeling.utils.molecular_ions import MolecularIonBuilder
from ifes_apt_tc_data_modeling.utils.molecular_ions import \
    PRACTICAL_ABUNDANCE, PRACTICAL_ABUNDANCE_PRODUCT, \
    PRACTICAL_MIN_HALF_LIFE, VERBOSE, SACRIFICE_ISOTOPIC_UNIQUENESS


def get_smart_chemical_symbols():
    """Organize element symbols such that search H does not match He."""
    priority_queue = []
    for symbol in chemical_symbols:
        if len(symbol) == 2:
            priority_queue.append(symbol)
    for symbol in chemical_symbols:
        if symbol != "X" and len(symbol) == 1:
            priority_queue.append(symbol)
    return priority_queue


def evaluate_env_range_line(line: str) -> dict:
    """Represent information content of a single range line."""
    # example line: ". 107.7240 108.0960 1 0 0 0 0 0 0 0 0 0 3 0 0 0"
    info: dict = {}
    info["identifier"] = None
    info["range"] = np.asarray([0., MQ_EPSILON], np.float64)
    info["atoms"] = []
    info["volume"] = np.float64(0.)
    info["color"] = ""
    info["name"] = ""

    tmp = line.split()

    # interpret zeroth token into a list of chemical symbols
    # interpret first token as inclusive left of m/q interval
    # interpret second token as inclusive right bound of m/q interval
    if len(tmp) < 3:
        print(f"WARNING::ENV file ranging definition {line} has insufficient information!")
        return None
    if is_range_significant(np.float64(tmp[1]), np.float64(tmp[2])) is True:
        info["range"] = np.asarray([tmp[1], tmp[2]], np.float64)
    else:
        print(f"WARNING::ENV file ranging definition {line} has insignificant range!")
        return None

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

    def __init__(self, filename: str):
        if (len(filename) <= 4) or (filename.lower().endswith(".env") is False):
            raise ImportError("WARNING::ENV file incorrect filename ending or file type!")
        self.filename = filename
        self.env: dict = {}
        self.env["ranges"] = {}
        self.env["ions"] = {}
        self.env["molecular_ions"] = []
        self.read_env()

    def read_env(self):
        """Read ENV system configuration and ranging definitions."""
        # GPM/Rouen ENV file format is neither standardized nor uses magic number
        with open(self.filename, mode="r", encoding="utf-8") as envf:
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

                m_ion = NxIon(isotope_vector=create_isotope_vector(
                    dct["atoms"]), charge_state=0)
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
                recovered_charge_state, m_ion_candidates = crawler.combinatorics(
                    m_ion.isotope_vector.typed_value,
                    m_ion.ranges.typed_value[0, 0],
                    m_ion.ranges.typed_value[0, 1])
                # print(f"{recovered_charge_state}")
                m_ion.charge_state = NxField(np.int8(recovered_charge_state), "")
                m_ion.update_human_readable_name()
                m_ion.add_charge_state_model({"min_abundance": PRACTICAL_ABUNDANCE,
                                              "min_abundance_product": PRACTICAL_ABUNDANCE_PRODUCT,
                                              "min_half_life": PRACTICAL_MIN_HALF_LIFE,
                                              "sacrifice_isotopic_uniqueness": SACRIFICE_ISOTOPIC_UNIQUENESS},
                                             m_ion_candidates)

                self.env["molecular_ions"].append(m_ion)
            print(f"{self.filename} parsed successfully")
