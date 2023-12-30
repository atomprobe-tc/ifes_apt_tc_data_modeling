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

"""Utilities for working with molecular ions in atom probe microscopy."""

# including convenience functions for translating human-readable ion names
# into the isotope_vector description proposed by Kühbach et al. in
# DOI: 10.1017/S1431927621012241

from typing import Tuple
import numpy as np

from ase.data import atomic_numbers, chemical_symbols
from ifes_apt_tc_data_modeling.utils.nist_isotope_data import isotopes
from ifes_apt_tc_data_modeling.utils.definitions import \
    MAX_NUMBER_OF_ATOMS_PER_ION, MQ_EPSILON


def isotope_to_hash(proton_number: int = 0,
                    neutron_number: int = 0) -> int:
    """Encode an isotope to a hashvalue."""
    if (0 <= proton_number < 256) and (0 <= neutron_number < 256):
        return int(np.uint16(proton_number) + (np.uint16(256) * np.uint16(neutron_number)))
    return 0


def hash_to_isotope(hashvalue: int = 0) -> Tuple[int, int]:
    """Decode a hashvalue to an isotope."""
    # assert isinstance(hashvalue, int), \
    #     "Argument hashvalue needs to be integer!"
    if 0 <= hashvalue <= int(np.iinfo(np.uint16).max):
        neutron_number = np.uint16(np.uint16(hashvalue) / np.uint16(256))
        proton_number = np.uint16(np.uint16(hashvalue) - neutron_number * np.uint16(256))
        return (int(proton_number), int(neutron_number))
    return (0, 0)


def create_isotope_vector(building_blocks: list) -> np.ndarray:
    """Create specifically-shaped array of isotope hashvalues."""
    # building_blocks are usually names of elements in the periodic table
    # if not we assume the ion is special such as user type or plain words
    # a typical expected test case is
    # create_isotope_vector(["Fe", "Fe", "O", "O", "O"])
    ivec = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,), np.uint16)
    if 0 < len(building_blocks) <= MAX_NUMBER_OF_ATOMS_PER_ION:
        symbol_to_proton_number = atomic_numbers
        hashvector = []
        for block in building_blocks:
            if isinstance(block, str) and block != "":
                if block.count("-") == 0:  # an element
                    if (block not in symbol_to_proton_number) or (block == "X"):
                        return ivec
                    hashvector.append(isotope_to_hash(
                        symbol_to_proton_number[block], 0))
                elif block.count("-") == 1:
                    symb_mass = block.split("-")
                    if (len(symb_mass) != 2) or (symb_mass[0] not in symbol_to_proton_number) or (symb_mass[0] == "X"):
                        print(f"WARNING:: {block} is not properly formatted <symbol>-<mass_number>!")
                        return ivec
                    proton_number = symbol_to_proton_number[symb_mass[0]]
                    mass_number = int(symb_mass[1])
                    neutron_number = mass_number - proton_number
                    if (proton_number in isotopes) and (mass_number in isotopes[proton_number]):
                        hashvector.append(isotope_to_hash(proton_number, neutron_number))
                    return ivec
                return ivec

        ivec[0:len(hashvector)] = np.sort(
            np.asarray(hashvector, np.uint16), kind="stable")[::-1]
    return ivec


def isotope_vector_to_nuclid_list(ivec: np.ndarray) -> np.ndarray:
    """Create a NeXus NXion nuclid list."""
    nuclid_list = np.zeros((2, MAX_NUMBER_OF_ATOMS_PER_ION), np.uint16)
    if np.shape(ivec) == (MAX_NUMBER_OF_ATOMS_PER_ION,):
        for idx in np.arange(0, MAX_NUMBER_OF_ATOMS_PER_ION):
            if ivec[idx] != 0:
                protons, neutrons = hash_to_isotope(int(ivec[idx]))
                nuclid_list[0, idx] = protons + neutrons
                nuclid_list[1, idx] = protons
        return nuclid_list
    print(f"WARNING:: Argument isotope_vector needs to be "
          f"shaped {MAX_NUMBER_OF_ATOMS_PER_ION},) !")
    return nuclid_list


def isotope_vector_to_dict_keyword(ivec: np.ndarray) -> str:
    """Create keyword for dictionary from isotope_vector."""
    if len(ivec) <= MAX_NUMBER_OF_ATOMS_PER_ION:
        lst = []
        for hashvalue in ivec:
            if hashvalue != 0:
                lst.append(f"{hashvalue}")
        if len(lst) > 0:
            return "_".join(lst)
    return "0"  # "_".join(np.asarray(np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,)), np.uint16))


def isotope_vector_to_human_readable_name(ivec: np.ndarray, charge_state: np.int8) -> str:
    """Get human-readable name from an isotope_vector."""
    if len(ivec) <= MAX_NUMBER_OF_ATOMS_PER_ION:
        human_readable = ""
        for hashvalue in ivec:
            if hashvalue != 0:
                protons, neutrons = hash_to_isotope(int(hashvalue))
                if neutrons > 0:
                    human_readable += f"{protons + neutrons}{chemical_symbols[protons]}"
                else:
                    human_readable += f"{chemical_symbols[protons]}"
                human_readable += " "
        if 0 < charge_state < 8:
            human_readable += "+" * charge_state
        elif -8 < charge_state < 0:
            human_readable += "-" * -charge_state
        else:
            human_readable = human_readable.rstrip()
        return human_readable
    return "unknown_iontype"


def is_range_overlapping(interval: np.ndarray,
                         interval_set: np.ndarray) -> bool:
    """Check if interval overlaps within with members of interval set."""
    if (np.shape(interval) == (2,)) and (np.shape(interval_set)[0] >= 1) \
            and (np.shape(interval_set)[1] == 2):
        # interval = np.array([53.789, 54.343])
        # interval_set = np.array([[27.778, 28.33]])  # for testing purposes
        left_and_right_delta = np.zeros([np.shape(interval_set)[0], 2], bool)
        left_and_right_delta[:, 0] = (interval_set[:, 0] - interval[1]) \
            > MQ_EPSILON
        left_and_right_delta[:, 1] = (interval[0] - interval_set[:, 1]) \
            > MQ_EPSILON
        no_overlap = np.array([np.any(i) for i in left_and_right_delta])
        if np.all(no_overlap):
            return False
        return True
    return False


def is_range_significant(left: np.float64, right: np.float64) -> bool:
    """Check if inclusive interval bounds [left, right] span a finite range."""
    if (np.float64(0.) <= left) and (np.float64(0.) <= right):
        if (right - left) > MQ_EPSILON:
            return True
    return False
