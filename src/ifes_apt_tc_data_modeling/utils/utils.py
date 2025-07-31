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

from typing import Tuple
import re
import numpy as np

from ase.data import atomic_numbers, chemical_symbols
from ifes_apt_tc_data_modeling.utils.nist_isotope_data import isotopes
from ifes_apt_tc_data_modeling.utils.definitions import (
    MAX_NUMBER_OF_ATOMS_PER_ION,
    MQ_EPSILON,
    NEUTRON_NUMBER_FOR_ELEMENT,
)
from ifes_apt_tc_data_modeling.utils.custom_logging import logger


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


def isotope_to_hash(
    proton_number: int = 0, neutron_number: int = NEUTRON_NUMBER_FOR_ELEMENT
) -> int:
    """Encode an isotope to a hashvalue."""
    # the acceptance of NXapm introduced a breaking change for this mapping
    # previously meaning the element, i.e., irrespective its isotopes used 0
    #   this case is degenerate though for the 1H isotope and hydrogen element
    # the new mapping therefore uses neutron_number = NEUTRON_NUMBER_FOR_ELEMENT
    # as the offset to, highlight when elements are meant instead of specific isotopes, i.e.
    # with NEUTRON_NUMBER_FOR_ELEMENT == 255
    # 1H is mapped to 1 + 0 * 256, while hydrogen as an element is mapped to 1 + 255 * 256
    if (0 <= proton_number < 256) and (0 <= neutron_number < 256):
        return int(
            np.uint16(proton_number) + (np.uint16(256) * np.uint16(neutron_number))
        )
    return 0


def hash_to_isotope(hashvalue: int = 0) -> Tuple[int, int]:
    """Decode a hashvalue to an isotope."""
    # see comment on isotope_to_hash
    if 0 <= hashvalue <= int(np.iinfo(np.uint16).max):
        neutron_number = np.uint16(np.uint16(hashvalue) / np.uint16(256))
        proton_number = np.uint16(
            np.uint16(hashvalue) - neutron_number * np.uint16(256)
        )
        return (int(proton_number), int(neutron_number))
    return (0, 0)


def create_nuclide_hash(building_blocks: list) -> np.ndarray:
    """Create specifically-shaped array of isotope hashvalues."""
    # building_blocks are usually names of elements in the periodic table
    # if not we assume the ion is special such as user type or plain words
    # a typical expected test case is
    # create_nuclide_hash(["Fe", "Fe", "O", "O", "O"])
    ivec = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,), np.uint16)
    if 0 < len(building_blocks) <= MAX_NUMBER_OF_ATOMS_PER_ION:
        symbol_to_proton_number = atomic_numbers
        hashvector = []
        for block in building_blocks:
            if isinstance(block, str) and block != "":
                if block.count("-") == 0:  # an element
                    if (block not in symbol_to_proton_number) or (block == "X"):
                        return ivec
                    hashvector.append(
                        isotope_to_hash(
                            symbol_to_proton_number[block], NEUTRON_NUMBER_FOR_ELEMENT
                        )
                    )
                elif block.count("-") == 1:
                    symb_mass = block.split("-")
                    if (
                        (len(symb_mass) != 2)
                        or (symb_mass[0] not in symbol_to_proton_number)
                        or (symb_mass[0] == "X")
                    ):
                        logger.warning(
                            f"{block} is not properly formatted <symbol>-<mass_number>"
                        )
                        return ivec
                    proton_number = symbol_to_proton_number[symb_mass[0]]
                    mass_number = int(symb_mass[1])
                    neutron_number = mass_number - proton_number
                    if (proton_number in isotopes) and (
                        mass_number in isotopes[proton_number]
                    ):
                        hashvector.append(
                            isotope_to_hash(proton_number, neutron_number)
                        )
        ivec[0 : len(hashvector)] = np.sort(
            np.asarray(hashvector, np.uint16), kind="stable"
        )[::-1]
    return ivec


def nuclide_hash_to_nuclide_list(ivec: np.ndarray) -> np.ndarray:
    """Create a NeXus NXion nuclide list."""
    nuclide_list = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION, 2), np.uint16)
    if np.shape(ivec) == (MAX_NUMBER_OF_ATOMS_PER_ION,):
        for idx in np.arange(0, MAX_NUMBER_OF_ATOMS_PER_ION):
            if ivec[idx] != 0:
                n_protons, n_neutrons = hash_to_isotope(int(ivec[idx]))
                if 0 < n_neutrons < 255:
                    # right bound must not be NEUTRON_NUMBER_FOR_ELEMENT but 255
                    # see breaking change isotope_to_hash
                    nuclide_list[idx, 0] = n_protons + n_neutrons
                nuclide_list[idx, 1] = n_protons
        return nuclide_list
    logger.warning(
        f"Argument nuclide_hash needs to be shaped ({MAX_NUMBER_OF_ATOMS_PER_ION},)"
    )
    return nuclide_list


def nuclide_hash_to_dict_keyword(ivec: np.ndarray) -> str:
    """Create keyword for dictionary from nuclide_hash."""
    if len(ivec) <= MAX_NUMBER_OF_ATOMS_PER_ION:
        lst = []
        for hashvalue in ivec:
            if hashvalue != 0:
                lst.append(f"{hashvalue}")
        if len(lst) > 0:
            return "_".join(lst)
    return (
        "0"  # "_".join(np.asarray(np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,)), np.uint16))
    )


def nuclide_hash_to_human_readable_name(ivec: np.ndarray, charge_state: np.int8) -> str:
    """Get human-readable name from an nuclide_hash."""
    if len(ivec) <= MAX_NUMBER_OF_ATOMS_PER_ION and -8 < charge_state < 8:
        human_readable = ""
        for hashvalue in ivec:
            if hashvalue != 0:
                protons, neutrons = hash_to_isotope(int(hashvalue))
                if 0 < neutrons < 255:
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
        if human_readable != "":
            return human_readable
    return "unknown_iontype"


def is_range_overlapping(interval: np.ndarray, interval_set: np.ndarray) -> bool:
    """Check if interval overlaps within with members of interval set."""
    if (
        (np.shape(interval) == (2,))
        and (np.shape(interval_set)[0] >= 1)
        and (np.shape(interval_set)[1] == 2)
    ):
        # interval = np.array([53.789, 54.343])
        # interval_set = np.array([[27.778, 28.33]])  # for testing purposes
        left_and_right_delta = np.zeros([np.shape(interval_set)[0], 2], bool)
        left_and_right_delta[:, 0] = (interval_set[:, 0] - interval[1]) > MQ_EPSILON
        left_and_right_delta[:, 1] = (interval[0] - interval_set[:, 1]) > MQ_EPSILON
        no_overlap = np.array([np.any(i) for i in left_and_right_delta])
        if np.all(no_overlap):
            return False
        return True
    return False


def is_range_significant(left: np.float64, right: np.float64) -> bool:
    """Check if inclusive interval bounds [left, right] span a finite range."""
    if (np.float64(0.0) <= left) and (np.float64(0.0) <= right):
        if (right - left) >= MQ_EPSILON:
            return True
    return False


def is_convertible_to_isotope_hash(symbol: str):
    """Check if human_readable symbol is convertible into nuclide hash tribool."""
    case = re.search("^([A-Z])([a-z])?(-)([0-9]+)$", symbol)
    if case is None:  # eventually element case e.g. "K"
        if not isinstance(symbol, str):
            raise ValueError("Argument symbol needs to be a string !")
        if symbol not in get_smart_chemical_symbols():
            raise ValueError(f"Symbol needs to be in {get_smart_chemical_symbols()}!")
        return 1
    # alternative case eventually specific nuclide e.g. "K-40"
    symb_mass = symbol.split("-")
    if len(symb_mass) != 2:
        raise TypeError(
            "Argument symbol is not properly formatted <symbol>-<mass_number>!"
        )
    if len(symb_mass[0]) != 1 and len(symb_mass[0]) != 2:
        raise ValueError(
            "Argument symbol is not properly formatted <symbol>-<mass_number>!"
        )
    if len(symb_mass[1]) <= 0:
        raise ValueError(
            f"Argument symbol {symb_mass[1]} needs to be a physical mass number!"
        )
    if symb_mass[0] not in get_smart_chemical_symbols():
        raise ValueError(
            f"{symb_mass[0]} is not a symbol in {get_smart_chemical_symbols()}!"
        )
    if int(symb_mass[1]) not in isotopes[atomic_numbers[symb_mass[0]]].keys():
        raise ValueError(
            f"No value for isotopes[atomic_numbers[{symb_mass[0]}][{int(symb_mass[1])}] exists!"
        )
    return 2


def element_or_nuclide_to_hash(symbol: str):
    """Converts an element symbol (e.g. K) or nuclide (K-40) to nuclide hash."""
    # consider moving this to the ifes_apt_tc_data_modeling library
    case = is_convertible_to_isotope_hash(symbol)
    if case == 1:
        return isotope_to_hash(atomic_numbers[symbol], NEUTRON_NUMBER_FOR_ELEMENT)
    if case == 2:
        symb_mass = symbol.split("-")
        return isotope_to_hash(
            atomic_numbers[symb_mass[0]],
            int(symb_mass[1]) - atomic_numbers[symb_mass[0]],
        )
    return isotope_to_hash(0, 0)


def symbol_lst_to_matrix_of_nuclide_vector(
    method: str = "resolve_all", symbol_lst: list = [], **kwargs
):
    """Parse human-readable element/isotope names to paraprobe."""
    # method = "resolve_all"
    # method = "resolve_unknown"
    # method = "resolve_element"
    # method = "resolve_isotope"
    # method = "resolve_ion"
    # symbol_lst = [["K", "C"], ["U-238", "H-2"]]
    # symbol_lst = [["K-40", "C-14"], ["U-238", "H-2"]]  # one list per mol. ion
    # kwargs_charge_lst = [+1, +3]  # one charge per ion
    # resolve_all and resolve_unknown are the special cases
    # where we do not need to specify a matrix of isotope vectors
    if method in ("resolve_all", "resolve_unknown"):
        return (None, None, None)

    supported = ("resolve_element", "resolve_ion", "resolve_isotope")
    if method not in supported:
        raise ValueError(f"Argument method needs to be in {supported} !")
    if not isinstance(symbol_lst, list) and all(
        isinstance(lst, list) for lst in symbol_lst
    ):
        raise TypeError("Argument symbol_lst must be a list of lists!")
    if not all(len(np.shape(lst)) == 1 for lst in symbol_lst):
        raise ValueError(
            "One list in argument symbol_lst is not a 1d list or an empty list!"
        )
    if not all(np.shape(lst)[0] >= 1 for lst in symbol_lst):
        raise ValueError(
            "One list in argument symbol_lst is not a 1d list or an empty list!"
        )
    matrix = np.zeros([len(symbol_lst), MAX_NUMBER_OF_ATOMS_PER_ION])
    charge = []
    if (method == "resolve_ion") and ("charge_lst" in kwargs):
        if not isinstance(kwargs["charge_lst"], list):
            raise TypeError("Keyword argument charge_lst must be a list of lists!")
        if not all(isinstance(val, int) for val in kwargs["charge_lst"]):
            raise ValueError("Keyword argument charge_lst needs to be a list of int !")
        if np.shape(symbol_lst)[0] != np.shape(kwargs["charge_lst"])[0]:
            raise ValueError(
                "Argument symbol_lst and keyword argument charge_lst need to have the same length !"
            )
    for idx, lst in enumerate(symbol_lst):
        ivec = np.zeros([1, MAX_NUMBER_OF_ATOMS_PER_ION])
        if lst == []:
            raise ValueError("Argument molecular ion must not be an empty list!")
        if len(lst) > MAX_NUMBER_OF_ATOMS_PER_ION:
            raise ValueError(
                f"Argument molecular ion must not contain more than "
                f"{MAX_NUMBER_OF_ATOMS_PER_ION} entries!"
            )
        jdx = 0
        for symbol in lst:
            if method == "resolve_element":
                if symbol in atomic_numbers:
                    ivec[0, jdx] = isotope_to_hash(
                        atomic_numbers[symbol], NEUTRON_NUMBER_FOR_ELEMENT
                    )  # do not encode isotope information
                    jdx += 1
                else:
                    # it might be that we have a specific nuclide e.g. K-40 try to parse element symbol
                    candidate = symbol.split("-", 1)[0]
                    if candidate not in atomic_numbers:
                        raise KeyError(
                            f"symbol_lst[{idx}] candidate does not specify an element!"
                        )
                    ivec[0, jdx] = isotope_to_hash(
                        atomic_numbers[candidate], NEUTRON_NUMBER_FOR_ELEMENT
                    )
                    jdx += 1
            else:  # "resolve_isotope", "resolve_ion":
                ivec[0, jdx] = element_or_nuclide_to_hash(symbol)
                jdx += 1
        matrix[idx, :] = ivec[0, :]
        if method == "resolve_ion":
            charge.append(kwargs["charge_lst"][idx])
    if charge == []:
        return (method, np.asarray(matrix, np.uint16), None)
    return (method, np.asarray(matrix, np.uint16), np.asarray(charge, np.int8))
