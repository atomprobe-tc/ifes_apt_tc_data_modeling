# Utilities for parsing data and molecular ions in atom probe microscopy.
#
# Also convenience functions are included which translate human-readable ion
# names into the isotope_vector description proposed by Kuehbach et al. in
# DOI: 10.1017/S1431927621012241 to the human-readable ion names which are use
# in P. Felfer et al.'s atom probe toolbox
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

import typing

from typing import Tuple

import mmap

import numpy as np

from ase.data import atomic_numbers, chemical_symbols

from ifes_apt_tc_data_modeling.utils.nist_isotope_data import isotopes

from ifes_apt_tc_data_modeling.utils.definitions import MAX_NUMBER_OF_ION_SPECIES
from ifes_apt_tc_data_modeling.utils.definitions import MAX_NUMBER_OF_ATOMS_PER_ION
from ifes_apt_tc_data_modeling.utils.definitions import MQ_EPSILON


def isotope_to_hash(proton_number: int = 0,
                 neutron_number: int = 0) -> int:
    """Encode an isotope to a hashvalue."""
    n_protons = np.uint16(proton_number)
    n_neutrons = np.uint16(neutron_number)
    assert (n_protons >= np.uint16(0)) and (n_protons < np.uint16(256)), \
        "Argument proton number on [0, 256) needed!"
    assert (n_neutrons >= np.uint16(0)) and (n_neutrons < np.uint16(256)), \
        "Argument neutron number on [0, 256) needed!"
    hashvalue = int(n_protons + (np.uint16(256) * n_neutrons))
    return hashvalue


def hash_to_isotope(hashvalue: int = 0) -> Tuple[int, int]:
    """Decode a hashvalue to an isotope."""
    # assert isinstance(hashvalue, int), \
    #     "Argument hashvalue needs to be integer!"
    val = np.uint16(hashvalue)
    assert val >= np.uint16(0), \
        "Argument hashvalue needs to be an unsigned integer!"
    assert val <= np.iinfo(np.uint16).max, \
        "Argument hashvalue needs to map on an uint16!"
    neutron_number = np.uint16(val / np.uint16(256))
    proton_number = np.uint16(val - neutron_number * np.uint16(256))
    return (int(proton_number), int(neutron_number))


def create_isotope_vector(building_blocks: list) -> np.ndarray:
    """Create specifically-shaped array of isotope hashvalues."""
    # building_blocks are usually names of elements in the periodic tables
    # if not we assume the ion is special, a user type

    # test cases:
    # create_isotope_vector(["Fe", "Fe", "O", "O", "O"])
    symbol_to_proton_number = atomic_numbers

    ivec = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,), np.uint16)
    hashvector = []
    assert len(building_blocks) <= MAX_NUMBER_OF_ATOMS_PER_ION, \
        "Faced an ion with an unsupported high complexity!"
    # MAX_NUMBER_OF_ATOMS_PER_ION can be modified to describe large fragments

    if len(building_blocks) == 0:  # special case unknown ion type
        return ivec

    for block in building_blocks:
        assert isinstance(block, str), \
            "Argument block has to be a string, symbol for an element like K or isotope K-40!"
        assert block != "", \
            "Argument block has to be a non-empty string!"
        if block.count("-") == 0:
            assert (block in symbol_to_proton_number) and (block != "X"), \
                block + " is not a valid chemical symbol!"
            proton_number = symbol_to_proton_number[block]
            neutron_number = 0
            hashvector.append(isotope_to_hash(proton_number, neutron_number))
        elif block.count("-") == 1:
            symb_mass = block.split("-")
            assert len(symb_mass) == 2, \
                block + " is not properly formatted <symbol>-<mass_number>!"
            assert (symb_mass[0] in symbol_to_proton_number) and (symb_mass[0] != "X"), \
                symb_mass[0] + " is not a valid chemical symbol!"
            proton_number = symbol_to_proton_number[symb_mass[0]]
            mass_number = int(symb_mass[1])
            neutron_number = mass_number - proton_number
            assert proton_number in isotopes.keys(), \
                "No isotopes for proton_number " + str(proton_number) + " via ase!"
            assert mass_number in isotopes[proton_number].keys(), \
                "No isotope for mass_number " + str(mass_number) + " via ase!"
            hashvector.append(isotope_to_hash(proton_number, neutron_number))
        else:
            print(f"WARNING: {block} does not specify a unique element name!")
            return ivec

    ivec[0:len(hashvector)] = np.sort(
        np.asarray(hashvector, np.uint16), kind="stable")[::-1]
    return ivec


def isotope_vector_to_nuclid_list(ivec: np.ndarray) -> np.ndarray:
    """Create a NeXus NXion nuclid list."""
    assert np.shape(ivec) == (MAX_NUMBER_OF_ATOMS_PER_ION,), \
        "Argument isotope_vector needs to be shaped (MAX_NUMBER_OF_ATOMS_PER_ION,) !"
    nuclid_list = np.zeros((2, MAX_NUMBER_OF_ATOMS_PER_ION), np.uint16)
    for idx in np.arange(0, MAX_NUMBER_OF_ATOMS_PER_ION):
        if ivec[idx] != 0:
            protons, neutrons = hash_to_isotope(int(ivec[idx]))
            if neutrons != 0:
                nuclid_list[0, idx] = protons + neutrons
            # else: leave only element name, i.e. number of protons
            nuclid_list[1, idx] = protons
    return nuclid_list


def isotope_vector_to_dict_keyword(ivec: np.ndarray) -> str:
    """Create keyword for dictionary from isotope_vector."""
    assert len(ivec) <= MAX_NUMBER_OF_ATOMS_PER_ION, \
        "Argument isotope_vector len <= MAX_NUMBER_OF_ATOMS_PER_ION !"
    lst = []
    for hashvalue in ivec:
        if hashvalue != 0:
            lst.append(str(hashvalue))
    if lst != []:
        return "_".join(lst)
    return "0"  # "_".join(np.asarray(np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,)), np.uint16))


def isotope_vector_to_human_readable_name(ivec: np.ndarray, charge: np.int8) -> str:
    """Get human-readable name from an isotope_vector."""
    assert len(ivec) <= MAX_NUMBER_OF_ATOMS_PER_ION, \
        "Argument isotope_vector len <= MAX_NUMBER_OF_ATOMS_PER_ION !"
    human_readable = ""
    if np.sum(ivec) != 0:
        for hashvalue in ivec:
            if hashvalue != 0:
                protons, neutrons = hash_to_isotope(int(hashvalue))
                if neutrons > 0:
                    human_readable += str(protons + neutrons) \
                        + chemical_symbols[protons]
                else:
                    human_readable += chemical_symbols[protons]
                human_readable += " "
        if charge > 0 and charge < 8:
            human_readable += "+" * charge
        elif charge > -8 and charge < 0:
            human_readable += "-" * -charge
        else:
            human_readable = human_readable.rstrip()
        return human_readable
    else:
        return "unknown_iontype"


def is_range_overlapping(interval: np.ndarray,
                        interval_set: np.float64) -> bool:
    """Check if interval overlaps within with members of interval set."""
    assert np.shape(interval) == (2,), "Interval needs to have two columns!"
    assert np.shape(interval_set)[1] == 2, \
        "Interval_set needs to have two columns!"
    # interval = np.array([53.789, 54.343])
    # interval_set = np.array([[27.778, 28.33]])  # for testing purposes
    if np.shape(interval_set)[0] >= 1:
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
    assert left >= np.float64(0.) and right >= np.float64(0.), \
        "Left and right bound have to be positive!"
    if (right - left) > MQ_EPSILON:
        return True
    return False
