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

"""Utility tool, analyze and combinatorics for (molecular) ions."""

# pylint: disable=too-many-arguments,too-many-instance-attributes,too-many-locals,too-many-nested-blocks,too-many-branches

import numpy as np
import radioactivedecay as rd

from ase.data import atomic_numbers, chemical_symbols
from ifes_apt_tc_data_modeling.utils.utils import \
    isotope_to_hash, hash_to_isotope, nuclide_hash_to_dict_keyword
from ifes_apt_tc_data_modeling.utils.nist_isotope_data import isotopes
from ifes_apt_tc_data_modeling.utils.definitions import \
    PRACTICAL_ABUNDANCE, PRACTICAL_ABUNDANCE_PRODUCT, \
    PRACTICAL_MIN_HALF_LIFE, SACRIFICE_ISOTOPIC_UNIQUENESS
VERBOSE = False


def get_chemical_symbols():
    """"Report only valid chemical symbols"""
    return chemical_symbols[1:]


class MolecularIonCandidate:
    """Define (molecular) ion build from nuclides."""

    def __init__(self,
                 ivec,
                 charge_state=0,
                 mass_sum=0.,
                 nat_abun_prod=0.,
                 min_half_life=np.inf):
        self.nuclide_hash = np.asarray(ivec, np.uint16)
        self.charge_state = np.int8(charge_state)
        self.mass = np.float64(mass_sum)
        self.abundance_product = np.float64(nat_abun_prod)
        self.shortest_half_life = np.float64(min_half_life)

    def unique_keyword(self):
        """Generate unique keyword."""
        keyword = f"{nuclide_hash_to_dict_keyword(np.sort(np.asarray(self.nuclide_hash, np.uint16), kind='stable')[::-1])}__{self.charge_state}"
        return keyword


class MolecularIonBuilder:
    """Class for holding properties of constructed molecular ions."""

    def __init__(self,
                 min_abundance=PRACTICAL_ABUNDANCE,
                 min_abundance_product=PRACTICAL_ABUNDANCE_PRODUCT,
                 min_half_life=PRACTICAL_MIN_HALF_LIFE,
                 sacrifice_uniqueness=SACRIFICE_ISOTOPIC_UNIQUENESS,
                 verbose=VERBOSE):
        self.nuclides = np.asarray([], np.uint16)
        self.element_isotopes = {}
        self.nuclide_mass = {}
        self.nuclide_abundance = {}
        self.nuclide_stable = {}  # observationally stable
        self.nuclide_unclear = {}  # unclear halflife
        self.nuclide_halflife = {}
        self.candidates = []
        self.parms = {"min_abundance": min_abundance,
                      "min_abundance_product": min_abundance_product,
                      "min_half_life": min_half_life,
                      "sacrifice_isotopic_uniqueness": sacrifice_uniqueness,
                      "verbose": verbose}

        for symbol, atomic_number in atomic_numbers.items():
            if symbol != "X":
                # assume that data from ase take preference
                # if half-life data are available in radioactive decay library
                # take these instead, if all fails mark an unclear_half_life == True
                element_isotopes = []
                for mass_number in isotopes[atomic_number]:
                    half_life = np.inf
                    observationally_stable = False
                    unclear_half_life = False

                    # test if half-life data available
                    trial_nuclide_name = f"{symbol}-{mass_number}"
                    try:
                        tmp = rd.Nuclide(trial_nuclide_name)
                    except ValueError:
                        tmp = None
                    if tmp is not None:
                        half_life = tmp.half_life()
                        if np.isinf(half_life):
                            observationally_stable = True
                            # these ions are always taken as they
                            # are most relevant for practical
                            # atom probe experiments
                        else:
                            if half_life < self.parms["min_half_life"]:
                                # ignore practically short living ions
                                continue
                    else:
                        continue
                        # do not consider exotic isotopes with unclear
                        # half-life as they are likely anyway irrelevant
                        # for practical atom probe experiments
                        # half_life = np.nan
                        # unclear_half_life = True

                    # get ase abundance data
                    n_protons = atomic_numbers[symbol]
                    n_neutrons = mass_number - n_protons
                    mass = isotopes[atomic_numbers[symbol]][mass_number]["mass"]
                    abundance = isotopes[atomic_numbers[symbol]][mass_number]["composition"]
                    hashvalue = isotope_to_hash(int(n_protons), int(n_neutrons))
                    if hashvalue != 0:
                        self.nuclides = np.append(self.nuclides, hashvalue)
                        self.nuclide_mass[hashvalue] = np.float64(mass)
                        self.nuclide_abundance[hashvalue] = np.float64(abundance)
                        self.nuclide_stable[hashvalue] = observationally_stable
                        self.nuclide_unclear[hashvalue] = unclear_half_life
                        self.nuclide_halflife[hashvalue] = half_life
                        element_isotopes = np.append(element_isotopes, hashvalue)
                self.element_isotopes[atomic_number] = np.sort(
                    np.asarray(element_isotopes, np.uint16), kind="stable")[::-1]
        self.nuclides = np.sort(self.nuclides, kind="stable")[::-1]
        if self.parms["verbose"] is True:
            print(f"MolecularIonBuilder initialized with {len(self.nuclides)} nuclides")

    def get_element_isotopes(self, hashvalue):
        """List of hashvalues all isotopes of element specified by hashvalue."""
        return self.element_isotopes[hash_to_isotope(hashvalue)[0]]

    def get_isotope_mass_sum(self, nuclide_arr):
        """Evaluate cumulated atomic_mass of isotopes in ivec."""
        # assuming no relativistic effects or other quantum effects
        # mass loss due to charge_state considered insignificant
        mass = 0.
        for hashvalue in nuclide_arr:
            if hashvalue != 0:
                mass += self.nuclide_mass[hashvalue]
        return mass

    def get_natural_abundance_product(self, nuclide_arr):
        """Get natural abundance product."""
        abun_prod = 1.
        for hashvalue in nuclide_arr:
            if hashvalue != 0:
                abun_prod *= self.nuclide_abundance[hashvalue]
        return abun_prod

    def get_shortest_half_life(self, nuclide_arr):
        """Get shortest half life for set of nuclides."""
        min_half_life = self.parms["min_half_life"]
        for hashvalue in nuclide_arr:
            if hashvalue != 0:
                if self.nuclide_halflife[hashvalue] != np.nan:
                    min_half_life = np.min((min_half_life, self.nuclide_halflife[hashvalue]))
                    # if the min_half_life is np.inf than we know that every nuclide in the
                    # molecular ion is observationally stable
                    # if not we know that considering this molecular ion might not be a good
                    # idea or only necessary in very few cases because even if we were
                    # to find such a combination of nuclides at least one would decay in observable
                    # time and this might possibly hint that studying this molecular is tricky
                else:
                    # if we do not information about the half-life it is very likely that
                    # this is an exotic nuclide likely never found in the wild
                    return np.nan
        return min_half_life

    def combinatorics(self, hash_arr, low, high):
        """Combinatorial analysis which (molecular) elements match within [low, high]."""
        # RNG/RRNG/ENV range files store (molecular) ion information for each range
        # BUT neither nuclide not charge state information, here we try to recover
        # both if possible or list all possible combinations
        # correspondingly this allows yield
        # only an nuclide_hash whose hashvalues have ALL in common
        # that the number of neutrons is 0, i.e. their hashvalue is the atomic_number
        # hash_arr is an nuclide_hash/ivec with such hashvalues
        # as an example the molecular ion C:2 H:1 will map to 6, 6, 1
        max_depth = 0  # number of non-zero entries in hash_arr
        for hashvalue in hash_arr:
            if hashvalue != 0:
                max_depth += 1
        if self.parms["verbose"] is True:
            print(f"Maximum recursion depth {max_depth}")
        self.candidates = []
        if self.parms["verbose"] is True:
            print(hash_arr)

        if max_depth > 0:
            depth = 0
            ith_nuclides = self.get_element_isotopes(hash_arr[depth])
            cand_arr_curr = []  # combinatorially add nuclides while recursing deeper
            self.iterate_molecular_ion(
                hash_arr, ith_nuclides, cand_arr_curr,
                depth, max_depth, low, high)
            if self.parms["verbose"] is True:
                print(f"Found {len(self.candidates)} candidates!")
                for obj in self.candidates:
                    print(f"{obj.nuclide_hash}, {obj.charge_state}, {obj.shortest_half_life}")
            return self.try_to_reduce_to_unique_solution()
            # will return a tuple of charge_state and list of relevant_candidates
        return (0, [])

    def iterate_molecular_ion(self,
                              hash_arr, jth_nuclides, cand_arr_prev,
                              i, max_n, low, high):
        """Recursive analysis of combinatorics on molecular ions."""
        if i < (max_n - 1):
            for nuclide in jth_nuclides:
                ixxth_nuclides = self.get_element_isotopes(hash_arr[i + 1])
                cand_arr_curr = np.append(cand_arr_prev, nuclide)
                self.iterate_molecular_ion(
                    hash_arr, ixxth_nuclides, cand_arr_curr,
                    i + 1, max_n, low, high)
        elif i == (max_n - 1):
            for nuclide in jth_nuclides:
                cand_arr_curr = np.append(cand_arr_prev, nuclide)
                # by this design the ivec does not necessarily remain ordered

                new_mass = self.get_isotope_mass_sum(cand_arr_curr)
                new_abun_prod = self.get_natural_abundance_product(cand_arr_curr)
                new_halflife = self.get_shortest_half_life(cand_arr_curr)

                for new_chrg in np.arange(1, 8):
                    mass_to_charge = new_mass / new_chrg
                    if mass_to_charge < low:
                        break
                    # we can break the entire charge state generation here already
                    # instead of continue because already the current mq is left out
                    # of interval [mqmin, mqmax]
                    # and all mq in the next iterations will result in even lower
                    # mass-to-charge you walk to the left increasing your
                    # distance to the left bound
                    if mass_to_charge > high:
                        # can be optimized and broken out of earlier
                        # if testing first chrg == 1
                        # and then chrg == APTMOLECULARION_MAX_CHANGE
                        continue
                        # must not be break here because with adding more charge
                        # we usually walk from right to left eventually into
                        # [low, high] !
                    # molecular ion is within user-specified bounds
                    self.candidates.append(
                        MolecularIonCandidate(cand_arr_curr,
                                              new_chrg,
                                              new_mass,
                                              new_abun_prod,
                                              new_halflife))
                # do not start a new recursion
        else:
            # break the recursion
            return

    def get_relevant(self):
        """Identify relevant candidates."""
        relevant = {}
        for cand in self.candidates:
            if cand.abundance_product >= self.parms["min_abundance_product"]:
                if not np.isnan(cand.shortest_half_life):
                    # don't dare to test np.isnan(cand.shortest_half_life) is False
                    if cand.shortest_half_life >= self.parms["min_half_life"]:
                        keyword = cand.unique_keyword()
                        if keyword not in relevant:
                            relevant[keyword] = cand

        if self.parms["verbose"] is True:
            print(f"Reduced set to {len(relevant.keys())} relevant candidates...")
            for key in relevant:
                print(key)
        relevant_candidates = []
        for key, obj in relevant.items():
            relevant_candidates.append(obj)
        return (relevant, relevant_candidates)

    def try_to_reduce_to_unique_solution(self):
        """Heuristics to identify if current candidates are unique."""
        if self.parms["verbose"] is True:
            print(f"Reduce set of {len(self.candidates)} candidates to a unique...")
        relevant, relevant_candidates = self.get_relevant()
        if len(relevant) == 0:
            if self.parms["verbose"] is True:
                print("WARNING::No relevant candidate meets all criteria!")
                print("WARNING::No solution possible for given criteria!")
            return (0, relevant_candidates)
        if len(relevant) == 1:
            if self.parms["verbose"] is True:
                print("One relevant candidate which meets all criteria")
            keywords = []
            for key in relevant:
                if isinstance(key, str):
                    keywords.append(key)
            assert len(keywords) >= 1, "List of relevant keywords is empty!"
            return (relevant[keywords[0]].charge_state, relevant_candidates)

        if self.parms["verbose"] is True:
            print("Multiple relevant candidates meet all selection criteria")
        keywords = []
        for key in relevant:
            if isinstance(key, str):
                keywords.append(key)
        assert len(keywords) >= 1, "List of relevant keywords is empty!"
        charge_state = relevant[keywords[0]].charge_state
        for key, val in relevant.items():
            if val.charge_state == charge_state:
                continue
            if self.parms["verbose"] is True:
                print("WARNING::Multiple relevant candidates differ in charge_state!")
                print("WARNING::No unique solution possible for given criteria!")
            return (0, relevant_candidates)
        # not returned yet, so all relevant candidates have the same charge_state
        if self.parms["verbose"] is True:
            print(f"Multiple relevant candidates have all the same charge_state {charge_state}")
        if self.parms["sacrifice_isotopic_uniqueness"] is True:
            return (charge_state, relevant_candidates)
        if self.parms["verbose"] is True:
            print("WARNING::Multiple relevant candidates differ in isotopes!")
            print("WARNING::But these have the same charge_state!")
            print("WARNING::No unique solution possible for given criteria!")
        return (0, relevant_candidates)
