# Utility tool to work with molecular ions in atom probe microscopy.
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

import typing

from typing import Tuple

import numpy as np

import radioactivedecay as rd

import ase
from ase.data import atomic_numbers, atomic_masses, atomic_names

from ifes_apt_tc_data_modeling.utils.utils import \
    isotope_to_hash, hash_to_isotope, isotope_vector_to_dict_keyword

from ifes_apt_tc_data_modeling.utils.nist_isotope_data import isotopes

# do not consider isotopes with a very low natural abundance
PRACTICAL_ABUNDANCE = 0.  # 1.0e-6
# do not consider candidate isotopically different molecular ions
# if their natural abundance product is too low
PRACTICAL_ABUNDANCE_PRODUCT = 0.  # 1.0e-12
# do consider too shortliving isotopes
PRACTICAL_MIN_HALF_LIFE = np.inf
# many examples of ranges are not constrainted strongly enough so that
# there are many isotopes (many of which admittedly hypothetical) ones
# which are within the range, this option lifts the constraint that
# there should be only one set of isotopically different molecular ions
# and if these have all these same charge it is assumed this is the
# charge of the molecular ion
# strictly speaking however one would have to rule out every possible
# molecular ion candidate (which is non-trivial and maybe not even
# with first principles theory possible...
SACRIFICE_ISOTOPIC_UNIQUENESS = True
VERBOSE=False


class MolecularIonCandidate:
    def __init__(self, ivec=[], charge=0, mass_sum=0., nat_abun_prod=0., min_half_life=np.inf):
        self.isotope_vector = np.asarray(np.atleast_1d(ivec), np.uint16)
        self.charge = np.int8(charge)
        self.mass = np.float64(mass_sum)
        self.abundance_product = np.float64(nat_abun_prod)
        self.shortest_half_life = np.float64(min_half_life)

    def unique_keyword(self):
        keyword = ""
        keyword += isotope_vector_to_dict_keyword(
            np.sort(np.asarray(self.isotope_vector, np.uint16), kind="stable")[::-1])
        keyword += "__" + str(self.charge)
        return keyword


class MolecularIonBuilder:
    def __init__(self,
                 min_abundance=1.0e-6,
                 min_abundance_product=1.0e-6,
                 min_half_life=np.inf,
                 sacrifice_uniqueness=True,
                 verbose=False):
        self.nuclids = np.asarray([], np.uint16)
        self.element_isotopes = {}
        self.nuclid_mass = {}
        self.nuclid_abundance = {}
        self.nuclid_stable = {}  # observationally stable
        self.nuclid_unclear = {}  # unclear halflife
        self.nuclid_halflife = {}
        self.candidates = []
        self.parms = {}
        self.parms["min_abundance"] = min_abundance
        self.parms["min_abundance_product"] = min_abundance_product
        self.parms["min_half_life"] = min_half_life
        self.parms["sacrifice_isotopic_uniqueness"] = sacrifice_uniqueness
        self.parms["verbose"] = verbose

        for symbol, atomic_number in atomic_numbers.items():
            if symbol != "X":
                # assume that data from ase take preference
                # if half-life data are available in rd take these if not mark
                # the isotope that it unclear_half_life == True
                element_isotopes = []
                for mass_number in isotopes[atomic_number]:
                    half_life = np.inf
                    observationally_stable = False
                    unclear_half_life = False

                    # test if half-life data available
                    trial_nuclid_name = symbol + "-" + str(mass_number)
                    try:
                        tmp = rd.Nuclide(trial_nuclid_name)
                    except:
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
                        half_life = np.nan
                        unclear_half_life = True

                    # get ase abundance data
                    n_protons = atomic_numbers[symbol]
                    n_neutrons = mass_number - n_protons
                    mass = isotopes[atomic_numbers[symbol]][mass_number]["mass"]
                    abundance = isotopes[atomic_numbers[symbol]][mass_number]["composition"]
                    # self.nuclids.append(
                    #     Isotope(mass,
                    #             abundance,
                    #             n_protons, n_neutrons,
                    #             observationally_stable,
                    #             unclear_half_life))
                    hashvalue = isotope_to_hash(int(n_protons), int(n_neutrons))
                    self.nuclids = np.asarray(np.append(self.nuclids, hashvalue), np.uint16)
                    self.nuclid_mass[hashvalue] = np.float64(mass)
                    self.nuclid_abundance[hashvalue] = np.float64(abundance)
                    self.nuclid_stable[hashvalue] = observationally_stable
                    self.nuclid_unclear[hashvalue] = unclear_half_life
                    self.nuclid_halflife[hashvalue] = half_life
                    assert np.uint16(0) not in self.nuclid_mass.keys(), \
                        "0 must not be a key in nuclid_mass dict!"
                    element_isotopes = np.append(element_isotopes, hashvalue)
                self.element_isotopes[atomic_number] = np.sort(
                    np.asarray(element_isotopes, np.uint16), kind="stable")[::-1]
        self.nuclids = np.sort(self.nuclids, kind="stable")[::-1]
        if self.parms["verbose"] is True:
            print(f"MolecularIonBuilder initialized with {len(self.nuclids)}")

    def get_element_isotopes(self, hashvalue):
        """List of hashvalues all isotopes of element specified by hashvalue."""
        # nuclids = np.asarray([], np.uint16)
        target_p, target_n = hash_to_isotope(hashvalue)
        return self.element_isotopes[target_p]
        # for source_hashvalue in self.nuclids:
        #     # naive search for now
        #     source_p, source_n = hash_to_isotope(source_hashvalue)
        #     if source_p != target_p:
        #         continue
        #     else:
        #         nuclids = np.append(nuclids, source_hashvalue)
        # return nuclids

    def get_isotope_mass_sum(self, nuclid_arr=[]):
        """Evaluate cumulated atomic_mass of isotopes in ivec."""
        # assuming no relativistic effects or other quantum effects
        # mass loss due to charge state considered insignificant
        mass = 0.
        for hashvalue in nuclid_arr:
            if hashvalue != 0:
                mass += self.nuclid_mass[hashvalue]
            # else:  # because ivec are systematically filled sorted in descending order !
            #   break
        return mass

    def get_natural_abundance_product(self, nuclid_arr=[]):
        abun_prod = 1.
        for hashvalue in nuclid_arr:
            if hashvalue != 0:
                abun_prod *= self.nuclid_abundance[hashvalue]
        return abun_prod

    def get_shortest_half_life(self, nuclid_arr=[]):
        min_half_life = self.parms["min_half_life"]
        for hashvalue in nuclid_arr:
            if hashvalue != 0:
                if self.nuclid_halflife[hashvalue] != np.nan:
                    min_half_life = np.min((min_half_life, self.nuclid_halflife[hashvalue]))
                    # if the min_half_life is np.inf than we know that every nuclid in the
                    # molecular ion is observationally stable
                    # if not we know that considering this molecular ion might not be a good
                    # idea or only necessary in very few cases because even if we were
                    # to find such a combination of nuclids at least one would decay in observable
                    # time and this might possibly hint that studying this molecular is tricky
                else:
                    # if we do not information about the half-life it is very likely that
                    # this is an exotic nuclids likely never found in the wild
                    return np.nan
        return min_half_life

    def combinatorics(self, element_arr, low, high):
        # RNG/RRNG range files do store element information for each range
        # BUT not isotope information, correspondingly this can yield
        # only an isotope_vector whose hashvalues have ALL in common
        # that the number of neutrons is 0, i.e. their hashvalue is the atomic_number
        # element_arr is an isotope_vector/ivec with such hashvalues
        # as an example the molecular ion C:2 H:1 will map to 6, 6, 1
        max_depth = 0  # number of non-zero entries in element_arr
        for hashvalue in element_arr:
            n_protons, n_neutrons = hash_to_isotope(hashvalue)
            assert n_neutrons == 0, \
                str(hashvalue) + " in element_arr_in has n_neutrons != 0 !"
            if hashvalue != 0:
                max_depth += 1
            # else:
            #     break
        if self.parms["verbose"] is True:
            print(f"Maximum recursion depth {max_depth}")
        self.candidates = []
        if self.parms["verbose"] is True:
            print(element_arr)

        if max_depth > 0:
            depth = 0
            ith_nuclids = self.get_element_isotopes(element_arr[depth])
            cand_arr_curr = []  # combinatorially add nuclids while recursing deeper
            self.iterate_molecular_ion(
                element_arr, ith_nuclids, cand_arr_curr,
                depth, max_depth, low, high)
            if self.parms["verbose"] is True:
                print(f"Found {len(self.candidates)} candidates!")
            return self.try_to_reduce_to_unique_solution()
            # will return a tuple of charge and list of relevant_candidates
        else:
            return (0, [])

    def iterate_molecular_ion(self,
                              element_arr, jth_nuclids, cand_arr_prev,
                              i, max_n, low, high):
        if i < (max_n - 1):
            for nuclid in jth_nuclids:

                ixxth_nuclids = self.get_element_isotopes(element_arr[i + 1])
                cand_arr_curr = np.append(cand_arr_prev, nuclid)
                self.iterate_molecular_ion(
                    element_arr, ixxth_nuclids, cand_arr_curr,
                    i + 1, max_n, low, high)
        elif i == (max_n - 1):
            for nuclid in jth_nuclids:
                cand_arr_curr = np.append(cand_arr_prev, nuclid)
                # by this design the ivec does not necessarily remain ordered

                new_mass = self.get_isotope_mass_sum(cand_arr_curr)
                new_abun_prod = self.get_natural_abundance_product(cand_arr_curr)
                new_halflife = self.get_shortest_half_life(cand_arr_curr)

                for new_chrg in np.arange(1, 8):
                    mass_to_charge = new_mass / new_chrg;
                    if mass_to_charge < low:
                        break
                    # we can break the entire charge state generation here already instead of continue
                    # because already the current mq is left out of interval [mqmin, mqmax]
                    # and all mq in the next iterations will result in even lower mass-to-charge
                    # you walk to the left increasing your distance to the left bound
                    if mass_to_charge > high:
                        # can be optimized and broken out of earlier if testing first chrg == 1
                        # and then chrg == APTMOLECULARION_MAX_CHANGE
                        continue;
                        # must not be break here because with adding more charge
                        # we usually walk from right to left eventually into [low, high] !
                    # molecular ion is within user-specified bounds
                    self.candidates.append(
                        MolecularIonCandidate(cand_arr_curr, new_chrg, new_mass, new_abun_prod, new_halflife))
                # do not start a new recursion
        else:
            # break the recursion
            return

    def try_to_reduce_to_unique_solution(self):
        """Hueristics to identify if current candidates are unique."""
        if self.parms["verbose"] is True:
            print(f"Reduce set of {len(self.candidates)} candidates to a unique...")
        self.relevant = {}
        for cand in self.candidates:
            if cand.abundance_product >= self.parms["min_abundance_product"]:
                if np.isnan(cand.shortest_half_life) == False:
                    if cand.shortest_half_life >= self.parms["min_half_life"]:
                        keyword = cand.unique_keyword()
                        if keyword not in self.relevant.keys():
                            self.relevant[keyword] = cand
        relevant_candidates = []
        if self.parms["verbose"] is True:
            print(f"Reduced set to {len(self.relevant.keys())} relevant candidates...")
            for key, val in self.relevant.items():
                print(key)
        for key, obj in self.relevant.items():
            relevant_candidates.append(obj)
        # print(type(relevant_candidates))

        if len(self.relevant) == 0:
            if self.parms["verbose"] is True:
                print("WARNING::No relevant candidate meets all criteria!")
                print("WARNING::No solution possible for given criteria!")
            return (0, relevant_candidates)
        elif len(self.relevant) == 1:
            if self.parms["verbose"] is True:
                print("One relevant candidate which meets all criteria")
            keywords = []
            for key in self.relevant.keys():
                if isinstance(key, str):
                    keywords.append(key)
            assert len(keywords) >= 1, "List of relevant keywords is empty!"
            return (self.relevant[keywords[0]].charge, relevant_candidates)
        else:
            if self.parms["verbose"] is True:
                print("Multiple relevant candidates meet all selection criteria")
            keywords = []
            for key in self.relevant.keys():
                if isinstance(key, str):
                    keywords.append(key)
            assert len(keywords) >= 1, "List of relevant keywords is empty!"
            charge = self.relevant[keywords[0]].charge
            for key, val in self.relevant.items():
                if val.charge == charge:
                    continue
                else:
                    if self.parms["verbose"] is True:
                        print("WARNING::Multiple relevant candidates differ in charge!")
                        print("WARNING::No unique solution possible for given criteria!")
                    return (0, relevant_candidates)
            # not returned yet, so all relevant candidates have the same charge
            if self.parms["verbose"] is True:
                print(f"Multiple relevant candidates have all the same charge {charge}")
            if self.parms["sacrifice_isotopic_uniqueness"] is True:
                return (charge, relevant_candidates)

            if self.parms["verbose"] is True:
                print("WARNING::Multiple relevant candidates differ in isotopes!")
                print("WARNING::But these have the same charge!")
                print("WARNING::No unique solution possible for given criteria!")
            return (0, relevant_candidates)
