# Set of utility tools for parsing file formats used by atom probe.
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

from ifes_apt_tc_data_modeling.utils.definitions import \
    MAX_NUMBER_OF_ATOMS_PER_ION
from ifes_apt_tc_data_modeling.utils.utils import \
    create_isotope_vector, isotope_vector_to_nuclid_list, \
    isotope_vector_to_human_readable_name, \
    is_range_overlapping, is_range_significant

from ifes_apt_tc_data_modeling.utils.molecular_ions import MolecularIonCandidate

from ifes_apt_tc_data_modeling.nexus.nx_field import NxField


class NxIon():
    """Representative of a NeXus base class NXion."""

    def __init__(self, *args, **kwargs):
        self.comment = NxField("", "")  # comment, use e.g. for label of custom ion types
        self.color = NxField("", "")  # color used by software which created the dataset
        self.volume = NxField("", "")  # volume value in range files
        self.ion_type = NxField("", "")
        if len(args) >= 1:
            assert isinstance(args[0], list), "args[0] needs to be a list !"
            self.isotope_vector = NxField(create_isotope_vector(args[0]), "")
        elif "isotope_vector" in kwargs:
            assert isinstance(kwargs["isotope_vector"], np.ndarray), \
                "kwargs isotope_vector needs to be an np.ndarray !"
            assert np.shape(kwargs["isotope_vector"]) \
                == (MAX_NUMBER_OF_ATOMS_PER_ION,), \
                "kwargs isotope_vector needs be a " \
                + "(" + str(MAX_NUMBER_OF_ATOMS_PER_ION) + ",) array!"
            self.isotope_vector \
                = NxField(np.asarray(kwargs["isotope_vector"], np.uint16), "")
        else:
            # the default UNKNOWN IONTYPE
            self.isotope_vector = NxField(create_isotope_vector([]), "")
        self.nuclid_list = NxField(
            isotope_vector_to_nuclid_list(self.isotope_vector.typed_value), "")
        if "charge" in kwargs:
            assert isinstance(kwargs["charge"], int), \
                "kwargs charge needs to be an int !"
            assert kwargs["charge"] > -8, \
                "kwargs charge needs to be at least -7 !"
            assert kwargs["charge"] < +8, \
                "kwargs charge needs to be at most +7 !"
            self.charge = NxField(np.int8(kwargs["charge"]), "")
        else:
            # try to identify the charge state, will return NxField
            # with typed_value for charge on [0, +7]
            # here 0 flags and warns that it was impossible to recover
            # the relevant charge which is usually a sign that the range
            # is not matching the theoretically expect peak location
            self.charge = NxField(np.int8(0), "")
        self.name = NxField(isotope_vector_to_human_readable_name(
            self.isotope_vector.typed_value, self.charge.typed_value))
        self.ranges = NxField(np.empty((0, 2), np.float64), "amu")

    def add_range(self, mqmin: np.float64, mqmax: np.float64):
        """Adding mass-to-charge-state ratio interval."""
        assert is_range_significant(mqmin, mqmax) is True, \
            "Refusing to add epsilon range!"
        # the following example shows that is_range_overlapping should not be checked for
        # like it was in the past
        # ion.add_range(10.0, 12.0), ion.add_range(12.0, 13.3)
        # is equivalent to ion.add_range(10.0, 13.3)
        # assert is_range_overlapping(np.asarray([mqmin, mqmax]),
        #                            self.ranges.typed_value) is False, \
        #                            "Refusing overlapping range!"
        self.ranges.typed_value = np.vstack(
            (self.ranges.typed_value, np.array([mqmin, mqmax])))

    def update_human_readable_name(self):
        """Reevaluate charge and isotope_vector for name."""
        self.name = NxField(isotope_vector_to_human_readable_name(
            self.isotope_vector.typed_value, self.charge.typed_value))

    def report(self):
        """Report values."""
        print("ion_type")
        print(self.ion_type.typed_value)
        print("isotope_vector")
        print(self.isotope_vector.typed_value)
        print("nuclid_list")
        print(self.nuclid_list.typed_value)
        print("human-readable name")
        print(self.name.typed_value)
        print("charge")
        print(self.charge.typed_value)
        print("ranges")
        print(self.ranges.typed_value)
        print("comment")
        print(self.comment.typed_value)
        print("color")
        print(self.color.typed_value)
        print("volume")
        print(self.volume.typed_value)

    def add_charge_model(self,
                         parameters={},
                         candidates=[]):
        """Add details about the model how self.charge was defined."""
        self.charge_model = {}
        assert "min_abundance" in parameters.keys(), \
            "Parameter min_abundance not defined!"
        assert "min_abundance_product" in parameters.keys(), \
            "Parameter min_abundance_product not defined!"
        assert "min_half_life" in parameters.keys(), \
            "Parameter min_half_life not defined!"
        assert "sacrifice_isotopic_uniqueness" in parameters.keys(), \
            "Parameter sacrifice_isotopic_uniqueness not defined!"
        self.charge_model = {
            "isotope_matrix": [],
            "charge_vector": [],
            "mass_vector": [],
            "nat_abun_prod_vector": [],
            "min_half_life_vector": []}
        for key, val in parameters.items():
            if key not in self.charge_model.keys():
                self.charge_model[key] = val
        n_cand = len(candidates)
        # print("n_cand " + str(n_cand))
        if n_cand > 0:
            self.charge_model["isotope_matrix"] = np.zeros(
                (n_cand, MAX_NUMBER_OF_ATOMS_PER_ION), np.uint16)
            self.charge_model["charge_vector"] = np.zeros(
                (n_cand, ), np.int8)
            self.charge_model["mass_vector"] = np.zeros(
                (n_cand, ), np.float64)
            self.charge_model["nat_abun_prod_vector"] = np.zeros(
                (n_cand, ), np.float64)
            self.charge_model["min_half_life_vector"] = np.zeros(
                (n_cand, ), np.float64)
            row = 0
            for cand in candidates:
                if isinstance(cand, MolecularIonCandidate):
                    self.charge_model["isotope_matrix"][row, 0:len(cand.isotope_vector)] \
                        = cand.isotope_vector
                    self.charge_model["charge_vector"][row] = cand.charge
                    self.charge_model["mass_vector"][row] = cand.mass
                    self.charge_model["nat_abun_prod_vector"][row] \
                        = cand.abundance_product
                    self.charge_model["min_half_life_vector"][row] \
                        = cand.shortest_half_life
                    row += 1
                else:
                    print(__name__ + " found cand which is not a MolecularIonCandidate!")
        # else:
        #     print("Not enough candidates to report as a charge model")
