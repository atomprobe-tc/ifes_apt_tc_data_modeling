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

"""Set of utility tools for parsing file formats used by atom probe."""

# pylint: disable=too-many-instance-attributes

import numpy as np

from ifes_apt_tc_data_modeling.utils.definitions import \
    MAX_NUMBER_OF_ATOMS_PER_ION
from ifes_apt_tc_data_modeling.utils.utils import \
    create_isotope_vector, isotope_vector_to_nuclid_list, \
    isotope_vector_to_human_readable_name, is_range_significant
from ifes_apt_tc_data_modeling.utils.molecular_ions import MolecularIonCandidate
from ifes_apt_tc_data_modeling.nexus.nx_field import NxField


class NxIon():
    """Representative of a NeXus base class NXion."""

    def __init__(self, *args, **kwargs):
        self.comment = NxField("", "")  # comment, use e.g. for label of custom ion types
        self.color = NxField("", "")  # color used by software which created the dataset
        self.volume = NxField("", "")  # volume value in range files
        self.ion_type = NxField("", "")
        self.charge_state_model = {}
        if len(args) >= 1:
            if isinstance(args[0], list) is False:
                raise ValueError("args[0] needs to be a list!")
            self.isotope_vector = NxField(create_isotope_vector(args[0]), "")
        elif "isotope_vector" in kwargs:
            if isinstance(kwargs["isotope_vector"], np.ndarray) is False:
                raise ValueError("kwargs isotope_vector needs to be an np.ndarray!")
            if np.shape(kwargs["isotope_vector"]) != (MAX_NUMBER_OF_ATOMS_PER_ION,):
                raise ValueError(
                    f"kwargs isotope_vector needs be a ({MAX_NUMBER_OF_ATOMS_PER_ION},) array!")
            self.isotope_vector = NxField(np.asarray(kwargs["isotope_vector"], np.uint16), "")
        else:
            # the default UNKNOWN IONTYPE
            self.isotope_vector = NxField(create_isotope_vector([]), "")
        self.nuclid_list = NxField(isotope_vector_to_nuclid_list(self.isotope_vector.values), "")
        if "charge_state" in kwargs:
            if isinstance(kwargs["charge_state"], int) \
                    and (-8 < kwargs["charge_state"] < +8):
                self.charge_state = NxField(np.int8(kwargs["charge_state"]), "")
        else:
            # charge_state 0 flags and warns that it was impossible to recover
            # the relevant charge which is usually a sign that the range
            # is not matching the theoretically expect peak location
            self.charge_state = NxField(np.int8(0), "")
        self.name = NxField(isotope_vector_to_human_readable_name(
            self.isotope_vector.values, self.charge_state.values))
        self.ranges = NxField(np.empty((0, 2), np.float64), "amu")

    def add_range(self, mqmin: np.float64, mqmax: np.float64):
        """Adding mass-to-charge-state ratio interval."""
        if is_range_significant(mqmin, mqmax) is False:
            raise ValueError(f"Refusing to add epsilon range [{mqmin}, {mqmax}] !")
        # the following example shows that is_range_overlapping should not be checked for
        # like it was in the past
        # ion.add_range(10.0, 12.0), ion.add_range(12.0, 13.3)
        # is equivalent to ion.add_range(10.0, 13.3)
        # assert is_range_overlapping(np.asarray([mqmin, mqmax]),
        #                            self.ranges.values) is False, \
        #                            "Refusing overlapping range!"
        self.ranges.values = np.vstack((self.ranges.values, np.array([mqmin, mqmax])))

    def update_human_readable_name(self):
        """Re-evaluate charge and isotope_vector for name."""
        self.name = NxField(isotope_vector_to_human_readable_name(
            self.isotope_vector.values, self.charge_state.values))

    def report(self):
        """Report values."""
        print(f"ion_type: {self.ion_type.values}\n"
              f"isotope_vector: {self.isotope_vector.values}\n"
              f"nuclid_list: {self.nuclid_list.values}\n"
              f"human-readable name: {self.name.values}\n"
              f"charge_state: {self.charge_state.values}\n"
              f"ranges: {self.ranges.values}\n"
              f"comment: {self.comment.values}\n"
              f"color: {self.color.values}\n"
              f"volume: {self.volume.values}\n")

    def add_charge_state_model(self,
                               parameters,
                               candidates):
        """Add details about the model how self.charge_state was defined."""
        self.charge_state_model = {}
        req_parms = ["min_abundance", "min_abundance_product",
                     "min_half_life", "sacrifice_isotopic_uniqueness"]
        for req in req_parms:
            if req in parameters:
                continue
            raise ValueError(f"Parameter {req} not defined in parameters dict!")
        self.charge_state_model = {"isotope_matrix": [],
                                   "charge_state_vector": [],
                                   "mass_vector": [],
                                   "nat_abun_prod_vector": [],
                                   "min_half_life_vector": []}
        for key, val in parameters.items():
            if key not in self.charge_state_model:
                self.charge_state_model[key] = val
        n_cand = len(candidates)
        if n_cand > 0:
            self.charge_state_model["isotope_matrix"] \
                = np.zeros((n_cand, MAX_NUMBER_OF_ATOMS_PER_ION), np.uint16)
            self.charge_state_model["charge_state_vector"] \
                = np.zeros((n_cand, ), np.int8)
            self.charge_state_model["mass_vector"] \
                = np.zeros((n_cand, ), np.float64)
            self.charge_state_model["nat_abun_prod_vector"] \
                = np.zeros((n_cand, ), np.float64)
            self.charge_state_model["min_half_life_vector"] \
                = np.zeros((n_cand, ), np.float64)
            row_idx = 0
            for cand in candidates:
                if isinstance(cand, MolecularIonCandidate):
                    self.charge_state_model["isotope_matrix"][row_idx, 0:len(cand.isotope_vector)] \
                        = cand.isotope_vector
                    self.charge_state_model["charge_state_vector"][row_idx] = cand.charge_state
                    self.charge_state_model["mass_vector"][row_idx] = cand.mass
                    self.charge_state_model["nat_abun_prod_vector"][row_idx] \
                        = cand.abundance_product
                    self.charge_state_model["min_half_life_vector"][row_idx] \
                        = cand.shortest_half_life
                    row_idx += 1
                else:
                    print(f"{__name__} found cand which is not a MolecularIonCandidate!")
        # else:
        #     print("Not enough candidates to report as a charge_state model")
