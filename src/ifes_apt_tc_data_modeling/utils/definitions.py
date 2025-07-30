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

"""Generic definitions when working with molecular ions useful for atom probe."""

import numpy as np

# default value for isotope_to_hash when the element is meant irrespective of its isotopes
# see also comment in ifes_apt_tc_data_modeling.utils.utils isotope_to_hash on a breaking
# change coming with the NIAC standardized NXapm uses 255 for marking the element while
# before the default value was 0
NEUTRON_NUMBER_FOR_ELEMENT = 255
# was 0 before breaking change, now required 255 as per NIAC standard, see NXatom nuclid_hash docstring
# restrict the number distinguished ion types
MAX_NUMBER_OF_ION_SPECIES = 256
# restrict number of atoms for molecular ion fragments
MAX_NUMBER_OF_ATOMS_PER_ION = 32
# practical and required minimum mass-resolution Da or atomic mass unit (amu)
MQ_EPSILON = np.float64(1.0 / 2000.0)


# do not consider isotopes with a very low natural abundance
PRACTICAL_ABUNDANCE = 1.0e-6  # 0.  # 1.0e-6
# do not consider candidate isotopically different molecular ions
# if their natural abundance product is too low
PRACTICAL_ABUNDANCE_PRODUCT = 0.0  # 1.0e-6  # 0.  # 1.0e-12
# do consider too shortliving isotopes
PRACTICAL_MIN_HALF_LIFE = np.inf
# many examples of ranges are not constrained strongly enough so that
# there are many isotopes (many of which admittedly hypothetical) ones
# which are within the range, this option lifts the constraint that
# there should be only one set of isotopically different molecular ions
# and if these have all these same charge it is assumed this is the
# charge of the molecular ion
# strictly speaking however one would have to rule out every possible
# molecular ion candidate (which is non-trivial and maybe not even
# with first principles theory possible...
SACRIFICE_ISOTOPIC_UNIQUENESS = True
