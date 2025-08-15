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

# default value for isotope_to_hash when the element is encoded irrespective of its isotopes
# see also comment in ifes_apt_tc_data_modeling.utils.utils isotope_to_hash on a breaking
# change coming with the NIAC standardized NXapm uses 255 for marking the element while
# before the default value was 0
NEUTRON_NUMBER_FOR_ELEMENT = 255
# restrict the number distinguished ion types
MAX_NUMBER_OF_ION_SPECIES = 256
# restrict number of atoms for molecular ion fragments
MAX_NUMBER_OF_ATOMS_PER_ION = 32
# practical and required minimum mass-resolution Da or atomic mass unit (amu)
MQ_EPSILON = np.float64(1.0 / 2000.0)

# three types of nuclids are distinguished based on their half life:
# i) half life infinite modeled as np.inf, stable never practically observed decaying,
# ii) half life finite, known half life, in many cases irrelevant except for geo
# iii) half life unclear, exotic nuclids modeled as np.nan, typically irrelevant for atom probe

# do not consider isotopes with a low natural abundance
PRACTICAL_ABUNDANCE = 0.0
# do not consider molecular ions whose product of natural abundance low
PRACTICAL_ABUNDANCE_PRODUCT = 0.0
# do consider too shortliving isotopes
# PRACTICAL_MIN_HALF_LIFE = 0.0  # this setting would give access to the most exhaustive search
# but in reality even a specimen with radioactive nuclides needs to be measureable which takes some
# time during which the fast decaying nuclides might have already experienced that many iterations that
# no such atom of the nuclid remains, let's assume that this practical half life is one minute, then
# even if one rushes setting up the experiment several dozen iterations will pass until
# that isotope is measured
PRACTICAL_MIN_HALF_LIFE = np.inf  # 60.0
# if false, accept only solutions whose charge state and combination of isotopes is unique
# if true, relax this constraint
SACRIFICE_ISOTOPIC_UNIQUENESS = True
