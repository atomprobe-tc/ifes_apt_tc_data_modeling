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

"""Function implementing combinatorial analysis for all ranging definitions."""

# shared functionality used for RNG/RRNG/ENV

import numpy as np

from ifes_apt_tc_data_modeling.nexus.nx_ion import NxField, NxIon
from ifes_apt_tc_data_modeling.utils.molecular_ions import \
    MolecularIonBuilder, PRACTICAL_ABUNDANCE, PRACTICAL_ABUNDANCE_PRODUCT, \
    PRACTICAL_MIN_HALF_LIFE, VERBOSE, SACRIFICE_ISOTOPIC_UNIQUENESS


def apply_combinatorics(m_ion: NxIon):
    """Apply specifically constrainted combinatorial analysis."""
    crawler = MolecularIonBuilder(min_abundance=PRACTICAL_ABUNDANCE,
                                  min_abundance_product=PRACTICAL_ABUNDANCE_PRODUCT,
                                  min_half_life=PRACTICAL_MIN_HALF_LIFE,
                                  sacrifice_uniqueness=SACRIFICE_ISOTOPIC_UNIQUENESS,
                                  verbose=VERBOSE)
    recovered_charge_state, m_ion_candidates \
        = crawler.combinatorics(m_ion.isotope_vector.values,
                                m_ion.ranges.values[0, 0],
                                m_ion.ranges.values[0, 1])
    # print(f"{recovered_charge_state}")
    m_ion.charge_state = NxField(np.int8(recovered_charge_state), "")
    m_ion.update_human_readable_name()
    m_ion.add_charge_state_model({"min_abundance": PRACTICAL_ABUNDANCE,
                                  "min_abundance_product": PRACTICAL_ABUNDANCE_PRODUCT,
                                  "min_half_life": PRACTICAL_MIN_HALF_LIFE,
                                  "sacrifice_isotopic_uniqueness": SACRIFICE_ISOTOPIC_UNIQUENESS},
                                 m_ion_candidates)
    return m_ion
