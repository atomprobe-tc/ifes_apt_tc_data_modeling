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

import pytest
import numpy as np
from ifes_apt_tc_data_modeling.nexus.nx_ion import NxIon, NxField
from ifes_apt_tc_data_modeling.utils.utils import create_nuclide_hash

# commented out expected test results are those returned when PRACTICAL_MINIMUM_HALF_LIFE = 0.0
@pytest.mark.parametrize(
    "nuclide_hash, left, right, expected",
    [
        (create_nuclide_hash(["H"]), 0.0, 10.0, "H"),
        (create_nuclide_hash(["Be"]), 5.0, 17.0, "Be +"),  # Be
        (create_nuclide_hash(["Tc"]), 84.0, 120.0, "Tc"), # Tc +
        (create_nuclide_hash(["Ra"]), 216.0, 236.0, "Ra"),
        (create_nuclide_hash(["U"]), 228.0, 249.0, "U"),
        (create_nuclide_hash(["Th"]), 222.0, 242.0, "Th"),
        (create_nuclide_hash(["Yb"]), 165.0, 206.0, "Yb +"),
        (create_nuclide_hash(["Fr"]), 222.0, 224.0, "Fr"),  # "Fr +"
        (create_nuclide_hash(["Cr", "Cr", "O"]), 57.819, 61.159, "Cr Cr O ++"),
    ],
    ids=[
        "hydrogen_landscape",
        "beryllium_landscape",
        "technetium_landscape",
        "radon_landscape",
        "uranium_landscape",
        "thorium_landscape",
        "ytterbium_landscape",
        "francium_landscape",
        "molecular_ion_cr_cr_o",
    ],
)
def test_combinatorial_analysis(
    nuclide_hash: np.uint16, left: np.float64, right: np.float64, expected: str):
    m_ion = NxIon(nuclide_hash=nuclide_hash, charge_state=1)
    m_ion.add_range(left, right)
    m_ion.comment = NxField("", "")
    m_ion.apply_combinatorics()
    m_ion.report()
    # print(m_ion.charge_state_model["nuclide_hash"])
    assert expected == m_ion.name.values


"""
class TestClass:
    def test_one(self):
        x = "this"
        assert "h" in x

    def test_two(self):
        x = "hello"
        assert hasattr(x, "check")
"""
