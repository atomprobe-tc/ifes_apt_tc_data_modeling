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


@pytest.mark.parametrize(
    "nuclide_hash, left, right",  # , expected",
    [
        # (create_nuclide_hash(["H"]), 0.0, 10.0),
        # (create_nuclide_hash(["Be"]), 5.0, 17.0),
        (create_nuclide_hash(["Tc"]), 84.0, 120.0),
        # (create_nuclide_hash(["Yb"]), 165.0, 206.0),
    ],
    ids=[
        # "hydrogen_landscape",
        # "beryllium_landscape",
        "technetium_landscape",
        # "ytterbium_landscape",
    ],
)
@pytest.mark.skip(
    reason=f"TODO needs reporting and comparison against common cases!"
    f"Failing on Tc expected cuz constrained to non-radioactive only!"
)
def test_combinatorial_analysis(
    nuclide_hash: np.uint16, left: np.float64, right: np.float64
):  # expected):
    m_ion = NxIon(nuclide_hash=nuclide_hash, charge_state=1)
    m_ion.add_range(left, right)
    m_ion.comment = NxField("", "")
    m_ion.apply_combinatorics()
    m_ion.report()
    print(m_ion.charge_state_model["nuclide_hash"])
    # print(m_ion.charge_state_model["shortest_half_life"])
    assert True


"""
class TestClass:
    def test_one(self):
        x = "this"
        assert "h" in x

    def test_two(self):
        x = "hello"
        assert hasattr(x, "check")
"""
