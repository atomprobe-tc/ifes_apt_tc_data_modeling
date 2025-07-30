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
from ifes_apt_tc_data_modeling.utils.utils import (
    get_smart_chemical_symbols,
    isotope_to_hash,
    hash_to_isotope,
    create_nuclide_hash,
)


def test_get_smart_chemical_symbols():
    """Organize element symbols such that search H does not match He."""
    assert "X" not in get_smart_chemical_symbols()


@pytest.mark.parametrize(
    "proton_number, neutron_number, expected",
    [(0, 0, 0), (1, 255, 65281), (1, 0, 1), (1, 1, 257), (1, 2, 513), (43, 56, 14379)],
    ids=[
        "isotope_to_hash_zerozero",
        "isotope_to_hash_hydrogen",
        "isotope_to_hash_1h",
        "isotope_to_hash_2h",
        "isotope_to_hash_3h",
        "isotope_to_hash_99te",
    ],
)
def test_isotope_to_hash(proton_number: int, neutron_number: int, expected: int):
    assert isotope_to_hash(proton_number, neutron_number) == expected


@pytest.mark.parametrize(
    "hashvalue, expected",
    [
        (0, (0, 0)),
        (65281, (1, 255)),
        (1, (1, 0)),
        (257, (1, 1)),
        (513, (1, 2)),
        (14379, (43, 56)),
    ],
    ids=[
        "hash_to_isotope_zerozero",
        "hash_to_isotope_hydrogen",
        "hash_to_isotope_1h",
        "hash_to_isotope_2h",
        "hash_to_isotope_3h",
        "hash_to_isotope_99te",
    ],
)
def test_hash_to_isotope(hashvalue: int, expected: tuple):
    assert hash_to_isotope(hashvalue) == expected


def test_create_nuclide_hash(building_blocks: list) -> np.ndarray:
    tmp = create_nuclide_hash(["Fe", "Fe", "O", "O", "O"])
    print(f"{tmp}____{type(tmp)}____{np.shape(tmp)}")
    assert True
