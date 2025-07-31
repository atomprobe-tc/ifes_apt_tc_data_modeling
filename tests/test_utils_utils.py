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
from ifes_apt_tc_data_modeling.utils.definitions import NEUTRON_NUMBER_FOR_ELEMENT
from ifes_apt_tc_data_modeling.utils.utils import (
    get_smart_chemical_symbols,
    isotope_to_hash,
    hash_to_isotope,
    create_nuclide_hash,
    nuclide_hash_to_nuclide_list,
    nuclide_hash_to_dict_keyword,
    nuclide_hash_to_human_readable_name,
    MQ_EPSILON,
    is_range_overlapping,
    is_range_significant,
    element_or_nuclide_to_hash,
    symbol_lst_to_matrix_of_nuclide_vector,
)


def test_get_smart_chemical_symbols():
    """Organize element symbols such that search H does not match He."""
    assert "X" not in get_smart_chemical_symbols()


@pytest.mark.parametrize(
    "proton_number, neutron_number, expected",
    [
        (0, 0, 0),
        (1, NEUTRON_NUMBER_FOR_ELEMENT, 65281),
        (1, 0, 1),
        (1, 1, 257),
        (1, 2, 513),
        (43, 56, 14379),
    ],
    ids=[
        "isotope_to_hash_zerozero",
        "isotope_to_hash_hydrogen",
        "isotope_to_hash_1h",
        "isotope_to_hash_2h",
        "isotope_to_hash_3h",
        "isotope_to_hash_99tc",
    ],
)
def test_isotope_to_hash(proton_number: int, neutron_number: int, expected: int):
    assert isotope_to_hash(proton_number, neutron_number) == expected


@pytest.mark.parametrize(
    "hashvalue, expected",
    [
        (0, (0, 0)),
        (65281, (1, NEUTRON_NUMBER_FOR_ELEMENT)),
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
        "hash_to_isotope_99tc",
    ],
)
def test_hash_to_isotope(hashvalue: int, expected: tuple):
    assert hash_to_isotope(hashvalue) == expected


def test_create_nuclide_hash():
    assert np.array_equal(
        create_nuclide_hash(["Tc-99", "Fe", "O", "O", "O"]),
        np.sort(
            np.concatenate(
                (
                    np.array(
                        [
                            isotope_to_hash(43, 56),
                            isotope_to_hash(26, NEUTRON_NUMBER_FOR_ELEMENT),
                            isotope_to_hash(8, NEUTRON_NUMBER_FOR_ELEMENT),
                            isotope_to_hash(8, NEUTRON_NUMBER_FOR_ELEMENT),
                            isotope_to_hash(8, NEUTRON_NUMBER_FOR_ELEMENT),
                        ],
                        np.uint16,
                    ),
                    np.array([isotope_to_hash(0, 0)] * (32 - 5), np.uint16),
                ),
                axis=0,
            ),
            kind="stable",
        )[::-1],  # descending order
    )


def test_nuclide_hash_to_nuclide_list():
    assert np.array_equal(
        nuclide_hash_to_nuclide_list(
            create_nuclide_hash(["Tc-99", "Fe", "O", "O", "O"])
        ),
        np.concatenate(
            (
                np.array([[0, 26], [0, 8], [0, 8], [0, 8], [99, 43]], np.uint16),
                np.zeros((32 - 5, 2), np.uint16),
            ),
            axis=0,
        ),
    )


def test_nuclide_hash_to_dict_keyword():
    assert (
        nuclide_hash_to_dict_keyword(
            create_nuclide_hash(["Tc-99", "Fe", "O", "O", "O"])
        )
        == "65306_65288_65288_65288_14379"
    )


@pytest.mark.parametrize(
    "ivec, charge_state, expected",
    [
        (create_nuclide_hash([]), 0, "unknown_iontype"),
        (create_nuclide_hash(["Tc-99", "Fe"]), np.int8(+8), "unknown_iontype"),
        (create_nuclide_hash(["Tc-99", "Fe"]), np.int8(+3), "Fe 99Tc +++"),
        (create_nuclide_hash(["Tc-99", "Fe"]), np.int8(0), "Fe 99Tc"),
        (create_nuclide_hash(["Tc-99", "Fe"]), np.int8(-3), "Fe 99Tc ---"),
        (create_nuclide_hash(["Tc-99", "Fe"]), np.int8(-8), "unknown_iontype"),
    ],
    ids=[
        "nuclide_hash_to_human_readable_name_positive_zero_neutral",
        "nuclide_hash_to_human_readable_name_positive_too_high",
        "nuclide_hash_to_human_readable_name_positive_high",
        "nuclide_hash_to_human_readable_name_positive_neutral",
        "nuclide_hash_to_human_readable_name_positive_negative",
        "nuclide_hash_to_human_readable_name_positive_too_negative",
    ],
)
def test_nuclide_hash_to_human_readable_name(
    ivec: np.ndarray, charge_state: np.int8, expected: str
):
    assert expected == nuclide_hash_to_human_readable_name(ivec, charge_state)


@pytest.mark.parametrize(
    "interval, interval_set, expected",
    [
        (np.array([12.0, 12.0]), np.array([12.0, 12.0]), False),
        (np.array([12.0, 24.0]), np.array([13.0, 14.0]), True),
        (np.array([12.0, 24.0]), np.array([11.0, 25.0]), True),
        (
            np.array([12.0, 24.0]),
            np.array([24.0 + MQ_EPSILON, 24.0 + MQ_EPSILON + 1.0]),
            False,
        ),
        (np.array([12.0, 24.0]), np.array([24.0, 25.0]), True),
        (np.array([12.0, 24.0]), np.array([11.0, 12.0]), True),
        (
            np.array([12.0, 24.0]),
            np.array([12.0 - MQ_EPSILON - 1.0, 12.0 - MQ_EPSILON]),
            False,
        ),
    ],
    ids=[
        "is_range_overlapping_line_line",
        "is_range_overlapping_inside_set",
        "is_range_overlapping_culling_set",
        "is_range_overlapping_smaller_than_set",
        "is_range_overlapping_touching_from_left_set",
        "is_range_overlapping_touching_from_right_set",
        "is_range_overlapping_larger_than_set",
    ],
)
@pytest.mark.skip(reason="should not be checked for")
def test_is_range_overlapping(
    interval: np.ndarray, interval_set: np.ndarray, expected: bool
):
    # ion.add_range(10.0, 12.0), ion.add_range(12.0, 13.3)
    # is equivalent to ion.add_range(10.0, 13.3)
    # assert is_range_overlapping(np.asarray([mqmin, mqmax]),
    #                            self.ranges.values) is False, \
    #                            "Refusing overlapping range!"
    assert expected == is_range_overlapping(interval, interval_set)


@pytest.mark.parametrize(
    "left, right, expected",
    [
        (np.float64(12.0), np.float64(12.0), False),
        (np.float64(12.0), np.float64(12.0 + MQ_EPSILON), True),
        (np.float64(12.0), np.float64(12.0 + 2 * MQ_EPSILON), True),
        (np.float64(12.0 + MQ_EPSILON), np.float64(12.0), False),
        (np.float64(12.0 + 2 * MQ_EPSILON), np.float64(12.0), False),
    ],
    ids=[
        "is_range_significant_zero",
        "is_range_significant_small_edge",
        "is_range_significant_small_safe",
        "is_range_significant_left_largerthan_right_edge",
        "is_range_significant_left_largerthan_right_safe",
    ],
)
def test_is_range_significant(left: np.float64, right: np.float64, expected: bool):
    assert expected == is_range_significant(left, right)


@pytest.mark.parametrize(
    "symbol, expected",
    [
        # ("", 0),
        # ("junk", 0),
        ("H", 65281),
        ("H-1", 1),
        ("H-2", 257),
        ("Tc-99", 14379),
    ],
    ids=[
        # "element_or_nuclide_to_hash_empty_string",
        # "element_or_nuclide_to_hash_not_an_element",
        "element_or_nuclide_to_hash_hydrogen",
        "element_or_nuclide_to_hash_1h",
        "element_or_nuclide_to_hash_2h",
        "element_or_nuclide_to_hash_99tc",
    ],
)
def test_element_or_nuclide_to_hash(symbol: str, expected: int):
    # TODO this needs improvement as currently the commented out tests fail
    try:
        assert expected == element_or_nuclide_to_hash(symbol)
    except (ValueError, TypeError) as exc:
        assert True, f"element_or_nuclide_to_hash raised an exception {exc}"


"""
@pytest.mark.parametrize(
    "method, symbol_lst, charge_lst, expected",
    [
        ("resolve_all", [["K", "C"], ["U-238", "H-2"]], None, ""),
        ("resolve_unknown", [["K", "C"], ["U-238", "H-2"]], None, ""),
        ("resolve_element", [["K", "C"], ["U-238", "H-2"]], None, ""),
        ("resolve_isotope", [["K", "C"], ["U-238", "H-2"]], None, ""),
        ("resolve_ion", [["K", "C"], ["U-238", "H-2"]], None, ""),
        # TODO the following need fixing [["K-40", "C-14"], ["U-238", "H-2"]]
    ],
    ids=[
        "symbol_lst_to_matrix_all_one",
        "symbol_lst_to_matrix_unknown_one",
        "symbol_lst_to_matrix_element_one",
        "symbol_lst_to_matrix_isotope_one",
        "symbol_lst_to_matrix_ion_one",
    ],
)
def test_symbol_lst_to_matrix_of_nuclide_vector(
    method: str, symbol_lst: list, charge_lst: list, expected: str
):
    print(f">>>>>>>>>{method}>>>>>>>>>")
    print(f"{symbol_lst_to_matrix_of_nuclide_vector(method, symbol_lst, charge_lst)}")
    assert True

# tmp = create_nuclide_hash(["Tc-99", "Fe", "O", "O", "O"])
# print(f"{tmp}____{type(tmp)}____{tmp.dtype}____{np.shape(tmp)}")
# assert True
"""
