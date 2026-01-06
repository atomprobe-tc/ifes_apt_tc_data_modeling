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

"""AMETEK APT(6) data exchange file reader used by atom probe microscopists."""

import numpy as np

APT_SECTION_NAME_MAX_LENGTH = 32
APT_SECTION_TYPE_MAX_LENGTH = 32


def np_uint16_to_string(uint16_array: np.ndarray) -> str:
    """Create string from array of uint16 numbers (UTF-16)."""
    str_parsed = ""
    for value in uint16_array:
        if value != 0:  # '\x00'
            str_parsed += chr(value)
    return f"{str_parsed}"


def string_to_typed_nparray(string: str, length: int, data_type: type) -> np.ndarray:
    """Create length long specifically typed numpy array from string."""
    if isinstance(data_type, type) and (len(string) <= length):
        numpy_array = np.zeros(length, dtype=data_type)  # type: ignore
        for value in np.arange(0, len(string)):
            numpy_array[value] = ord(string[value])
        return numpy_array
    raise ValueError(
        f"{data_type} is either not a type or {string} is not <= {length}."
    )
