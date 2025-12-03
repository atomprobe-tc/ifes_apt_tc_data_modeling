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

import re
from ifes_apt_tc_data_modeling.utils.molecular_ions import get_chemical_symbols


def get_all_elements_with_given_number_of_characters(
    number_characters: int = 0,
) -> list:
    selection = []
    for symbol in get_chemical_symbols():
        if len(symbol) == number_characters:
            selection.append(symbol)
    return selection


def parse_elements(element: str) -> list[str]:
    parsed = []
    elements_with_one_character_symbol = (
        get_all_elements_with_given_number_of_characters(1)
    )
    elements_with_two_character_symbol = (
        get_all_elements_with_given_number_of_characters(2)
    )
    idx = 0
    while idx < len(element):
        if f"{element[idx : idx + 2]}" in elements_with_two_character_symbol:
            match = re.search(r"^[0-9]+", element[idx + 2 : :])
            if match:
                parsed += [element[idx : idx + 2]] * int(match.group(0))
                idx += 1 + len(f"{match.group(0)}")
            else:
                parsed.append(f"{element[idx : idx + 2]}")
                idx += 2
        elif f"{element[idx]}" in elements_with_one_character_symbol:
            match = re.search(r"^[0-9]+", element[idx + 1 : :])
            if match:
                parsed += [f"{element[idx]}"] * int(match.group(0))
                idx += 1 + len(f"{match.group(0)}")
            else:
                if idx + 2 < len(element):
                    if (
                        f"{element[idx : idx + 2]}"
                        in elements_with_two_character_symbol
                    ):
                        match = re.search(r"^[0-9]+", element[idx + 2 : :])
                        if match:
                            parsed += [f"{element[idx]}"] * int(match.group(0))
                            idx += 2 + len(f"{match.group(0)}")
                        else:
                            parsed.append(f"{element[idx : idx + 2]}")
                            idx += 2
                    else:
                        parsed.append(f"{element[idx]}")
                        idx += 1
                else:
                    parsed.append(f"{element[idx]}")
                    idx += 1
        else:
            idx += 1
    return parsed
