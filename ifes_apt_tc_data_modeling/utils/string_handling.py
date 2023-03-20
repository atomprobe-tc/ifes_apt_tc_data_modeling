# Utility string mangling when parsing data in atom probe microscopy.
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

# pylint: disable=E1101

import typing

from typing import Tuple


def rchop(string: str = "", suffix: str = "") -> str:
    """Right-chop a string."""
    if suffix and string.endswith(suffix):
        return string[:-len(suffix)]
    return string
