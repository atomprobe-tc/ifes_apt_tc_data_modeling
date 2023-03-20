# Utility for parsing files via memory mapping.
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

import mmap

import numpy as np


@typing.no_type_check
def get_memory_mapped_data(file_name: str, data_type: str, oset: int,
                           strd: int, shp: int):
    """Memory-maps file plus offset strided read of typed data."""
    # https://stackoverflow.com/questions/60493766/ \
    #       read-binary-flatfile-and-skip-bytes for I/O access details

    with open(file_name, "rb") as file_handle, \
            mmap.mmap(file_handle.fileno(), length=0, access=mmap.ACCESS_READ) as memory_mapped:
        return np.ndarray(buffer=memory_mapped, dtype=data_type,
                          offset=oset, strides=strd, shape=shp).copy()
    return None
