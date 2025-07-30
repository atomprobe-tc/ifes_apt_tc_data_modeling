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

"""Utility for parsing files via memory mapping."""

import typing
import mmap
import numpy as np


@typing.no_type_check
def get_memory_mapped_data(fpath: str, dtyp: str, oset: int, strd: int, shp: int):
    """Read typed data from memory-mapped file from offset with stride."""
    # https://stackoverflow.com/questions/60493766/ \
    #       read-binary-flatfile-and-skip-bytes for I/O access details

    with (
        open(fpath, "rb") as fp,
        mmap.mmap(fp.fileno(), length=0, access=mmap.ACCESS_READ) as memory_mapped,
    ):
        return np.ndarray(
            buffer=memory_mapped, dtype=dtyp, offset=oset, strides=strd, shape=shp
        ).copy()
    return None
