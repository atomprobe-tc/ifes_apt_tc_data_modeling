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

"""PoSAP ops file format reader used by atom probe microscopists."""

# pylint: disable=duplicate-code

# ops is the acquisition format for raw data of the legacy Oxford Position-sensitive Atom Probe
# PoSAP instrument that paved the way for the modern LEAP type systems
# this is a Python implementation of the legacy PoSAP ops reader from libatom-probe
# https://sourceforge.net/p/apttools/libatomprobe/ci/default/tree/src/io/dataFiles.cpp#l1172
# https://www.repository.cam.ac.uk/items/2af37a7a-65d3-421f-ab06-a0ae2400b5f7 for example data

import os

import numpy as np

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data
from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg


class ReadOpsFileFormat:
    """Read *.ops file format."""

    def __init__(self, file_path: str, verbose: bool = False):
        """Initialize the reader."""
        self.supported = False
        if not file_path.lower().endswith(".ops"):
            logger.warning(f"{file_path} is likely not a PoSAP ops file")
            return
        self.supported = True
        self.file_path = file_path
        self.verbose = verbose

    def parse(self):
        """Interpret ops file."""
        pass
