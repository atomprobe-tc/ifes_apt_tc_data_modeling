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

"""Set of utility tools for parsing file formats used by atom probe."""


class NxField():
    """Representative of a NeXus field."""

    def __init__(self, values=None, unit: str = ""):
        self.parent = None
        self.is_a = None  # ontology reference concept ID e.g.
        self.values = values
        self.unit_category = None
        self.unit = unit
        self.attributes = None

    def get_value(self):
        """Get value."""
        return self.values

    def get_unit(self):
        """Get unit."""
        return self.unit
