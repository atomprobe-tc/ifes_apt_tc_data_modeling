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

from ifes_apt_tc_data_modeling.utils.string_handling import rchop


@pytest.mark.parametrize(
    "string, suffix, expected",
    [("nochop", "", "nochop"), ("nochop", "chop", "no"), ("chop_me", "_me", "chop")],
    ids=["rchop-nochop", "rchop-no-underscore", "rchop-underscore"],
)
def test_string_handling_rchop(string: str, suffix: str, expected: str):
    assert rchop(string, suffix) == expected
