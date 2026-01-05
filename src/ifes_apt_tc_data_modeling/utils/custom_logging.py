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

import logging

# https://docs.python.org/3/howto/logging.html

try:
    from pynxtools_apm.utils.custom_logging import logger
except ImportError as exc:
    DEFAULT_LOGGER_NAME = "ifes_apt_tc_data_modeling"
    logger = logging.getLogger(DEFAULT_LOGGER_NAME)
    logging.basicConfig(
        filename=f"{DEFAULT_LOGGER_NAME}.log",
        filemode="w",  # use "a" to collect all in a session, use "w" to overwrite
        format="%(levelname)s %(asctime)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
        encoding="utf-8",
        level=logging.DEBUG,
    )
