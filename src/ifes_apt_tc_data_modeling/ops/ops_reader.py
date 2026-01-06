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

OPS_LINE_DASH = 0
OPS_LINE_OTHER_OR_NONE = 1


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

    def parse(self, number_delay_lines: int = 2, strict_mode: bool = True):
        """Interpret ops file."""
        with open(self.file_path, encoding="utf8") as ops_fp:
            txt = ops_fp.read()
        txt = txt.replace("\r\n", "\n")  # windows to unix EOL conversion
        txt = txt.replace(",", ".")  # use decimal dots instead of comma
        txt_stripped = [line
            for line in txt.split("\n")
            if line.strip() != ""
        ]
        del txt

        instrument: dict = {}
        event_data_keywords = ["voltage", "pulse_number", "time_of_flight"]  # det_x", "det_y"]
        event_data: dict[str, list] = {}  # collects all valid
        event_cache: dict[str, list] = {}
        voltage_keywords = ["voltage", "pulse_voltage", "beta"]
        voltage_cache: dict[str, float] = {}
        for keyword in event_data_keywords:
            event_data[keyword] = []
            event_cache[keyword] = []

        # need to parse sequence of pairs of a dash line that ought to be followed
        # by exactly one S line holding event data, ones such a pair is found it is
        # stored in the event_cache and only if the pair is valid that cache content
        # gets added to event_data otherwise the cache is cleared
        # lines store metadata and raw data of different type, the first character
        # identifies what content a line stores
        last_line = OPS_LINE_OTHER_OR_NONE
        hit_group = 0
        for line in txt_stripped:
            if line.startswith("*") or line == "":
                continue

            parts = [value.strip() for value in line.strip().split()]
            if parts[0] == "C":
                if len(parts) == 5:
                    for idx, parameter in enumerate(["flight_path", "alpha", "beta", "t_zero"]):
                        if parameter not in instrument:
                            instrument[parameter] = parts[idx + 1]
                            last_line = OPS_LINE_OTHER_OR_NONE
                        else:
                            logger.warning(f"OPS_READER_FORMAT_DUPLICATE_C {parameter}")
                            return
                else:
                    logger.warning("OPS_READER_FORMAT_C")
                    return
            elif parts[0] == "I" or parts[0] == "IR":
                if len(parts) == 2:
                    if "detector_radius" not in instrument:
                        try:
                            instrument["detector_radius"] = float(parts[1])
                            last_line = OPS_LINE_OTHER_OR_NONE
                        except ValueError:
                            logger.warning(f"OPS_READER_FORMAT_I value error {parts[1]}")
                            return
                    else:
                        logger.warning("OPS_READER_FORMAT_DUPLICATE_I")
                        return
                else:
                    logger.warning("OPS_READER_FORMAT_I")
                    return
            elif parts[0] == "P":
                if len(parts) == 2:
                    if "detector_channels" not in instrument:
                        instrument["detector_channels"] = int(parts[1])
                        last_line = OPS_LINE_OTHER_OR_NONE
                    else:
                        logger.warning("OPS_READER_FORMAT_DUPLICATE_P")
                        return
                else:
                    logger.warning("OPS_READER_FORMAT_P")
                    return
            elif parts[0] == "T":
                logger.warning("OPS_READER_FORMAT_T found but ignored")
                break
            elif parts[0] == "V":
                if 3 <= len(parts) <= 4:
                    event_cache["voltage"] = {}
                    event_cache["voltage"]["voltage"] = int(parts[1])
                    event_cache["voltage"]["pulse_voltage"] = int(parts[2])
                    if len(parts) == 4:
                        event_cache["voltage"]["beta"] = float(parts[3])
                    else:
                        if "beta" in instrument:
                            event_cache["voltage"]["beta"] = float(instrument["beta"])
                        else:
                            logger.warning("OPS_READER_FORMAT_V_NO_BETA")
                            return
                    event_cache["voltage"]["next_hit_group_offset"] = len(event_data[])  # WHICH_ONE ??
                    # ?? .append(voltage_data)
                    last_line = OPS_LINE_OTHER_OR_NONE
                else:
                    logger.warning("OPS_READER_FORMAT_V")
                    return
            elif parts[0].startswith("-"):
                if last_line == OPS_LINE_DASH:
                    logger.warning("OPS_READER_FORMAT_DOUBLE_DASH")
                    return
                if len(parts) < 2:
                    logger.warning("OPS_READER_FORMAT_DASH")
                    return

                pulse_delta = float(parts[0].replace("-", ""))
                # TODO what to do with pulse_number
                event_cache["pulse_number"].append(pulse_delta + 1)

                for idx in range(1, len(parts)):
                    event_cache["time_of_flight"].append(int(parts[idx]))

                last_line = OPS_LINE_DASH
            elif parts[0].startswith("S"):
                if last_line != OPS_LINE_DASH:
                    if strict_mode:
                        logger.warning("OPS_READER_FORMAT_SLINE_PREFIX_ERR")
                        return
                    else:
                        continue  # hope for the best with the next line

                # TODO handle less positions than TOF, pop

                number_of_events = int(parts[0].replace("S", ""))
                # TODO handle number_of_events pop

                if number_of_events != (len(parts) - 1) / (number_delay_lines + 1):
                    logger.warning("OPS_READER_FORMAT_S_EVENT_COUNT")
                    return

                if number_of_events < len(event_cache["time_of_flight"]):
                    if strict_mode:
                        logger.warning("OPS_READER_FORMAT_S_EVENT_COUNT")
                        return
                    else:
                        # TODO event_count pop
                        last_line = OPS_LINE_OTHER_OR_NONE

                # TODO
            else:
                logger.warning("OPS_READER_FORMAT line with an unknown format")
                return

        if last_line == OPS_LINE_DASH:
            if strict_mode:
                logger.warning("OPS_READER_FORMAT_TRAILING_DASH_ERR")
            discard = event_data.pop()

        return
        # TODO clean data