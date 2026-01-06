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
# https://doi.org/10.1063/1.2709758

import os
from typing import Any

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
        txt_stripped = [line for line in txt.split("\n") if line.strip() != ""]
        del txt

        instrument: dict[str, str] = {}
        event_data: dict[int, dict[str, list[float]]] = {}
        # first-level key is hit_group: int
        # second-level keys "tof", "det_x", "det_y", all: str
        voltage_keywords = [
            "standing_voltage",  # ureg.volt
            "pulse_voltage",  # ureg.volt
            "beta",  # coupling coefficient, ureg.dimensionless
        ]
        voltage_data: dict[int, dict[str, str]] = {}

        # need to parse sequence of pairs of a dash line that ought to be followed
        # by exactly one S line holding event data, ones such a pair is found it is
        # stored in the event_cache and only if the pair is valid that cache content
        # gets added to event_data otherwise the cache is cleared
        # lines store metadata and raw data of different type, the first character
        # identifies what content a line stores
        # prettify values after these have been parsed
        last_line = OPS_LINE_OTHER_OR_NONE
        voltage_sequence = 0
        hit_sequence = 0
        # hit_group_to_pulse_number_look_up: dict[int, int] = {}
        for line in txt_stripped:
            if line.startswith("*") or line == "":
                continue

            parts = []
            for value in line.strip().split():
                if value != "":
                    parts.append(value)

            if parts[0] == "C":
                if len(parts) == 5:
                    for idx, parameter in enumerate(
                        [
                            "flight_path",  # ureg.millimeter
                            "alpha",  # ureg.dimensionless
                            "beta",  # ureg.dimensionless
                            "t_zero",  # ureg.nanosecond
                        ]
                    ):
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
                if parts[0] == "IR":
                    instrument["reflectron"] = "yes"
                if len(parts) == 2:
                    if "detector_radius" not in instrument:
                        instrument["detector_radius"] = parts[1]
                        last_line = OPS_LINE_OTHER_OR_NONE
                    else:
                        logger.warning("OPS_READER_FORMAT_DUPLICATE_I")
                        return
                else:
                    logger.warning("OPS_READER_FORMAT_I")
                    return
            elif parts[0] == "P":
                if len(parts) == 2:
                    if "detector_channels" not in instrument:
                        instrument["detector_channels"] = parts[1]
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
                    voltage_sequence += 1
                    voltage_data[voltage_sequence] = {
                        "standing_voltage": parts[1],
                        "pulse_voltage": parts[2],
                    }
                    if len(parts) == 4:
                        voltage_data[voltage_sequence]["beta"] = parts[3]
                    else:
                        if "beta" in instrument:
                            voltage_data[voltage_sequence]["beta"] = instrument["beta"]
                        else:
                            logger.warning("OPS_READER_FORMAT_V_NO_BETA")
                            return
                    voltage_data[voltage_sequence]["next_hit_group_offset"] = (
                        f"{len(event_data.keys())}"
                    )
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

                # TODO::parsing of pulse_delta does not get dereferenced
                # pulse_delta = parts[0].replace("-", "")
                hit_sequence += 1
                event_data[hit_sequence] = {"tof": [], "det_x": [], "det_y": []}
                for idx in range(1, len(parts)):
                    try:
                        event_data[hit_sequence]["tof"].append(float(parts[idx]))
                    except ValueError:
                        logger.warning("OPS_READER_FORMAT_DASH_TOF_VALUE_ERROR")
                        return

                last_line = OPS_LINE_DASH
            elif parts[0].startswith("S"):
                if last_line != OPS_LINE_DASH:
                    if strict_mode:
                        logger.warning("OPS_READER_FORMAT_SLINE_PREFIX_ERR")
                        return
                    else:
                        continue  # hope for the best with the next line

                # TODO handle less positions than TOF, pop
                if ((len(parts) - 1) / (number_delay_lines + 1)) < len(
                    event_data[hit_sequence]["tof"]
                ):
                    if hit_sequence in event_data:
                        del event_data[hit_sequence]
                    last_line = OPS_LINE_OTHER_OR_NONE
                    continue

                try:
                    number_of_events = int(parts[0].replace("S", ""))
                except ValueError:
                    logger.warning("OPS_READER_FORMAT_S_NUMBER_OF_EVENTS")
                    if strict_mode:
                        return
                    if hit_sequence in event_data:
                        del event_data[hit_sequence]
                    last_line = OPS_LINE_OTHER_OR_NONE
                    continue

                if number_of_events != (len(parts) - 1) / (number_delay_lines + 1):
                    # -1 to ignore parts[0], +1 for two position values + tof map index
                    logger.warning("OPS_READER_FORMAT_S_EVENT_COUNT")
                    return

                if number_of_events < len(event_data[hit_sequence]["tof"]):
                    logger.warning("OPS_READER_FORMAT_S_EVENT_COUNT")
                    if strict_mode:
                        return
                    else:
                        if hit_sequence in event_data:
                            del event_data[hit_sequence]
                        last_line = OPS_LINE_OTHER_OR_NONE

                tof_values: list[float] = [0.0] * number_of_events
                det_x: list[float] = []
                det_y: list[float] = []
                # det_y = ["n/a"] * number_of_events
                for idx in range(0, len(event_data[hit_sequence]["tof"])):
                    tof_values[idx] = event_data[hit_sequence]["tof"][idx]

                healthy = True
                timing_index = number_of_events * 2 + 1
                for idx in range(1, timing_index, 2):
                    # when idx is 1, 3, 5, ... like here
                    # then offset = idx /2 is nothing but 0, 1, 2, ...
                    # so we can append
                    try:
                        det_x.append(float(parts[idx]))
                    except ValueError:
                        if strict_mode:
                            logger.warning("OPS_READER_FORMAT_S_DET_X_VALUE_ERROR")
                            return
                        else:
                            healthy = False
                            break
                    try:
                        det_y.append(float(parts[idx + 1]))
                    except ValueError:
                        if strict_mode:
                            logger.warning("OPS_READER_FORMAT_S_DET_Y_VALUE_ERROR")
                            return
                        else:
                            healthy = False
                            break
                if not healthy:
                    if hit_sequence in event_data:
                        del event_data[hit_sequence]
                        last_line = OPS_LINE_OTHER_OR_NONE

                time_map: list[int] = []
                for idx in range(timing_index, len(parts)):
                    try:
                        time_map.append(int(parts[idx]))
                    except ValueError:
                        if strict_mode:
                            logger.warning(f"OPS_READER_FORMAT_S_MAP_ENTRY {line}")
                            return
                        else:
                            healthy = False
                            break
                if not healthy:
                    if hit_sequence in event_data:
                        del event_data[hit_sequence]
                        last_line = OPS_LINE_OTHER_OR_NONE

                if len(time_map) != number_of_events:
                    logger.warning(f"OPS_READER_FORMAT_S_ASSERT")
                    return

                event_data[hit_sequence]["tof"] = [0.0] * number_of_events
                for idx in range(0, number_of_events):
                    if time_map[idx] != 0:
                        event_data[hit_sequence]["tof"][idx] = tof_values[time_map[idx] - 1]
                event_data[hit_sequence]["det_x"] = det_x.copy()
                event_data[hit_sequence]["det_y"] = det_y.copy()
                last_line = OPS_LINE_OTHER_OR_NONE
            else:
                logger.warning("OPS_READER_FORMAT line with an unknown format")
                return

        if last_line == OPS_LINE_DASH:
            if strict_mode:
                logger.warning("OPS_READER_FORMAT_TRAILING_DASH")
                return
            if hit_sequence in event_data:
                del event_data[hit_sequence]

        return
        # TODO clean data


# from datetime import datetime
# from zoneinfo import ZoneInfo
# local_tz = ZoneInfo("America/New_York")
# s = "Tue Sep  2 15:06:14 2003"
# dt = datetime.strptime(s, "%a %b %d %H:%M:%S %Y")
# dt_local = dt.replace(tzinfo=local_tz)
# dt_utc = dt_local.astimezone(ZoneInfo("UTC"))
# iso_utc = dt_utc.isoformat()
