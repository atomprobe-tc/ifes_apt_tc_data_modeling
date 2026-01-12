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

import re
from datetime import datetime
from zoneinfo import ZoneInfo

import numpy as np

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
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
        self.parse()

    def parse(
        self,
        number_delay_lines: int = 2,
        strict_mode: bool = True,
        local_time_zone=ZoneInfo("Europe/London"),  # as PoSAP was located in Oxford
    ):
        """Interpret ops file."""
        with open(self.file_path, encoding="utf8") as ops_fp:
            txt = ops_fp.read()
        txt = txt.replace("\r\n", "\n")  # windows to unix EOL conversion
        txt = txt.replace(",", ".")  # use decimal dots instead of comma
        lines = [line.strip() for line in txt.split("\n")]
        del txt

        self.instrument: dict[str, str] = {}
        self.event_data: dict[int, dict[str, list[float]]] = {}
        # first-level key is hit_group: int
        # second-level keys "tof", "det_x", "det_y", all: str
        # voltage_keywords = [
        #     "standing_voltage",  # ureg.volt
        #     "pulse_voltage",  # ureg.volt
        #     "beta",  # coupling coefficient, ureg.dimensionless
        # ]
        self.voltage_data: dict[int, dict[str, str]] = {}
        # first-level key is hit_group: int
        # second-level key standing_voltage, pulse_voltage, beta, all: str
        self.pulse_number: dict[int, int] = {}
        # first-level key is hit_group: int

        # need to parse sequence of pairs of a dash line that ought to be followed
        # by exactly one S line holding event data, ones such a pair is found it is
        # stored in the event_cache and only if the pair is valid that cache content
        # gets added to event_data otherwise the cache is cleared
        # lines store metadata and raw data of different type, the first character
        # identifies what content a line stores
        # prettify values after these have been parsed
        last_line = OPS_LINE_OTHER_OR_NONE
        self.voltage_sequence: int = 0
        self.hit_sequence: int = 0
        self.last_pulse_number: int = 0
        # hit_group_to_pulse_number_look_up: dict[int, int] = {}

        # check if first header line gives us date information
        pattern = re.compile(
            r"""
            ^
            \*
            \s+
            (?:Mon|Tue|Wed|Thu|Fri|Sat|Sun)
            \s+
            (?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)
            \s+
            \d{1,2}
            \s+
            \d{2}:\d{2}:\d{2}
            \s+
            \d{4}
            """,
            re.VERBOSE,
        )
        match = re.match(pattern, lines[0])
        if match:
            dt = datetime.strptime(
                match.group(0).replace("*", "").strip(), "%a %b %d %H:%M:%S %Y"
            )
            self.instrument["time_stamp"] = (
                dt.replace(tzinfo=local_time_zone)
                .astimezone(ZoneInfo("UTC"))
                .isoformat()
            )
        for jdx, line in enumerate(lines):
            if line.startswith("*") or line == "":
                continue

            parts = []
            for value in line.split():
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
                        if parameter not in self.instrument:
                            self.instrument[parameter] = parts[idx + 1]
                        else:
                            logger.warning(f"OPS_READER_FORMAT_DUPLICATE_C {parameter}")
                            return
                else:
                    logger.warning("OPS_READER_FORMAT_C")
                    return
                last_line = OPS_LINE_OTHER_OR_NONE
            elif parts[0] == "I" or parts[0] == "IR":
                if parts[0] == "I":
                    self.instrument["reflectron"] = "no"
                else:
                    self.instrument["reflectron"] = "yes"

                if len(parts) == 2:
                    if "detector_radius" not in self.instrument:
                        self.instrument["detector_radius"] = parts[1]
                    else:
                        logger.warning("OPS_READER_FORMAT_DUPLICATE_I")
                        return
                else:
                    logger.warning("OPS_READER_FORMAT_I")
                    return
                last_line = OPS_LINE_OTHER_OR_NONE
            elif parts[0] == "P":
                if len(parts) == 2:
                    if "detector_channels" not in self.instrument:
                        self.instrument["detector_channels"] = parts[1]
                    else:
                        logger.warning("OPS_READER_FORMAT_DUPLICATE_P")
                        return
                else:
                    logger.warning("OPS_READER_FORMAT_P")
                    return
                last_line = OPS_LINE_OTHER_OR_NONE
            elif parts[0] == "T":
                logger.warning("OPS_READER_FORMAT_T found but ignored")
                last_line = OPS_LINE_OTHER_OR_NONE
            elif parts[0] == "V":
                if len(parts) != 3 and len(parts) != 4:
                    logger.warning("OPS_READER_FORMAT_V")
                    return
                else:
                    self.voltage_sequence += 1
                    self.voltage_data[self.voltage_sequence] = {
                        "standing_voltage": parts[1],
                        "pulse_voltage": parts[2],
                    }
                    if len(parts) == 4:
                        self.voltage_data[self.voltage_sequence]["beta"] = parts[3]
                    else:
                        if "beta" in self.instrument:
                            self.voltage_data[self.voltage_sequence]["beta"] = (
                                self.instrument["beta"]
                            )
                        else:
                            logger.warning("OPS_READER_FORMAT_V_NO_BETA")
                            return
                    self.voltage_data[self.voltage_sequence][
                        "next_hit_group_offset"
                    ] = f"{len(self.event_data.keys())}"
                    last_line = OPS_LINE_OTHER_OR_NONE
            elif parts[0].startswith("-"):
                if last_line == OPS_LINE_DASH:
                    logger.warning("OPS_READER_FORMAT_DOUBLE_DASH")
                    return
                if len(parts) < 2:
                    logger.warning("OPS_READER_FORMAT_DASH")
                    return

                self.hit_sequence += 1
                # pulse_delta are parsed by libatomprobe but thereafter
                try:
                    pulse_delta = int(parts[0].replace("-", ""))
                    self.last_pulse_number += pulse_delta + 1
                    self.pulse_number[self.hit_sequence] = self.last_pulse_number
                except ValueError:
                    logger.warning(f"OPS_READER_FORMAT_DASH pulse_delta {parts[0]}")
                    return

                self.event_data[self.hit_sequence] = {
                    "tof": [],
                    "det_x": [],
                    "det_y": [],
                }
                for idx in range(1, len(parts)):
                    try:
                        self.event_data[self.hit_sequence]["tof"].append(
                            float(parts[idx])
                        )
                    except ValueError:
                        logger.warning("OPS_READER_FORMAT_DASH_TOF_VALUE_ERROR")
                        return

                last_line = OPS_LINE_DASH
            elif parts[0].startswith("S"):
                if last_line != OPS_LINE_DASH:
                    logger.warning("OPS_READER_FORMAT_SLINE_PREFIX_ERR")
                    if strict_mode:
                        return
                    else:
                        last_line = OPS_LINE_OTHER_OR_NONE
                        continue  # hope for the best with the next line

                # TODO handle less positions than TOF, pop
                if ((len(parts) - 1) / (number_delay_lines + 1)) < len(
                    self.event_data[self.hit_sequence]["tof"]
                ):
                    if self.hit_sequence in self.event_data:
                        del self.event_data[self.hit_sequence]
                    last_line = OPS_LINE_OTHER_OR_NONE
                    continue

                try:
                    number_of_events = int(parts[0].replace("S", ""))
                except ValueError:
                    logger.warning("OPS_READER_FORMAT_S_NUMBER_OF_EVENTS")
                    if strict_mode:
                        return
                    if self.hit_sequence in self.event_data:
                        del self.event_data[self.hit_sequence]
                    last_line = OPS_LINE_OTHER_OR_NONE
                    continue

                if number_of_events != (len(parts) - 1) / (number_delay_lines + 1):
                    # -1 to ignore parts[0], +1 for two position values + tof map index
                    logger.warning("OPS_READER_FORMAT_S_EVENT_COUNT")
                    return

                if number_of_events < len(self.event_data[self.hit_sequence]["tof"]):
                    logger.warning("OPS_READER_FORMAT_S_EVENT_COUNT")
                    if strict_mode:
                        return
                    else:
                        if self.hit_sequence in self.event_data:
                            del self.event_data[self.hit_sequence]
                        last_line = OPS_LINE_OTHER_OR_NONE
                        continue

                tof_values: list[float] = [0.0] * number_of_events
                det_x: list[float] = []
                det_y: list[float] = []
                for idx in range(0, len(self.event_data[self.hit_sequence]["tof"])):
                    tof_values[idx] = self.event_data[self.hit_sequence]["tof"][idx]

                healthy = True
                timing_index = number_of_events * 2 + 1
                for idx in range(1, timing_index, 2):
                    # when idx is 1, 3, 5, ... like here
                    # then offset = idx /2 is nothing but 0, 1, 2, ...
                    # so we can append
                    try:
                        det_x.append(float(parts[idx]))
                    except ValueError:
                        logger.warning("OPS_READER_FORMAT_S_DET_X_VALUE_ERROR")
                        if strict_mode:
                            return
                        else:
                            healthy = False
                            break  # no need to continue on that loop go to healthy check
                    try:
                        det_y.append(float(parts[idx + 1]))
                    except ValueError:
                        logger.warning("OPS_READER_FORMAT_S_DET_Y_VALUE_ERROR")
                        if strict_mode:
                            return
                        else:
                            healthy = False
                            break
                if not healthy:
                    if self.hit_sequence in self.event_data:
                        del self.event_data[self.hit_sequence]
                    last_line = OPS_LINE_OTHER_OR_NONE
                    continue

                time_map: list[int] = []
                for idx in range(timing_index, len(parts)):
                    try:
                        time_map.append(int(parts[idx]))
                    except ValueError:
                        logger.warning(f"OPS_READER_FORMAT_S_MAP_ENTRY {line}")
                        if strict_mode:
                            return
                        else:
                            healthy = False
                            break
                if not healthy:
                    if self.hit_sequence in self.event_data:
                        del self.event_data[self.hit_sequence]
                    last_line = OPS_LINE_OTHER_OR_NONE
                    continue

                if len(time_map) != number_of_events:
                    logger.warning(
                        "OPS_READER_FORMAT_S_TIME_MAP_NUMBER_OF_EVENTS_ASSERT"
                    )
                    return

                self.event_data[self.hit_sequence]["tof"] = [0.0] * number_of_events
                for idx in range(0, number_of_events):
                    if time_map[idx] > 0:
                        # > 0 because for Python using int, while in C++ was unsigned int
                        # also > 0 required cuz index time_map[idx] - 1 needs to be >= 0
                        self.event_data[self.hit_sequence]["tof"][idx] = tof_values[
                            time_map[idx] - 1
                        ]
                    # no else statement as values were already set to zero
                self.event_data[self.hit_sequence]["det_x"] = det_x.copy()
                self.event_data[self.hit_sequence]["det_y"] = det_y.copy()
                last_line = OPS_LINE_OTHER_OR_NONE
            else:
                logger.warning("OPS_READER_FORMAT line with an unknown format")
                return

        if last_line == OPS_LINE_DASH:
            if strict_mode:
                logger.warning("OPS_READER_FORMAT_TRAILING_DASH")
                return
            if self.hit_sequence in self.event_data:
                del self.event_data[self.hit_sequence]

        print(f"Parsed {len(lines)}, now normalizing...")
        print(f"instrument, {self.instrument}")
        print(f"len(event_data.keys()), {len(self.event_data.keys())}")
        print(f"len(voltage_data.keys()), {len(self.voltage_data.keys())}")
        print(f"len(pulse_number.keys()), {len(self.pulse_number.keys())}")
        print(f"last_pulse_number, {self.last_pulse_number}")
        print(f"voltage_sequence, {self.voltage_sequence}")
        print(f"hit_sequence, {self.hit_sequence}")

        parameter_names = [
            ("tof", ureg.nanosecond),
            ("det_x", ureg.dimensionless),  # not 100% sure, but should be relative
            ("det_y", ureg.dimensionless),  # not 100% sure
        ]
        event_data_stats = {}
        for parameter_name, unit in parameter_names:
            event_data_stats[parameter_name] = 0
            for hit_sequence, event_dict in self.event_data.items():
                if parameter_name in event_dict:
                    event_data_stats[parameter_name] += len(event_dict[parameter_name])
            print(
                f"event_data_stats[{parameter_name}], {event_data_stats[parameter_name]}"
            )

        # build flattened event data into pint quantities
        self.events: dict[str, ureg.Quantity] = {}
        for parameter_name, unit in parameter_names:
            numpy_array = np.zeros((event_data_stats[parameter_name],), np.float32)
            idx = 0
            for hit_sequence, event_dict in self.event_data.items():
                if parameter_name in event_dict:
                    for jdx in range(0, len(event_dict[parameter_name])):
                        numpy_array[idx + jdx] = event_dict[parameter_name][jdx]
                    idx += len(event_dict[parameter_name])
            print(f"numpy_array.dtype {numpy_array.dtype}")
            print(f"np.shape(numpy_array) {np.shape(numpy_array)}")
            self.events[parameter_name] = ureg.Quantity(numpy_array, unit)
            del numpy_array
            print(
                f"event_data_flattened {parameter_name}, {self.events[parameter_name]}"
            )

        # type convert and build flattened voltage data into pint quantities
        self.voltages: dict[str, ureg.Quantity] = {}
        parameter_names = [
            ("standing_voltage", ureg.volt),
            ("pulse_voltage", ureg.volt),
            ("beta", ureg.dimensionless),
            ("next_hit_group_offset", ureg.dimensionless),
        ]
        voltage_data_stats = {}
        for parameter_name, unit in parameter_names:
            voltage_data_stats[parameter_name] = len(self.voltage_data.keys())
            if parameter_name != "next_hit_group_offset":
                numpy_array = np.zeros(
                    (voltage_data_stats[parameter_name],), np.float32
                )
            else:
                numpy_array = np.zeros((voltage_data_stats[parameter_name],), np.int64)
            idx = 0
            for voltage_sequence, voltage_dict in self.voltage_data.items():
                # will report in insertion order
                try:
                    if parameter_name != "next_hit_group_offset":
                        numpy_array[idx] = float(voltage_dict[parameter_name])
                    else:
                        numpy_array[idx] = int(voltage_dict[parameter_name])
                    idx += 1
                except ValueError:
                    logger.warning("ValueError during voltage data conversion")
                    return

            print(f"numpy_array.dtype {numpy_array.dtype}")
            print(f"np.shape(numpy_array) {np.shape(numpy_array)}")
            self.voltages[parameter_name] = ureg.Quantity(numpy_array, unit)
            del numpy_array
            print(
                f"voltage_data_flattened {parameter_name}, {self.voltages[parameter_name]}"
            )

        for parameter_name, unit in [
            ("flight_path", ureg.millimeter),
            ("alpha", ureg.dimensionless),
            ("beta", ureg.dimensionless),
            ("t_zero", ureg.nanosecond),
            ("detector_radius", ureg.millimeter),  # possibly incorrect?
            ("detector_channels", ureg.dimensionless),
        ]:
            if parameter_name in self.instrument:
                self.instrument[parameter_name] = ureg.Quantity(
                    self.instrument[parameter_name], unit
                )
        print(f"instrument, {self.instrument}")

        # finally only expose flattened event data
        del self.event_data
        del self.voltage_data
