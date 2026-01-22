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

"""AMETEK APT(6) data exchange file reader used by atom probe microscopists."""

import os

import numpy as np
import pandas as pd

from ifes_apt_tc_data_modeling.apt.apt6_headers import AptFileHeaderMetadata
from ifes_apt_tc_data_modeling.apt.apt6_sections import AptFileSectionMetadata
from ifes_apt_tc_data_modeling.apt.apt6_sections_branches import EXPECTED_SECTIONS
from ifes_apt_tc_data_modeling.apt.apt6_utils import np_uint16_to_string
from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data
from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg


class ReadAptFileFormat:
    """Read AMETEK's open exchange *.apt file format."""

    def __init__(self, file_path: str, verbose: bool = False):
        self.supported = False
        if not file_path.lower().endswith(".apt"):
            logger.warning(f"{file_path} is likely not a Cameca APT file")
            return
        self.supported = True
        self.file_path = file_path
        self.verbose = verbose
        self.file_size = os.path.getsize(self.file_path)
        logger.debug(f"Reading {self.file_path} which is {self.file_size} B")

        self.header_section = None
        self.byte_offsets: dict = {}
        self.available_sections: dict = {}
        # where do sections start bytes from beginning of the file?
        self.apt: dict = {}

        self.parse_file_structure()

    def parse_file_structure(self):
        """Parse APT file header plus flat collection of metadata/data pairs.

        Each pair has a so-called section header and a corresponding data block.
        Section headers detail the content of the immediately trailing data block.
        An APT file can store none, some, or all of the possible sections.
        Furthermore, the file can contain additional pieces of information
        which this parser currently cannot read-out because the APT format as
        not yet all details of the APT file format specification have been shared
        by AMETEK/Cameca. Indeed, the IVAS/APSuite source code is currently the only
        reliable source of information about which content the sections encode
        and how these get formatted when exporting an APT file from APSuite
        for a specific version and build number and type of experiment plus
        combinations of settings.
        Parse header of the file and check which parsable sections the file
        contains get the byte offsets of the sections from the beginning
        relative to the start/beginning of the file.
        """
        self.byte_offsets = {}
        self.header_section = None
        self.available_sections = {}

        with open(self.file_path, "rb") as fp:
            self.dummy_header = AptFileHeaderMetadata()
            found_header = np.fromfile(
                fp, self.dummy_header.get_numpy_struct(), count=1
            )

            if not self.dummy_header.matches(found_header):
                raise TypeError(
                    "Found an unexpectedly formatted header. Create an issue to help us fix this."
                )

            logger.info(f"File describes {found_header['llIonCount'][0]} ions")

            self.header_section = found_header
            self.byte_offsets["header"] = np.uint64(fp.tell())
            logger.debug(f"Currently at byte_offset {self.byte_offsets['header']} B")

            end_of_file_not_reached = b"yes"
            while end_of_file_not_reached != b"":
                # probe for end of file
                end_of_file_not_reached = fp.read(1)
                if end_of_file_not_reached != b"":
                    fp.seek(-1, os.SEEK_CUR)
                else:
                    logger.debug(f"End of file at {fp.tell()} B")
                    break

                dummy_section = AptFileSectionMetadata()
                found_section = np.fromfile(
                    fp, dummy_section.get_numpy_struct(), count=1
                )
                keyword = np_uint16_to_string(found_section["wcSectionType"][0])

                logger.debug(f"keyword: {keyword}, found_section: {found_section}")
                if keyword in self.available_sections:
                    raise ValueError("Found a duplicate of an already parsed section.")
                if keyword not in ["Delta Pulse", "Epos ToF"]:
                    if keyword not in EXPECTED_SECTIONS:
                        raise ValueError(
                            "Found an unknown section, seems like an unknown/new branch."
                        )
                    metadata_section = EXPECTED_SECTIONS[keyword]
                    if metadata_section.matches(found_section):
                        self.available_sections[keyword] = metadata_section
                else:
                    logger.warning(
                        f"Found uninterpretable non-registered section, create an issue to help us fix this, parsing continues at llByteCount {found_section['llByteCount'][0]} B"
                    )

                self.byte_offsets[keyword] = np.uint64(fp.tell())
                if keyword == "Position":
                    # special case six IEEE 32-bit floats preceding raw data
                    self.byte_offsets[keyword] += np.uint64(6 * 4)
                self.byte_offsets[keyword] += np.uint64(found_section["llByteCount"][0])
                logger.debug(
                    f"Byte offset for reading data for section: {keyword}"
                    f" {self.byte_offsets[keyword]} B"
                )
                fp.seek(self.byte_offsets[keyword], os.SEEK_SET)

    # one convenience reader function for every known section
    # is useful because it structures the parsers, enables reading the file
    # partially and reduces main memory consumption during full parsing
    def get_header(self):
        """Report metadata in the header."""
        metadata_dict = {
            "cSignature": np_uint16_to_string(self.header_section["cSignature"][0]),
            "iHeaderSize": np.int32(self.header_section["iHeaderSize"][0]),
            "iHeaderVersion": np.int32(self.header_section["iHeaderVersion"][0]),
            "wcFilename": np_uint16_to_string(self.header_section["wcFilename"][0]),
            "ftCreationTime": np.uint64(self.header_section["ftCreationTime"][0]),
            "llIonCount": np.uint64(self.header_section["llIonCount"][0]),
        }
        # check e.g. https://gist.github.com/Mostafa-Hamdy-Elgiar/
        # 9714475f1b3bc224ea063af81566d873 repo
        # for converting Windows/MSDN time to Python time
        for key, value in iter(metadata_dict.items()):
            logger.debug(f"{key}: {value}")

    def get_metadata(self, keyword: str):
        """Report available metadata for quantity if it exists."""
        if (keyword in self.available_sections) and (keyword in self.byte_offsets):
            metadata_dict = self.available_sections[keyword].get_metadata()
            for key, value in iter(metadata_dict.items()):
                logger.debug(f"{key}: {value}")

    def get_metadata_table(self):
        """Create table from all metadata for each section."""
        column_names = ["section"]  # header
        if "Mass" not in self.available_sections:
            raise ValueError(
                "Mass section not available to guide creation of the table header."
            )
        for key in self.available_sections["Mass"].get_metadata().keys():
            column_names.append(key)
        data_frame = pd.DataFrame(columns=column_names)

        for keyword, value in self.available_sections.items():
            row_dct = {"section": keyword}
            row_dct = {**row_dct, **value.get_metadata()}
            # logger.debug(value.get_metadata())
            row_df = pd.DataFrame(row_dct, index=[0])
            data_frame = pd.concat([data_frame, row_df], ignore_index=True)

        data_frame.style.format(precision=3, thousands=",", decimal=".").format_index(
            str.upper, axis=1
        )
        return data_frame

    def get_named_quantity(self, keyword: str):
        """Read quantity with name in keyword from APT file if it exists."""
        if (keyword in self.available_sections) and (keyword in self.byte_offsets):
            byte_position_start = (
                self.byte_offsets[keyword]
                - self.available_sections[keyword].get_ametek_size()
            )
            logger.info(f"Reading section {keyword} at {byte_position_start}")

            dtype = self.available_sections[keyword].get_ametek_type()
            offset = byte_position_start
            stride = np.uint64(
                self.available_sections[keyword].meta["i_data_type_size"] / 8
            )
            count = self.available_sections[keyword].get_ametek_count()
            shape = tuple(
                [
                    int(extent)
                    for extent in self.available_sections[keyword].get_ametek_shape()
                ]
            )
            logger.info(
                f"dtype {dtype}, offset {offset}, stride {stride}, count {count}, shape {shape}, file_size {self.file_size}"
            )
            data = get_memory_mapped_data(
                self.file_path, dtype, offset, stride, count
            )  # type ignore
            if data is not None:
                shape = self.available_sections[keyword].get_ametek_shape()
                unit = self.available_sections[keyword].meta["wc_data_unit"]
                # be careful with reshaping, above variable data is a 1d np.ndarray
                if len(shape) == 2:
                    if f"{np_uint16_to_string(unit)}" != "%/100":
                        clean_unit = f"{np_uint16_to_string(unit)}"
                    else:
                        clean_unit = "percent_per_100"
                    if int(shape[1]) == 1:
                        # e.g. "Mass" section, memory mapping yields (1*n,) should remain (n,)
                        # do not unnecessarily promote 1d arrays to 2d as this caused
                        # that previous versions of the library required a flattening
                        # of the ureg.magnitude return value which is unnecessary
                        return ureg.Quantity(np.asarray(data), clean_unit)
                    else:
                        # e.g. "Position" section, memory mapping yields (3*n,) but needs (n, 3)
                        return ureg.Quantity(
                            np.reshape(data, shape=(int(shape[0]), int(shape[1]))),
                            clean_unit,
                        )
                else:
                    raise ValueError("len(get_ametek_shape()) > 2 is not supported")
            else:
                logger.warning(f"Unable to get_named_quantity {keyword}")
        else:
            logger.warning(f"Unable to get_named_quantity {keyword}")

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge."""
        return self.get_named_quantity("Mass")

    def get_reconstructed_positions(self):
        """Read reconstructed positions."""
        return self.get_named_quantity("Position")
