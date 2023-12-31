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

from ifes_apt_tc_data_modeling.apt.apt6_utils import np_uint16_to_string
from ifes_apt_tc_data_modeling.apt.apt6_headers import AptFileHeaderMetadata
from ifes_apt_tc_data_modeling.apt.apt6_sections import AptFileSectionMetadata
from ifes_apt_tc_data_modeling.apt.apt6_sections_branches import EXPECTED_SECTIONS
from ifes_apt_tc_data_modeling.nexus.nx_field import NxField
from ifes_apt_tc_data_modeling.utils.mmapped_io import get_memory_mapped_data


class ReadAptFileFormat():
    """Read AMETEK's open exchange *.apt file format."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 4) or (file_path.lower().endswith(".apt") is False):
            raise ImportError("WARNING::APT file incorrect file_path ending or file type!")
        self.file_path = file_path
        self.file_size = os.path.getsize(self.file_path)
        print(f"Reading {self.file_path} which is {self.file_size} B")

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
            found_header = np.fromfile(fp,
                                       self.dummy_header.get_numpy_struct(),
                                       count=1)

            assert self.dummy_header.matches(found_header), \
                "Found an unexpectedly formatted header. Create an issue to help us fix this!"
            print(f"File describes {found_header['llIonCount'][0]} ions")

            self.header_section = found_header
            self.byte_offsets["header"] = np.uint64(fp.tell())
            print(f"Currently at byte_offset {self.byte_offsets['header']} B")

            end_of_file_not_reached = b"yes"
            while end_of_file_not_reached != b"":
                # probe for end of file
                end_of_file_not_reached = fp.read(1)
                if end_of_file_not_reached != b"":
                    fp.seek(-1, os.SEEK_CUR)
                else:
                    print(f"End of file at {fp.tell()} B")
                    break

                dummy_section = AptFileSectionMetadata()
                found_section = np.fromfile(fp, dummy_section.get_numpy_struct(), count=1)
                keyword = np_uint16_to_string(found_section["wcSectionType"][0])

                print(f"keyword: {keyword}, found_section: {found_section}")
                if keyword in self.available_sections:
                    raise ValueError("Found a duplicate of an already parsed section! "
                                     "Create an issue to help us fix this!")

                if keyword not in ["Delta Pulse", "Epos ToF"]:
                    if keyword not in EXPECTED_SECTIONS:
                        raise ValueError("Found an unknown section, seems like an unknown/new "
                                         "branch! Create an issue to help us fix this!")
                    metadata_section = EXPECTED_SECTIONS[keyword]
                    if metadata_section.matches(found_section) is True:
                        self.available_sections[keyword] = metadata_section
                else:
                    print(f"Found an uninterpretable non-registered section."
                          f"Create an issue to help us fix this!, Parsing continues"
                          f"llByteCount {found_section['llByteCount'][0]} B")

                self.byte_offsets[keyword] = np.uint64(fp.tell())
                if keyword == "Position":
                    # special case six IEEE 32-bit floats preceeding raw data
                    self.byte_offsets[keyword] += np.uint64(6 * 4)
                self.byte_offsets[keyword] += np.uint64(found_section["llByteCount"][0])
                print(f"Byte offset for reading data for section: {keyword}"
                      f" {self.byte_offsets[keyword]} B")
                fp.seek(self.byte_offsets[keyword], os.SEEK_SET)

    # one convenience reader function for every known section
    # is useful because it structures the parsers, enables reading the file
    # partially and reduces main memory consumption during full parsing
    def get_header(self):
        """Report metadata in the header."""
        metadata_dict = {
            "cSignature":
                np_uint16_to_string(self.header_section["cSignature"][0]),
            "iHeaderSize":
                np.int32(self.header_section["iHeaderSize"][0]),
            "iHeaderVersion":
                np.int32(self.header_section["iHeaderVersion"][0]),
            "wcFilename":
                np_uint16_to_string(self.header_section["wcFilename"][0]),
            "ftCreationTime":
                np.uint64(self.header_section["ftCreationTime"][0]),
            "llIonCount":
                np.uint64(self.header_section["llIonCount"][0])}
        # check e.g. https://gist.github.com/Mostafa-Hamdy-Elgiar/
        # 9714475f1b3bc224ea063af81566d873 repo
        # for converting Windows/MSDN time to Python time
        for key, value in iter(metadata_dict.items()):
            print(f"{key}: {value}")

    def get_metadata(self, keyword: str):
        """Report available metadata for quantity if it exists."""
        if (keyword in self.available_sections) and (keyword in self.byte_offsets):
            metadata_dict = self.available_sections[keyword].get_metadata()
            for key, value in iter(metadata_dict.items()):
                print(f"{key}: {value}")

    def get_metadata_table(self):
        """Create table from all metadata for each section."""
        column_names = ["section"]  # header
        if "Mass" not in self.available_sections:
            raise ValueError("Mass section not available to guide "
                             " creation of the table header!")
        for key in self.available_sections["Mass"].get_metadata().keys():
            column_names.append(key)
        data_frame = pd.DataFrame(columns=column_names)

        for keyword, value in self.available_sections.items():
            row_dct = {"section": keyword}
            row_dct = {**row_dct, **value.get_metadata()}
            # print(value.get_metadata())
            row_df = pd.DataFrame(row_dct, index=[0])
            data_frame = pd.concat([data_frame, row_df], ignore_index=True)

        data_frame.style.format(precision=3, thousands=",", decimal=".") \
            .format_index(str.upper, axis=1)
        return data_frame

    def get_named_quantity(self, keyword: str):
        """Read quantity with name in keyword from APT file if it exists."""
        if (keyword in self.available_sections) and (keyword in self.byte_offsets):
            byte_position_start = self.byte_offsets[keyword] \
                - self.available_sections[keyword].get_ametek_size()
            print(f"Reading section {keyword} at {byte_position_start}")

            dtype = self.available_sections[keyword].get_ametek_type()
            offset = byte_position_start
            stride = np.uint64(self.available_sections[keyword].meta["i_data_type_size"] / 8)
            count = self.available_sections[keyword].get_ametek_count()
            data = get_memory_mapped_data(self.file_path, dtype, offset, stride, count)
            shape = tuple(self.available_sections[keyword].get_ametek_shape())
            unit = self.available_sections[keyword].meta["wc_data_unit"]
            return NxField(np.reshape(data, newshape=shape), np_uint16_to_string(unit))

        return NxField()

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge."""
        return self.get_named_quantity("Mass")

    def get_reconstructed_positions(self):
        """Read reconstructed positions."""
        return self.get_named_quantity("Position")
