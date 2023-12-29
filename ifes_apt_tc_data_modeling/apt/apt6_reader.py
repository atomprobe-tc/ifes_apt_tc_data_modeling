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

# pylint: disable=no-member,duplicate-code

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

    def __init__(self, filename: str):
        if (len(filename) <= 4) or (filename.lower().endswith(".apt") is False):
            raise ImportError("WARNING::APT file incorrect filename ending or file type!")
        self.filename = filename

        self.filesize = os.path.getsize(self.filename)
        print(f"Reading {self.filename} which is {self.filesize} B")

        self.header_section = None
        self.byte_offsets: dict = {}
        self.available_sections: dict = {}
        # where do sections start bytes from beginning of the file?
        self.apt: dict = {}

        self.parse_file_structure()

    def parse_file_structure(self):
        """Parse APT file header plus flat collection of metadata/data pairs.

        Each pair has a so-called section header and a corresponding raw data
        block. Section headers detail the content of the immediately trailing
        raw data block.
        An APT file can store none, some, or all of the possible sections.
        Furthermore the file can contain additional pieces of information
        which this parser cannot read-out because the APT format is maintained
        by AMETEK. The AMETEK source code is the only reliable source of
        information about which content the sections encode and how these
        get formatted when exporting an APT file from APSuite for a specific
        version and build number and type of experiment plus
        combinations of settings.
        Parse header of the file and check which parsable sections the file
        contains get the byte offsets of the sections from the beginning
        relative to the start/beginning of the file.
        """
        self.byte_offsets = {}
        self.header_section = None
        self.available_sections = {}

        with open(self.filename, 'rb') as file_handle:
            self.dummy_header = AptFileHeaderMetadata()
            found_header = np.fromfile(file_handle,
                                       self.dummy_header.get_numpy_struct(),
                                       count=1)

            assert self.dummy_header.matches(found_header), \
                f"Found an unexpectedly formatted/versioned header." \
                f"Create an issue to help us fix this!"
            print(f"File describes {found_header['llIonCount'][0]} ions")

            self.header_section = found_header
            self.byte_offsets['header'] = np.uint64(file_handle.tell())
            print(self.byte_offsets['header'])

            end_of_file_not_reached = b'yes'
            while end_of_file_not_reached != b'':
                # probe for end of file
                end_of_file_not_reached = file_handle.read(1)
                if end_of_file_not_reached != b'':
                    file_handle.seek(-1, os.SEEK_CUR)
                else:
                    print(f"End of file at {file_handle.tell()} B")
                    break

                dummy_section = AptFileSectionMetadata()
                found_section = np.fromfile(file_handle,
                                            dummy_section.get_numpy_struct(),
                                            count=1)
                keyword = np_uint16_to_string(
                    found_section['wcSectionType'][0])

                print(keyword)
                print(found_section)
                assert keyword not in self.available_sections, \
                    f"Found a duplicate of an already parsed section!" \
                    f"Create an issue to help us fix this!"

                if keyword not in ['Delta Pulse', 'Epos ToF']:
                    assert keyword in EXPECTED_SECTIONS, \
                        f"Found an unknown section, seems like an unknown/new branch!" \
                        f"Create an issue to help us fix this!"

                    metadata_section = EXPECTED_SECTIONS[keyword]
                    if metadata_section.matches(found_section) is True:
                        # assert metadata_section.matches(found_section), \
                        #     'Found an uninterpretable section! Please contact the \
                        #     development team to help us fixing this.'
                        self.available_sections[keyword] = metadata_section
                else:
                    print(f"Found an uninterpretable non-registered section."
                          f"Create an issue to help us fix this!, Parsing continues"
                          f"llByteCount {found_section['llByteCount'][0]} B")

                self.byte_offsets[keyword] = np.uint64(file_handle.tell())
                if keyword == 'Position':
                    # special case six IEEE 32-bit floats preceeding raw data
                    self.byte_offsets[keyword] += np.uint64(6 * 4)
                self.byte_offsets[keyword] += np.uint64(
                    found_section['llByteCount'][0])
                print(f"Byte offset for reading data for section: {keyword}"
                      f" {self.byte_offsets[keyword]} B")
                # print(file_handle.tell())
                file_handle.seek(self.byte_offsets[keyword], os.SEEK_SET)
                # print(file_handle.tell())

    # one convenience reader function for every known section
    # is useful because it structures the parsers, enables reading the file
    # partially and reduces main memory consumption during full parsing
    def get_header(self):
        """Report metadata in the header."""
        metadata_dict = {
            'cSignature':
                np_uint16_to_string(self.header_section['cSignature'][0]),
            'iHeaderSize':
                np.int32(self.header_section['iHeaderSize'][0]),
            'iHeaderVersion':
                np.int32(self.header_section['iHeaderVersion'][0]),
            'wcFilename':
                np_uint16_to_string(self.header_section['wcFilename'][0]),
            'ftCreationTime':
                np.uint64(self.header_section['ftCreationTime'][0]),
            'llIonCount':
                np.uint64(self.header_section['llIonCount'][0])}
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
        column_names = ['section']  # header
        assert 'Mass' in self.available_sections, \
            "Mass section not available to guide creation of the table header!"
        for key in self.available_sections['Mass'].get_metadata().keys():
            column_names.append(key)
        data_frame = pd.DataFrame(columns=column_names)

        for keyword, value in self.available_sections.items():
            row_dct = {'section': keyword}
            # print(f"{keyword}")
            row_dct = {**row_dct, **value.get_metadata()}
            # print(value.get_metadata())
            row_df = pd.DataFrame(row_dct, index=[0])
            # print(row_df)
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
            stride = np.uint64(
                self.available_sections[keyword].meta['i_data_type_size'] / 8)
            count = self.available_sections[keyword].get_ametek_count()

            data = get_memory_mapped_data(
                self.filename, dtype, offset, stride, count)

            shape = tuple(self.available_sections[keyword].get_ametek_shape())
            unit = self.available_sections[keyword].meta['wc_data_unit']

            return NxField(
                np.reshape(data, newshape=shape), np_uint16_to_string(unit))

        return NxField()

    def get_mass_to_charge_state_ratio(self):
        """Read mass-to-charge."""
        return self.get_named_quantity('Mass')

    def get_reconstructed_positions(self):
        """Read reconstructed positions."""
        return self.get_named_quantity('Position')

# test cases how to use the parser
# TEST_FILE_NAME = '70_50_50.apt'  # Xuyang Zhou's (MPIE) \
# apt = ReadAptFileFormat(TEST_FILE_NAME)
# print(apt.get_metadata_table())
# print(apt.get_header())
# xyz = apt.get_reconstructed_positions()
# equivalent to
# xyz = apt.get_metadata('Position')
# mq = apt.get_mass_to_charge_state_ratios()
# equivalent to
# mq = parsedFile.get_named_quantity('Mass')
