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

"""Reader for *.analysisset XML/text-serialized objects from likely APSuite/IVAS."""

# This is an example how to extract pieces of information from *.analysisset files
# The example below shows how to extract ranging definitions.

# import xml.etree.ElementTree as ET
import html
import re

import flatdict as fd
import numpy as np
import xmltodict
from ase.data import atomic_numbers

from ifes_apt_tc_data_modeling.utils.custom_logging import logger
from ifes_apt_tc_data_modeling.utils.definitions import (
    MAX_NUMBER_OF_ATOMS_PER_ION,
    NEUTRON_NUMBER_FOR_ELEMENT,
)
from ifes_apt_tc_data_modeling.utils.molecular_ions import (
    get_chemical_symbols,
    isotope_to_hash,
)
from ifes_apt_tc_data_modeling.utils.nx_ion import NxIon


class ReadAnalysissetFileFormat:
    """Read *.analysisset file (format), extract ranging definitions as an example."""

    def __init__(self, file_path: str, verbose: bool = False):
        self.supported = False
        if not file_path.lower().endswith(".analysisset"):
            logger.warning(f"{file_path} is likely not an analysisset file")
            return
        self.supported = True
        self.file_path = file_path
        self.verbose = verbose
        self.analysisset: dict = {"ranges": {}, "ions": {}, "molecular_ions": []}
        self.read_analysisset_ranging_definitions()

    def read_analysisset_ranging_definitions(self):
        with open(self.file_path, encoding="utf-8") as fp:
            txt = fp.read()

            # example instances of *.analysisset files show that
            # the file is in general valid XML but internal XML block content
            # is text serialized; this in-block payload needs to be extracted first
            pattern = r"<Content>(.*?)</Content>"
            match = re.search(pattern, txt, re.DOTALL)
            del txt  # continue to work only with the payload in the Content block

            if not match:
                logger.warning(
                    f"{self.file_path} does not contain XML/text payload in a Content block"
                )
                return

            content = match.group(1)
            del pattern, match

            # unescape HTML entities to get proper XML
            content_xml = html.unescape(content)

            # example how to write to and read in again
            # UTF-16 file with BOM for debugging purposes
            # with open("content_xml.xml", "w", encoding="utf-16") as out_fp:
            #     out_fp.write(content_xml)
            # with open("content_output.xml", encoding="utf-16") as in_fp:
            #     xml = xmltodict.parse(in_fp.read())

            # parse unescaped XML to dictionary without I/O
            content_dict = xmltodict.parse(content_xml)
            del content_xml

            flattened = fd.FlatDict(content_dict, "/")
            del content_dict

            ion_type_info_record: dict[int, fd.FlatDict] = {}
            ion_type_info_range: dict[int, dict | list[dict]] = {}
            ion_idx = 0
            for key, obj in flattened["TopLevelRoiNodeContainer"].items():
                if not key.startswith("Ions"):
                    continue
                else:
                    if isinstance(obj, list):
                        for entry in obj:
                            ion_type_info_record[ion_idx] = fd.FlatDict(entry, "/")
                            ion_idx += 1
            ion_idx = 0
            for key, obj in flattened["TopLevelRoiNodeContainer"].items():
                if not key.startswith("RangesPerIon"):
                    continue
                else:
                    if isinstance(obj, list):
                        for entry in obj:
                            if isinstance(entry, dict):
                                if "RangeEntry" in entry:
                                    if isinstance(entry["RangeEntry"], list):
                                        ion_type_info_range[ion_idx] = entry[
                                            "RangeEntry"
                                        ]
                                        ion_idx += 1
                                    elif isinstance(entry["RangeEntry"], dict):
                                        ion_type_info_range[ion_idx] = entry[
                                            "RangeEntry"
                                        ]
                                        ion_idx += 1
            del flattened
            if ion_type_info_record.keys() != ion_type_info_range.keys():
                logger.warning(
                    f"{self.file_path} unable to extract valid ion definition and ranging data"
                )
                return

            for idx in range(0, len(ion_type_info_record)):
                key_name = "Formula/Element/Name"
                ivec = []  # used for all ranges ion_type_info_range[idx]
                if key_name in ion_type_info_record[idx]:  # element, e.g., "Zn"
                    data = ion_type_info_record[idx][key_name]
                    if isinstance(data, str):
                        if data in get_chemical_symbols():
                            proton_number = atomic_numbers[data]
                            neutron_number = NEUTRON_NUMBER_FOR_ELEMENT
                            key_count = "Formula/Element/Count"
                            if key_count in ion_type_info_record[idx]:
                                multiplier = int(ion_type_info_record[idx][key_count])
                                if multiplier > 0:
                                    ivec.extend(
                                        [isotope_to_hash(proton_number, neutron_number)]
                                        * multiplier
                                    )
                                else:
                                    logger.warning(
                                        f"{self.file_path} idx {idx} Formula/Element/Count not finite"
                                    )
                                    continue  # continue parsing at least some if possible
                            else:
                                logger.warning(
                                    f"{self.file_path} idx {idx} Formula/Element/Count not found"
                                )
                                continue
                        else:
                            logger.warning(
                                f"{self.file_path} idx {idx} {data} not a chemical symbol"
                            )
                            continue
                elif "Formula/Element" in ion_type_info_record[idx]:
                    data = ion_type_info_record[idx]["Formula/Element"]
                    if isinstance(data, list):
                        # e.g. molecular ions  [{'Name': 'Si', 'Count': '1'}, {'Name': 'Fe', 'Count': '1'}]
                        if all(
                            tuple(dictionary.keys()) == ("Name", "Count")
                            for dictionary in data
                        ):
                            for dictionary in data:
                                if dictionary["Name"] in get_chemical_symbols():
                                    proton_number = atomic_numbers[dictionary["Name"]]
                                    neutron_number = NEUTRON_NUMBER_FOR_ELEMENT
                                    multiplier = int(dictionary["Count"])
                                    if multiplier > 0:
                                        ivec.extend(
                                            [
                                                isotope_to_hash(
                                                    proton_number, neutron_number
                                                )
                                            ]
                                            * multiplier
                                        )
                                    else:
                                        logger.warning(
                                            f"{self.file_path} idx {idx} Formula/Element/Count not finite"
                                        )
                                        break  # must not continue parsing molecular ion partly
                                else:
                                    logger.warning(
                                        f"{self.file_path} idx {idx} {dictionary['Name']} not a chemical symbol"
                                    )
                                    break
                        else:
                            logger.warning(
                                f"{self.file_path} idx {idx} Formula/Element dictionary malformed"
                            )
                            continue
                    else:
                        logger.warning(
                            f"{self.file_path} idx {idx} Formula/Element malformed"
                        )
                        continue
                else:
                    logger.warning(
                        f"{self.file_path} idx {idx} {ion_type_info_record[idx]} malformed"
                    )
                    continue
                if len(ivec) > 0:
                    ivec = np.sort(np.asarray(ivec, np.uint16))[::-1]
                    ivector = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,), np.uint16)
                    ivector[0 : len(ivec)] = ivec
                    if isinstance(ion_type_info_range[idx], dict):
                        if tuple(ion_type_info_range[idx].keys()) == ("Max", "Min"):
                            mq = [
                                float(ion_type_info_range[idx]["Min"]),
                                float(ion_type_info_range[idx]["Max"]),
                            ]
                            m_ion = NxIon(nuclide_hash=ivector, charge_state=0)
                            m_ion.add_range(mq[0], mq[1])
                            m_ion.apply_combinatorics()
                            # m_ion.report()
                            self.analysisset["molecular_ions"].append(m_ion)
                        else:
                            logger.warning(
                                f"{self.file_path} idx {idx} ion_type_info_range dictionary malformed"
                            )
                            continue
                    elif isinstance(ion_type_info_range[idx], list):
                        if all(
                            tuple(dictionary.keys()) == ("Max", "Min")
                            for dictionary in ion_type_info_range[idx]
                        ):
                            for dictionary in ion_type_info_range[idx]:
                                mq = [
                                    float(dictionary["Min"]),
                                    float(dictionary["Max"]),
                                ]
                                m_ion = NxIon(nuclide_hash=ivector, charge_state=0)
                                m_ion.add_range(mq[0], mq[1])
                                m_ion.apply_combinatorics()
                                # m_ion.report()
                                self.analysisset["molecular_ions"].append(m_ion)
                        else:
                            logger.warning(
                                f"{self.file_path} idx {idx} ion_type_info_range list malformed"
                            )
                            continue
                    else:
                        logger.warning(
                            f"{self.file_path} idx {idx} ion_type_info_range neither dict nor list"
                        )
                        continue
                else:
                    logger.warning(f"{self.file_path} idx {idx} ivec empty")
                    continue
