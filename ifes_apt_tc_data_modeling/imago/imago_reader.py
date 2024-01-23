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

"""Reader for *.analysis JSON-serialized objects from IVAS/Imago."""

# This is an example how to extract pieces of information from *.analysis files
# of Imago's IVAS analysis. Imago is the forerunner company of AMETEK/Cameca
# The example below shows how to extract ranging definitions.

import re
import xmltodict
import flatdict as fd
import numpy as np

from ase.data import atomic_numbers, chemical_symbols
from ifes_apt_tc_data_modeling.nexus.nx_ion import NxField, NxIon
from ifes_apt_tc_data_modeling.utils.definitions import \
    MAX_NUMBER_OF_ATOMS_PER_ION
from ifes_apt_tc_data_modeling.utils.molecular_ions import \
    get_chemical_symbols, isotope_to_hash


class ReadImagoAnalysisFileFormat():
    """Read *.analysis file (format), extract ranging definitions as an example."""

    def __init__(self, file_path: str):
        if (len(file_path) <= 9) or (file_path.lower().endswith(".analysis") is False):
            raise ImportError("WARNING::ANALYSIS file incorrect file_path ending or file type!")
        self.file_path = file_path
        self.imago: dict = {"ranges": {},
                            "ions": {},
                            "molecular_ions": []}
        self.read_imago_analysis_ranging_definitions()

    def read_imago_analysis_ranging_definitions(self):
        """Read *.analysis range file content."""
        # https://wilfried-grupe.de/Java_XMLEncoder.html
        # IVAS, the integrated visualization and analysis software/suite for
        # atom probe was a forerunner of the product that is now called APSuite
        # IVAS is a java application that uses JSON/XML object serialization
        # to store state files, many atom probe groups have still the resultant
        # state files within their projects, this function shows an example
        # how one can extract information from these files

        # unfortunately the serialized XML is in a format that when represented
        # in python mixes recurrent nested lists within dictionaries. Such a data
        # structure is not directly flattenable and hence many checks are required
        # I expect that this parser does not work out of the box for many
        # examples given that the formatting of object serialized content
        # depends heavily on the implementation of the host application (IVAS)
        # the main idea behind the example is to show that information can
        # be extracted and to motivate that nowadays one should use data
        # structures that are more conveniently parsable
        with open(self.file_path, "r", encoding="utf-8") as xmlf:
            xml = xmltodict.parse(xmlf.read())
            flt = fd.FlatDict(xml, "/")
            for entry in flt["java/object/void"]:
                # strategy is, walk the data structure and try to discard non-ranging content as early as possible
                for key, val in fd.FlatDict(entry, "/").items():
                    if (not isinstance(val, list)) or (key != "object/void"):
                        continue
                    for member in val:
                        if (not isinstance(member, dict)) \
                                or ("@method" not in member.keys()) \
                                or (member["@method"] != "add"):
                            continue
                        # print(">>>>>>>>>>>At the level of a molecular ion that can be so simple that it is just an element ion")
                        cand_dct = fd.FlatDict(member, "/")
                        # print(f">>>>> {cand_dct}")
                        all_reqs_exist = True
                        reqs = ["@method", "object/@id", "object/@class", "object/string", "object/boolean", "object/void"]
                        for req in reqs:
                            if req not in cand_dct.keys():
                                all_reqs_exist = False
                        if all_reqs_exist == False:
                            continue

                        if (not cand_dct["object/@id"].startswith("AtomDataRealRange")) \
                                or (cand_dct["object/@class"] != "com.imago.core.atomdata.AtomDataRealRange") \
                                or (not isinstance(cand_dct["object/void"], list)):
                            continue
                        for lst in cand_dct["object/void"]:
                            rng = fd.FlatDict(lst, "/")
                            element_symbol = []
                            mq = []
                            if "object/void/string" in rng.keys():
                                if rng["object/void/string"] in get_chemical_symbols():
                                    if "object/double" in rng.keys():
                                        mq = rng["object/double"][0:2]
                                        element_symbol.append(rng["object/void/string"])  # assuming multiplicity is one !
                            else:
                                if "object/void" in rng.keys():
                                    if isinstance(rng["object/void"], list):
                                        mq = rng["object/double"][0:2]
                                        element_symbol = []
                                        for block in rng["object/void"]:
                                            if isinstance(block, dict):
                                                if "@method" in block.keys() and "string" in block.keys() and "double" in block.keys():
                                                    if block["string"] in chemical_symbols:
                                                        for mult in np.arange(0, int(block["double"].split('.')[0])):
                                                            element_symbol.append(block["string"])
                            if (len(element_symbol) >= 1) and (len(mq) == 2):
                                # print(f"------------>{element_symbol}, {mq}")
                                ivec = []
                                for isotope in element_symbol:
                                    if isotope != "":
                                        prefix = re.findall("^[0-9]+", isotope)
                                        mass_number = 0
                                        if len(prefix) == 1:
                                            if int(prefix[0]) > 0:
                                                mass_number = int(prefix[0])
                                        suffix = re.findall("[0-9]+$", isotope)
                                        multiplier = 1
                                        if len(suffix) == 1:
                                            multiplier = int(suffix[0])
                                        symbol = isotope.replace(
                                            f"{mass_number}", "").replace(f"{multiplier}", "").replace(" ", "")
                                        if symbol in get_chemical_symbols():
                                            proton_number = atomic_numbers[symbol]
                                            neutron_number = 0
                                            if mass_number != 0:
                                                neutron_number = mass_number - proton_number
                                            ivec.extend([isotope_to_hash(proton_number, neutron_number)] * multiplier)
                                ivec = np.sort(np.asarray(ivec, np.uint16))[::-1]
                                ivector = np.zeros((MAX_NUMBER_OF_ATOMS_PER_ION,), np.uint16)
                                ivector[0:len(ivec)] = ivec

                                m_ion = NxIon(nuclide_hash=ivector, charge_state=0)
                                m_ion.add_range(float(mq[0]), float(mq[1]))
                                m_ion.comment = NxField(" ".join(element_symbol), "")
                                m_ion.apply_combinatorics()
                                m_ion.report()

                                self.imago["molecular_ions"].append(m_ion)
        print(f"{self.file_path} parsed successfully")
