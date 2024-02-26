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

from ifes_apt_tc_data_modeling.apt.apt6_sections import AptFileSectionMetadata


# define which header and sections to expect in an *.apt file
# the sections below are referred to as branches in the commercial software
# APSuite enables users select prior I/O which of the sections to write
# normally a range file from based on the exchange with AMETEK and

# in addition there is usually a sextet of preceeding IEEE 32-bit floats
# to the position section. This sextett encodes the min,max bounds of the point
# cloud in x, y, z direction respectively
# F. M. M. de Oliveira reported cases where *.apt from mere flat test runs,
# i.e. software sessions during which no reconstruction is performed
# this sextett is absent

EXPECTED_SECTIONS = {}

# most tested sections

EXPECTED_SECTIONS["Position"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Position"].set_section_name("Position")
EXPECTED_SECTIONS["Position"].set_c_signature()
EXPECTED_SECTIONS["Position"].set_i_header_size(148 + 6 * 4)
# six 4 byte IEEE 32-bit floats are following immediately after a "Position"
# section trailing the actual position value array
EXPECTED_SECTIONS["Position"].set_i_header_version(2)
EXPECTED_SECTIONS["Position"].set_wc_section_type("Position")
EXPECTED_SECTIONS["Position"].set_i_section_version(1)
EXPECTED_SECTIONS["Position"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Position"].set_e_record_type(1)
EXPECTED_SECTIONS["Position"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Position"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Position"].set_i_record_size(12)
EXPECTED_SECTIONS["Position"].set_wc_data_unit("nm")
EXPECTED_SECTIONS["Position"].set_accepted_units(["nm"])

EXPECTED_SECTIONS["Mass"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Mass"].set_section_name("Mass")
EXPECTED_SECTIONS["Mass"].set_c_signature()
EXPECTED_SECTIONS["Mass"].set_i_header_size(148)
EXPECTED_SECTIONS["Mass"].set_i_header_version(2)
EXPECTED_SECTIONS["Mass"].set_wc_section_type("Mass")
EXPECTED_SECTIONS["Mass"].set_i_section_version(1)
EXPECTED_SECTIONS["Mass"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Mass"].set_e_record_type(1)
EXPECTED_SECTIONS["Mass"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Mass"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Mass"].set_i_record_size(4)
EXPECTED_SECTIONS["Mass"].set_wc_data_unit("Da")
EXPECTED_SECTIONS["Mass"].set_accepted_units(["da", "Da", "amu"])

# M. K\"uhbach detected there are at least APSuite versions which write
# files with Da and some with da as wcDataUnits argument, it seems there
# was a source code change, we use pint to detect these issues and warn

# sections which demand more testing with the parser

EXPECTED_SECTIONS["z"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["z"].set_section_name("z")
EXPECTED_SECTIONS["z"].set_c_signature()
EXPECTED_SECTIONS["z"].set_i_header_size(148)
EXPECTED_SECTIONS["z"].set_i_header_version(2)
EXPECTED_SECTIONS["z"].set_wc_section_type("z")
EXPECTED_SECTIONS["z"].set_i_section_version(1)
EXPECTED_SECTIONS["z"].set_e_relationship_type(1)
EXPECTED_SECTIONS["z"].set_e_record_type(1)
EXPECTED_SECTIONS["z"].set_e_record_data_type(1)
EXPECTED_SECTIONS["z"].set_i_data_type_size(64)
EXPECTED_SECTIONS["z"].set_i_record_size(8)
EXPECTED_SECTIONS["z"].set_wc_data_unit("ions")
EXPECTED_SECTIONS["z"].set_accepted_units(["ions"])

EXPECTED_SECTIONS["tof"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["tof"].set_section_name("tof")
EXPECTED_SECTIONS["tof"].set_c_signature()
EXPECTED_SECTIONS["tof"].set_i_header_size(148)
EXPECTED_SECTIONS["tof"].set_i_header_version(2)
EXPECTED_SECTIONS["tof"].set_wc_section_type("tof")
EXPECTED_SECTIONS["tof"].set_i_section_version(1)
EXPECTED_SECTIONS["tof"].set_e_relationship_type(1)
EXPECTED_SECTIONS["tof"].set_e_record_type(1)
EXPECTED_SECTIONS["tof"].set_e_record_data_type(3)
EXPECTED_SECTIONS["tof"].set_i_data_type_size(32)
EXPECTED_SECTIONS["tof"].set_i_record_size(4)
EXPECTED_SECTIONS["tof"].set_wc_data_unit("ns")
EXPECTED_SECTIONS["tof"].set_accepted_units(["ns"])

EXPECTED_SECTIONS["Voltage"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Voltage"].set_section_name("Voltage")
EXPECTED_SECTIONS["Voltage"].set_c_signature()
EXPECTED_SECTIONS["Voltage"].set_i_header_size(148)
EXPECTED_SECTIONS["Voltage"].set_i_header_version(2)
EXPECTED_SECTIONS["Voltage"].set_wc_section_type("Voltage")
EXPECTED_SECTIONS["Voltage"].set_i_section_version(1)
EXPECTED_SECTIONS["Voltage"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Voltage"].set_e_record_type(1)
EXPECTED_SECTIONS["Voltage"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Voltage"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Voltage"].set_i_record_size(4)
EXPECTED_SECTIONS["Voltage"].set_wc_data_unit("")

EXPECTED_SECTIONS["pulse"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["pulse"].set_section_name("pulse")
EXPECTED_SECTIONS["pulse"].set_c_signature()
EXPECTED_SECTIONS["pulse"].set_i_header_size(148)
EXPECTED_SECTIONS["pulse"].set_i_header_version(2)
EXPECTED_SECTIONS["pulse"].set_wc_section_type("pulse")
EXPECTED_SECTIONS["pulse"].set_i_section_version(1)
EXPECTED_SECTIONS["pulse"].set_e_relationship_type(1)
EXPECTED_SECTIONS["pulse"].set_e_record_type(1)
EXPECTED_SECTIONS["pulse"].set_e_record_data_type(3)
EXPECTED_SECTIONS["pulse"].set_i_data_type_size(32)
EXPECTED_SECTIONS["pulse"].set_i_record_size(4)
EXPECTED_SECTIONS["pulse"].set_wc_data_unit("")

EXPECTED_SECTIONS["freq"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["freq"].set_section_name("freq")
EXPECTED_SECTIONS["freq"].set_c_signature()
EXPECTED_SECTIONS["freq"].set_i_header_size(148)
EXPECTED_SECTIONS["freq"].set_i_header_version(2)
EXPECTED_SECTIONS["freq"].set_wc_section_type("freq")
EXPECTED_SECTIONS["freq"].set_i_section_version(1)
EXPECTED_SECTIONS["freq"].set_e_relationship_type(1)
EXPECTED_SECTIONS["freq"].set_e_record_type(1)
EXPECTED_SECTIONS["freq"].set_e_record_data_type(3)
EXPECTED_SECTIONS["freq"].set_i_data_type_size(32)
EXPECTED_SECTIONS["freq"].set_i_record_size(4)
EXPECTED_SECTIONS["freq"].set_wc_data_unit("Hz")
EXPECTED_SECTIONS["freq"].set_accepted_units(["Hz"])

EXPECTED_SECTIONS["tElapsed"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["tElapsed"].set_section_name("tElapsed")
EXPECTED_SECTIONS["tElapsed"].set_c_signature()
EXPECTED_SECTIONS["tElapsed"].set_i_header_size(148)
EXPECTED_SECTIONS["tElapsed"].set_i_header_version(2)
EXPECTED_SECTIONS["tElapsed"].set_wc_section_type("tElapsed")
EXPECTED_SECTIONS["tElapsed"].set_i_section_version(1)
EXPECTED_SECTIONS["tElapsed"].set_e_relationship_type(1)
EXPECTED_SECTIONS["tElapsed"].set_e_record_type(1)
EXPECTED_SECTIONS["tElapsed"].set_e_record_data_type(3)
EXPECTED_SECTIONS["tElapsed"].set_i_data_type_size(32)
EXPECTED_SECTIONS["tElapsed"].set_i_record_size(4)
EXPECTED_SECTIONS["tElapsed"].set_wc_data_unit("")

EXPECTED_SECTIONS["erate"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["erate"].set_section_name("erate")
EXPECTED_SECTIONS["erate"].set_c_signature()
EXPECTED_SECTIONS["erate"].set_i_header_size(148)
EXPECTED_SECTIONS["erate"].set_i_header_version(2)
EXPECTED_SECTIONS["erate"].set_wc_section_type("erate")
EXPECTED_SECTIONS["erate"].set_i_section_version(1)
EXPECTED_SECTIONS["erate"].set_e_relationship_type(1)
EXPECTED_SECTIONS["erate"].set_e_record_type(1)
EXPECTED_SECTIONS["erate"].set_e_record_data_type(3)
EXPECTED_SECTIONS["erate"].set_i_data_type_size(32)
EXPECTED_SECTIONS["erate"].set_i_record_size(4)
EXPECTED_SECTIONS["erate"].set_wc_data_unit("%/100")
EXPECTED_SECTIONS["erate"].set_accepted_units(["%/100"])

EXPECTED_SECTIONS["xstage"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["xstage"].set_section_name("xstage")
EXPECTED_SECTIONS["xstage"].set_c_signature()
EXPECTED_SECTIONS["xstage"].set_i_header_size(148)
EXPECTED_SECTIONS["xstage"].set_i_header_version(2)
EXPECTED_SECTIONS["xstage"].set_wc_section_type("xstage")
EXPECTED_SECTIONS["xstage"].set_i_section_version(1)
EXPECTED_SECTIONS["xstage"].set_e_relationship_type(1)
EXPECTED_SECTIONS["xstage"].set_e_record_type(1)
EXPECTED_SECTIONS["xstage"].set_e_record_data_type(1)
EXPECTED_SECTIONS["xstage"].set_i_data_type_size(32)
EXPECTED_SECTIONS["xstage"].set_i_record_size(4)
EXPECTED_SECTIONS["xstage"].set_wc_data_unit("")

EXPECTED_SECTIONS["ystage"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["ystage"].set_section_name("ystage")
EXPECTED_SECTIONS["ystage"].set_c_signature()
EXPECTED_SECTIONS["ystage"].set_i_header_size(148)
EXPECTED_SECTIONS["ystage"].set_i_header_version(2)
EXPECTED_SECTIONS["ystage"].set_wc_section_type("ystage")
EXPECTED_SECTIONS["ystage"].set_i_section_version(1)
EXPECTED_SECTIONS["ystage"].set_e_relationship_type(1)
EXPECTED_SECTIONS["ystage"].set_e_record_type(1)
EXPECTED_SECTIONS["ystage"].set_e_record_data_type(1)
EXPECTED_SECTIONS["ystage"].set_i_data_type_size(32)
EXPECTED_SECTIONS["ystage"].set_i_record_size(4)
EXPECTED_SECTIONS["ystage"].set_wc_data_unit("")

EXPECTED_SECTIONS["zstage"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["zstage"].set_section_name("zstage")
EXPECTED_SECTIONS["zstage"].set_c_signature()
EXPECTED_SECTIONS["zstage"].set_i_header_size(148)
EXPECTED_SECTIONS["zstage"].set_i_header_version(2)
EXPECTED_SECTIONS["zstage"].set_wc_section_type("zstage")
EXPECTED_SECTIONS["zstage"].set_i_section_version(1)
EXPECTED_SECTIONS["zstage"].set_e_relationship_type(1)
EXPECTED_SECTIONS["zstage"].set_e_record_type(1)
EXPECTED_SECTIONS["zstage"].set_e_record_data_type(1)
EXPECTED_SECTIONS["zstage"].set_i_data_type_size(32)
EXPECTED_SECTIONS["zstage"].set_i_record_size(4)
EXPECTED_SECTIONS["zstage"].set_wc_data_unit("")

EXPECTED_SECTIONS["tstage"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["tstage"].set_section_name("tstage")
EXPECTED_SECTIONS["tstage"].set_c_signature()
EXPECTED_SECTIONS["tstage"].set_i_header_size(148)
EXPECTED_SECTIONS["tstage"].set_i_header_version(2)
EXPECTED_SECTIONS["tstage"].set_wc_section_type("tstage")
EXPECTED_SECTIONS["tstage"].set_i_section_version(1)
EXPECTED_SECTIONS["tstage"].set_e_relationship_type(1)
EXPECTED_SECTIONS["tstage"].set_e_record_type(1)
EXPECTED_SECTIONS["tstage"].set_e_record_data_type(2)
EXPECTED_SECTIONS["tstage"].set_i_data_type_size(16)
EXPECTED_SECTIONS["tstage"].set_i_record_size(2)
EXPECTED_SECTIONS["tstage"].set_wc_data_unit("")

EXPECTED_SECTIONS["TargetErate"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["TargetErate"].set_section_name("TargetErate")
EXPECTED_SECTIONS["TargetErate"].set_c_signature()
EXPECTED_SECTIONS["TargetErate"].set_i_header_size(148)
EXPECTED_SECTIONS["TargetErate"].set_i_header_version(2)
EXPECTED_SECTIONS["TargetErate"].set_wc_section_type("TargetErate")
EXPECTED_SECTIONS["TargetErate"].set_i_section_version(1)
EXPECTED_SECTIONS["TargetErate"].set_e_relationship_type(1)
EXPECTED_SECTIONS["TargetErate"].set_e_record_type(1)
EXPECTED_SECTIONS["TargetErate"].set_e_record_data_type(3)
EXPECTED_SECTIONS["TargetErate"].set_i_data_type_size(32)
EXPECTED_SECTIONS["TargetErate"].set_i_record_size(4)
EXPECTED_SECTIONS["TargetErate"].set_wc_data_unit("%/100")
EXPECTED_SECTIONS["TargetErate"].set_accepted_units(["%/100"])

EXPECTED_SECTIONS["TargetFlux"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["TargetFlux"].set_section_name("TargetFlux")
EXPECTED_SECTIONS["TargetFlux"].set_c_signature()
EXPECTED_SECTIONS["TargetFlux"].set_i_header_size(148)
EXPECTED_SECTIONS["TargetFlux"].set_i_header_version(2)
EXPECTED_SECTIONS["TargetFlux"].set_wc_section_type("TargetFlux")
EXPECTED_SECTIONS["TargetFlux"].set_i_section_version(1)
EXPECTED_SECTIONS["TargetFlux"].set_e_relationship_type(1)
EXPECTED_SECTIONS["TargetFlux"].set_e_record_type(1)
EXPECTED_SECTIONS["TargetFlux"].set_e_record_data_type(3)
EXPECTED_SECTIONS["TargetFlux"].set_i_data_type_size(32)
EXPECTED_SECTIONS["TargetFlux"].set_i_record_size(4)
EXPECTED_SECTIONS["TargetFlux"].set_wc_data_unit("")

EXPECTED_SECTIONS["pulseDelta"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["pulseDelta"].set_section_name("pulseDelta")
EXPECTED_SECTIONS["pulseDelta"].set_c_signature()
EXPECTED_SECTIONS["pulseDelta"].set_i_header_size(148)
EXPECTED_SECTIONS["pulseDelta"].set_i_header_version(2)
EXPECTED_SECTIONS["pulseDelta"].set_wc_section_type("pulseDelta")
EXPECTED_SECTIONS["pulseDelta"].set_i_section_version(1)
EXPECTED_SECTIONS["pulseDelta"].set_e_relationship_type(1)
EXPECTED_SECTIONS["pulseDelta"].set_e_record_type(1)
EXPECTED_SECTIONS["pulseDelta"].set_e_record_data_type(1)
EXPECTED_SECTIONS["pulseDelta"].set_i_data_type_size(16)
EXPECTED_SECTIONS["pulseDelta"].set_i_record_size(2)
EXPECTED_SECTIONS["pulseDelta"].set_wc_data_unit("")

EXPECTED_SECTIONS["Pres"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Pres"].set_section_name("Pres")
EXPECTED_SECTIONS["Pres"].set_c_signature()
EXPECTED_SECTIONS["Pres"].set_i_header_size(148)
EXPECTED_SECTIONS["Pres"].set_i_header_version(2)
EXPECTED_SECTIONS["Pres"].set_wc_section_type("Pres")
EXPECTED_SECTIONS["Pres"].set_i_section_version(1)
EXPECTED_SECTIONS["Pres"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Pres"].set_e_record_type(1)
EXPECTED_SECTIONS["Pres"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Pres"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Pres"].set_i_record_size(4)
EXPECTED_SECTIONS["Pres"].set_wc_data_unit("torr")
EXPECTED_SECTIONS["Pres"].set_accepted_units(["torr"])

EXPECTED_SECTIONS["VAnodeMon"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["VAnodeMon"].set_section_name("VAnodeMon")
EXPECTED_SECTIONS["VAnodeMon"].set_c_signature()
EXPECTED_SECTIONS["VAnodeMon"].set_i_header_size(148)
EXPECTED_SECTIONS["VAnodeMon"].set_i_header_version(2)
EXPECTED_SECTIONS["VAnodeMon"].set_wc_section_type("VAnodeMon")
EXPECTED_SECTIONS["VAnodeMon"].set_i_section_version(1)
EXPECTED_SECTIONS["VAnodeMon"].set_e_relationship_type(1)
EXPECTED_SECTIONS["VAnodeMon"].set_e_record_type(1)
EXPECTED_SECTIONS["VAnodeMon"].set_e_record_data_type(3)
EXPECTED_SECTIONS["VAnodeMon"].set_i_data_type_size(32)
EXPECTED_SECTIONS["VAnodeMon"].set_i_record_size(4)
EXPECTED_SECTIONS["VAnodeMon"].set_wc_data_unit("")

EXPECTED_SECTIONS["Temp"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Temp"].set_section_name("Temp")
EXPECTED_SECTIONS["Temp"].set_c_signature()
EXPECTED_SECTIONS["Temp"].set_i_header_size(148)
EXPECTED_SECTIONS["Temp"].set_i_header_version(2)
EXPECTED_SECTIONS["Temp"].set_wc_section_type("Temp")
EXPECTED_SECTIONS["Temp"].set_i_section_version(1)
EXPECTED_SECTIONS["Temp"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Temp"].set_e_record_type(1)
EXPECTED_SECTIONS["Temp"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Temp"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Temp"].set_i_record_size(4)
EXPECTED_SECTIONS["Temp"].set_wc_data_unit("K")
EXPECTED_SECTIONS["Temp"].set_accepted_units(["K"])


EXPECTED_SECTIONS["AmbTemp"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["AmbTemp"].set_section_name("AmbTemp")
EXPECTED_SECTIONS["AmbTemp"].set_c_signature()
EXPECTED_SECTIONS["AmbTemp"].set_i_header_size(148)
EXPECTED_SECTIONS["AmbTemp"].set_i_header_version(2)
EXPECTED_SECTIONS["AmbTemp"].set_wc_section_type("AmbTemp")
EXPECTED_SECTIONS["AmbTemp"].set_i_section_version(1)
EXPECTED_SECTIONS["AmbTemp"].set_e_relationship_type(1)
EXPECTED_SECTIONS["AmbTemp"].set_e_record_type(1)
EXPECTED_SECTIONS["AmbTemp"].set_e_record_data_type(3)
EXPECTED_SECTIONS["AmbTemp"].set_i_data_type_size(32)
EXPECTED_SECTIONS["AmbTemp"].set_i_record_size(4)
EXPECTED_SECTIONS["AmbTemp"].set_wc_data_unit("C")
EXPECTED_SECTIONS["AmbTemp"].set_accepted_units(["C"])

EXPECTED_SECTIONS["laserx"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["laserx"].set_section_name("laserx")
EXPECTED_SECTIONS["laserx"].set_c_signature()
EXPECTED_SECTIONS["laserx"].set_i_header_size(148)
EXPECTED_SECTIONS["laserx"].set_i_header_version(2)
EXPECTED_SECTIONS["laserx"].set_wc_section_type("laserx")
EXPECTED_SECTIONS["laserx"].set_i_section_version(1)
EXPECTED_SECTIONS["laserx"].set_e_relationship_type(1)
EXPECTED_SECTIONS["laserx"].set_e_record_type(1)
EXPECTED_SECTIONS["laserx"].set_e_record_data_type(1)
EXPECTED_SECTIONS["laserx"].set_i_data_type_size(32)
EXPECTED_SECTIONS["laserx"].set_i_record_size(4)
EXPECTED_SECTIONS["laserx"].set_wc_data_unit("")

EXPECTED_SECTIONS["lasery"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["lasery"].set_section_name("lasery")
EXPECTED_SECTIONS["lasery"].set_c_signature()
EXPECTED_SECTIONS["lasery"].set_i_header_size(148)
EXPECTED_SECTIONS["lasery"].set_i_header_version(2)
EXPECTED_SECTIONS["lasery"].set_wc_section_type("lasery")
EXPECTED_SECTIONS["lasery"].set_i_section_version(1)
EXPECTED_SECTIONS["lasery"].set_e_relationship_type(1)
EXPECTED_SECTIONS["lasery"].set_e_record_type(1)
EXPECTED_SECTIONS["lasery"].set_e_record_data_type(1)
EXPECTED_SECTIONS["lasery"].set_i_data_type_size(32)
EXPECTED_SECTIONS["lasery"].set_i_record_size(4)
EXPECTED_SECTIONS["lasery"].set_wc_data_unit("")

EXPECTED_SECTIONS["laserz"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["laserz"].set_section_name("laserz")
EXPECTED_SECTIONS["laserz"].set_c_signature()
EXPECTED_SECTIONS["laserz"].set_i_header_size(148)
EXPECTED_SECTIONS["laserz"].set_i_header_version(2)
EXPECTED_SECTIONS["laserz"].set_wc_section_type("laserz")
EXPECTED_SECTIONS["laserz"].set_i_section_version(1)
EXPECTED_SECTIONS["laserz"].set_e_relationship_type(1)
EXPECTED_SECTIONS["laserz"].set_e_record_type(1)
EXPECTED_SECTIONS["laserz"].set_e_record_data_type(1)
EXPECTED_SECTIONS["laserz"].set_i_data_type_size(32)
EXPECTED_SECTIONS["laserz"].set_i_record_size(4)
EXPECTED_SECTIONS["laserz"].set_wc_data_unit("")

EXPECTED_SECTIONS["laserpower"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["laserpower"].set_section_name("laserpower")
EXPECTED_SECTIONS["laserpower"].set_c_signature()
EXPECTED_SECTIONS["laserpower"].set_i_header_size(148)
EXPECTED_SECTIONS["laserpower"].set_i_header_version(2)
EXPECTED_SECTIONS["laserpower"].set_wc_section_type("laserpower")
EXPECTED_SECTIONS["laserpower"].set_i_section_version(1)
EXPECTED_SECTIONS["laserpower"].set_e_relationship_type(1)
EXPECTED_SECTIONS["laserpower"].set_e_record_type(1)
EXPECTED_SECTIONS["laserpower"].set_e_record_data_type(3)
EXPECTED_SECTIONS["laserpower"].set_i_data_type_size(32)
EXPECTED_SECTIONS["laserpower"].set_i_record_size(4)
EXPECTED_SECTIONS["laserpower"].set_wc_data_unit("")

EXPECTED_SECTIONS["FractureGuard"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["FractureGuard"].set_section_name("FractureGuard")
EXPECTED_SECTIONS["FractureGuard"].set_c_signature()
EXPECTED_SECTIONS["FractureGuard"].set_i_header_size(148)
EXPECTED_SECTIONS["FractureGuard"].set_i_header_version(2)
EXPECTED_SECTIONS["FractureGuard"].set_wc_section_type("FractureGuard")
EXPECTED_SECTIONS["FractureGuard"].set_i_section_version(1)
EXPECTED_SECTIONS["FractureGuard"].set_e_relationship_type(1)
EXPECTED_SECTIONS["FractureGuard"].set_e_record_type(1)
EXPECTED_SECTIONS["FractureGuard"].set_e_record_data_type(2)
EXPECTED_SECTIONS["FractureGuard"].set_i_data_type_size(16)
EXPECTED_SECTIONS["FractureGuard"].set_i_record_size(2)
EXPECTED_SECTIONS["FractureGuard"].set_wc_data_unit("")

EXPECTED_SECTIONS["Noise"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Noise"].set_section_name("Noise")
EXPECTED_SECTIONS["Noise"].set_c_signature()
EXPECTED_SECTIONS["Noise"].set_i_header_size(148)
EXPECTED_SECTIONS["Noise"].set_i_header_version(2)
EXPECTED_SECTIONS["Noise"].set_wc_section_type("Noise")
EXPECTED_SECTIONS["Noise"].set_i_section_version(1)
EXPECTED_SECTIONS["Noise"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Noise"].set_e_record_type(1)
EXPECTED_SECTIONS["Noise"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Noise"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Noise"].set_i_record_size(4)
EXPECTED_SECTIONS["Noise"].set_wc_data_unit("ions")
EXPECTED_SECTIONS["Noise"].set_accepted_units(["", "ions"])

EXPECTED_SECTIONS["Uniformity"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Uniformity"].set_section_name("Uniformity")
EXPECTED_SECTIONS["Uniformity"].set_c_signature()
EXPECTED_SECTIONS["Uniformity"].set_i_header_size(148)
EXPECTED_SECTIONS["Uniformity"].set_i_header_version(2)
EXPECTED_SECTIONS["Uniformity"].set_wc_section_type("Uniformity")
EXPECTED_SECTIONS["Uniformity"].set_i_section_version(1)
EXPECTED_SECTIONS["Uniformity"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Uniformity"].set_e_record_type(1)
EXPECTED_SECTIONS["Uniformity"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Uniformity"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Uniformity"].set_i_record_size(4)
EXPECTED_SECTIONS["Uniformity"].set_wc_data_unit("")

EXPECTED_SECTIONS["tofc"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["tofc"].set_section_name("tofc")
EXPECTED_SECTIONS["tofc"].set_c_signature()
EXPECTED_SECTIONS["tofc"].set_i_header_size(148)
EXPECTED_SECTIONS["tofc"].set_i_header_version(2)
EXPECTED_SECTIONS["tofc"].set_wc_section_type("tofc")
EXPECTED_SECTIONS["tofc"].set_i_section_version(1)
EXPECTED_SECTIONS["tofc"].set_e_relationship_type(1)
EXPECTED_SECTIONS["tofc"].set_e_record_type(1)
EXPECTED_SECTIONS["tofc"].set_e_record_data_type(3)
EXPECTED_SECTIONS["tofc"].set_i_data_type_size(32)
EXPECTED_SECTIONS["tofc"].set_i_record_size(4)
EXPECTED_SECTIONS["tofc"].set_wc_data_unit("ns")
EXPECTED_SECTIONS["tofc"].set_accepted_units(["ns"])

EXPECTED_SECTIONS["tofb"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["tofb"].set_section_name("tofb")
EXPECTED_SECTIONS["tofb"].set_c_signature()
EXPECTED_SECTIONS["tofb"].set_i_header_size(148)
EXPECTED_SECTIONS["tofb"].set_i_header_version(2)
EXPECTED_SECTIONS["tofb"].set_wc_section_type("tofb")
EXPECTED_SECTIONS["tofb"].set_i_section_version(1)
EXPECTED_SECTIONS["tofb"].set_e_relationship_type(1)
EXPECTED_SECTIONS["tofb"].set_e_record_type(1)
EXPECTED_SECTIONS["tofb"].set_e_record_data_type(3)
EXPECTED_SECTIONS["tofb"].set_i_data_type_size(32)
EXPECTED_SECTIONS["tofb"].set_i_record_size(4)
EXPECTED_SECTIONS["tofb"].set_wc_data_unit("ns")
EXPECTED_SECTIONS["tofb"].set_accepted_units(["ns"])

EXPECTED_SECTIONS["xs"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["xs"].set_section_name("xs")
EXPECTED_SECTIONS["xs"].set_c_signature()
EXPECTED_SECTIONS["xs"].set_i_header_size(148)
EXPECTED_SECTIONS["xs"].set_i_header_version(2)
EXPECTED_SECTIONS["xs"].set_wc_section_type("xs")
EXPECTED_SECTIONS["xs"].set_i_section_version(1)
EXPECTED_SECTIONS["xs"].set_e_relationship_type(1)
EXPECTED_SECTIONS["xs"].set_e_record_type(1)
EXPECTED_SECTIONS["xs"].set_e_record_data_type(3)
EXPECTED_SECTIONS["xs"].set_i_data_type_size(32)
EXPECTED_SECTIONS["xs"].set_i_record_size(4)
EXPECTED_SECTIONS["xs"].set_wc_data_unit("nm")
EXPECTED_SECTIONS["xs"].set_accepted_units(["nm"])

EXPECTED_SECTIONS["ys"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["ys"].set_section_name("ys")
EXPECTED_SECTIONS["ys"].set_c_signature()
EXPECTED_SECTIONS["ys"].set_i_header_size(148)
EXPECTED_SECTIONS["ys"].set_i_header_version(2)
EXPECTED_SECTIONS["ys"].set_wc_section_type("ys")
EXPECTED_SECTIONS["ys"].set_i_section_version(1)
EXPECTED_SECTIONS["ys"].set_e_relationship_type(1)
EXPECTED_SECTIONS["ys"].set_e_record_type(1)
EXPECTED_SECTIONS["ys"].set_e_record_data_type(3)
EXPECTED_SECTIONS["ys"].set_i_data_type_size(32)
EXPECTED_SECTIONS["ys"].set_i_record_size(4)
EXPECTED_SECTIONS["ys"].set_wc_data_unit("nm")
EXPECTED_SECTIONS["ys"].set_accepted_units(["nm"])

EXPECTED_SECTIONS["zs"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["zs"].set_section_name("zs")
EXPECTED_SECTIONS["zs"].set_c_signature()
EXPECTED_SECTIONS["zs"].set_i_header_size(148)
EXPECTED_SECTIONS["zs"].set_i_header_version(2)
EXPECTED_SECTIONS["zs"].set_wc_section_type("zs")
EXPECTED_SECTIONS["zs"].set_i_section_version(1)
EXPECTED_SECTIONS["zs"].set_e_relationship_type(1)
EXPECTED_SECTIONS["zs"].set_e_record_type(1)
EXPECTED_SECTIONS["zs"].set_e_record_data_type(3)
EXPECTED_SECTIONS["zs"].set_i_data_type_size(32)
EXPECTED_SECTIONS["zs"].set_i_record_size(4)
EXPECTED_SECTIONS["zs"].set_wc_data_unit("nm")
EXPECTED_SECTIONS["zs"].set_accepted_units(["nm"])

EXPECTED_SECTIONS["rTip"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["rTip"].set_section_name("rTip")
EXPECTED_SECTIONS["rTip"].set_c_signature()
EXPECTED_SECTIONS["rTip"].set_i_header_size(148)
EXPECTED_SECTIONS["rTip"].set_i_header_version(2)
EXPECTED_SECTIONS["rTip"].set_wc_section_type("rTip")
EXPECTED_SECTIONS["rTip"].set_i_section_version(1)
EXPECTED_SECTIONS["rTip"].set_e_relationship_type(1)
EXPECTED_SECTIONS["rTip"].set_e_record_type(1)
EXPECTED_SECTIONS["rTip"].set_e_record_data_type(3)
EXPECTED_SECTIONS["rTip"].set_i_data_type_size(32)
EXPECTED_SECTIONS["rTip"].set_i_record_size(4)
EXPECTED_SECTIONS["rTip"].set_wc_data_unit("nm")
EXPECTED_SECTIONS["rTip"].set_accepted_units(["nm"])

EXPECTED_SECTIONS["zApex"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["zApex"].set_section_name("zApex")
EXPECTED_SECTIONS["zApex"].set_c_signature()
EXPECTED_SECTIONS["zApex"].set_i_header_size(148)
EXPECTED_SECTIONS["zApex"].set_i_header_version(2)
EXPECTED_SECTIONS["zApex"].set_wc_section_type("zApex")
EXPECTED_SECTIONS["zApex"].set_i_section_version(1)
EXPECTED_SECTIONS["zApex"].set_e_relationship_type(1)
EXPECTED_SECTIONS["zApex"].set_e_record_type(1)
EXPECTED_SECTIONS["zApex"].set_e_record_data_type(3)
EXPECTED_SECTIONS["zApex"].set_i_data_type_size(32)
EXPECTED_SECTIONS["zApex"].set_i_record_size(4)
EXPECTED_SECTIONS["zApex"].set_wc_data_unit("nm")
EXPECTED_SECTIONS["zApex"].set_accepted_units(["nm"])

EXPECTED_SECTIONS["zSphereCorr"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["zSphereCorr"].set_section_name("zSphereCorr")
EXPECTED_SECTIONS["zSphereCorr"].set_c_signature()
EXPECTED_SECTIONS["zSphereCorr"].set_i_header_size(148)
EXPECTED_SECTIONS["zSphereCorr"].set_i_header_version(2)
EXPECTED_SECTIONS["zSphereCorr"].set_wc_section_type("zSphereCorr")
EXPECTED_SECTIONS["zSphereCorr"].set_i_section_version(1)
EXPECTED_SECTIONS["zSphereCorr"].set_e_relationship_type(1)
EXPECTED_SECTIONS["zSphereCorr"].set_e_record_type(1)
EXPECTED_SECTIONS["zSphereCorr"].set_e_record_data_type(3)
EXPECTED_SECTIONS["zSphereCorr"].set_i_data_type_size(32)
EXPECTED_SECTIONS["zSphereCorr"].set_i_record_size(4)
EXPECTED_SECTIONS["zSphereCorr"].set_wc_data_unit("nm")
EXPECTED_SECTIONS["zSphereCorr"].set_accepted_units(["nm"])

# the next three seem to (have been /were used for AMETEK development purposes
EXPECTED_SECTIONS["Position_0"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Position_0"].set_section_name("Position_0")
EXPECTED_SECTIONS["Position_0"].set_c_signature()
EXPECTED_SECTIONS["Position_0"].set_i_header_size(148 + 6 * 4)
EXPECTED_SECTIONS["Position_0"].set_i_header_version(2)
EXPECTED_SECTIONS["Position_0"].set_wc_section_type("Position_0")
EXPECTED_SECTIONS["Position_0"].set_i_section_version(1)
EXPECTED_SECTIONS["Position_0"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Position_0"].set_e_record_type(1)
EXPECTED_SECTIONS["Position_0"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Position_0"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Position_0"].set_i_record_size(12)
EXPECTED_SECTIONS["Position_0"].set_wc_data_unit("")

EXPECTED_SECTIONS["Position_1"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Position_1"].set_section_name("Position_1")
EXPECTED_SECTIONS["Position_1"].set_c_signature()
EXPECTED_SECTIONS["Position_1"].set_i_header_size(148 + 6 * 4)
EXPECTED_SECTIONS["Position_1"].set_i_header_version(2)
EXPECTED_SECTIONS["Position_1"].set_wc_section_type("Position_1")
EXPECTED_SECTIONS["Position_1"].set_i_section_version(1)
EXPECTED_SECTIONS["Position_1"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Position_1"].set_e_record_type(1)
EXPECTED_SECTIONS["Position_1"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Position_1"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Position_1"].set_i_record_size(12)
EXPECTED_SECTIONS["Position_1"].set_wc_data_unit("")

EXPECTED_SECTIONS["Position_2"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Position_2"].set_section_name("Position_2")
EXPECTED_SECTIONS["Position_2"].set_c_signature()
EXPECTED_SECTIONS["Position_2"].set_i_header_size(148 + 6 * 4)
EXPECTED_SECTIONS["Position_2"].set_i_header_version(2)
EXPECTED_SECTIONS["Position_2"].set_wc_section_type("Position_2")
EXPECTED_SECTIONS["Position_2"].set_i_section_version(1)
EXPECTED_SECTIONS["Position_2"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Position_2"].set_e_record_type(1)
EXPECTED_SECTIONS["Position_2"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Position_2"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Position_2"].set_i_record_size(12)
EXPECTED_SECTIONS["Position_2"].set_wc_data_unit("")

EXPECTED_SECTIONS["XDet_mm"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["XDet_mm"].set_section_name("XDet_mm")
EXPECTED_SECTIONS["XDet_mm"].set_c_signature()
EXPECTED_SECTIONS["XDet_mm"].set_i_header_size(148)
EXPECTED_SECTIONS["XDet_mm"].set_i_header_version(2)
EXPECTED_SECTIONS["XDet_mm"].set_wc_section_type("XDet_mm")
EXPECTED_SECTIONS["XDet_mm"].set_i_section_version(1)
EXPECTED_SECTIONS["XDet_mm"].set_e_relationship_type(1)
EXPECTED_SECTIONS["XDet_mm"].set_e_record_type(1)
EXPECTED_SECTIONS["XDet_mm"].set_e_record_data_type(3)
EXPECTED_SECTIONS["XDet_mm"].set_i_data_type_size(32)
EXPECTED_SECTIONS["XDet_mm"].set_i_record_size(4)
EXPECTED_SECTIONS["XDet_mm"].set_wc_data_unit("mm")
EXPECTED_SECTIONS["XDet_mm"].set_accepted_units(["mm"])

EXPECTED_SECTIONS["YDet_mm"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["YDet_mm"].set_section_name("YDet_mm")
EXPECTED_SECTIONS["YDet_mm"].set_c_signature()
EXPECTED_SECTIONS["YDet_mm"].set_i_header_size(148)
EXPECTED_SECTIONS["YDet_mm"].set_i_header_version(2)
EXPECTED_SECTIONS["YDet_mm"].set_wc_section_type("YDet_mm")
EXPECTED_SECTIONS["YDet_mm"].set_i_section_version(1)
EXPECTED_SECTIONS["YDet_mm"].set_e_relationship_type(1)
EXPECTED_SECTIONS["YDet_mm"].set_e_record_type(1)
EXPECTED_SECTIONS["YDet_mm"].set_e_record_data_type(3)
EXPECTED_SECTIONS["YDet_mm"].set_i_data_type_size(32)
EXPECTED_SECTIONS["YDet_mm"].set_i_record_size(4)
EXPECTED_SECTIONS["YDet_mm"].set_wc_data_unit("mm")
EXPECTED_SECTIONS["YDet_mm"].set_accepted_units(["mm"])

EXPECTED_SECTIONS["Multiplicity"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Multiplicity"].set_section_name("Multiplicity")
EXPECTED_SECTIONS["Multiplicity"].set_c_signature()
EXPECTED_SECTIONS["Multiplicity"].set_i_header_size(148)
EXPECTED_SECTIONS["Multiplicity"].set_i_header_version(2)
EXPECTED_SECTIONS["Multiplicity"].set_wc_section_type("Multiplicity")
EXPECTED_SECTIONS["Multiplicity"].set_i_section_version(1)
EXPECTED_SECTIONS["Multiplicity"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Multiplicity"].set_e_record_type(1)
EXPECTED_SECTIONS["Multiplicity"].set_e_record_data_type(1)
EXPECTED_SECTIONS["Multiplicity"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Multiplicity"].set_i_record_size(4)
EXPECTED_SECTIONS["Multiplicity"].set_wc_data_unit("")

EXPECTED_SECTIONS["Vap"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Vap"].set_section_name("Vap")
EXPECTED_SECTIONS["Vap"].set_c_signature()
EXPECTED_SECTIONS["Vap"].set_i_header_size(148)
EXPECTED_SECTIONS["Vap"].set_i_header_version(2)
EXPECTED_SECTIONS["Vap"].set_wc_section_type("Vap")
EXPECTED_SECTIONS["Vap"].set_i_section_version(1)
EXPECTED_SECTIONS["Vap"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Vap"].set_e_record_type(1)
EXPECTED_SECTIONS["Vap"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Vap"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Vap"].set_i_record_size(4)
EXPECTED_SECTIONS["Vap"].set_wc_data_unit("V")
EXPECTED_SECTIONS["Vap"].set_accepted_units(["V"])

EXPECTED_SECTIONS["Detector Coordinates"] = AptFileSectionMetadata()
EXPECTED_SECTIONS["Detector Coordinates"].set_section_name(
    "Detector Coordinates")
EXPECTED_SECTIONS["Detector Coordinates"].set_c_signature()
EXPECTED_SECTIONS["Detector Coordinates"].set_i_header_size(148)
EXPECTED_SECTIONS["Detector Coordinates"].set_i_header_version(2)
EXPECTED_SECTIONS["Detector Coordinates"].set_wc_section_type(
    "Detector Coordinates")
EXPECTED_SECTIONS["Detector Coordinates"].set_i_section_version(1)
EXPECTED_SECTIONS["Detector Coordinates"].set_e_relationship_type(1)
EXPECTED_SECTIONS["Detector Coordinates"].set_e_record_type(1)
EXPECTED_SECTIONS["Detector Coordinates"].set_e_record_data_type(3)
EXPECTED_SECTIONS["Detector Coordinates"].set_i_data_type_size(32)
EXPECTED_SECTIONS["Detector Coordinates"].set_i_record_size(8)
EXPECTED_SECTIONS["Detector Coordinates"].set_wc_data_unit("mm")
EXPECTED_SECTIONS["Detector Coordinates"].set_accepted_units(["mm"])

# there is at least a hint that some time ago this section"s wcDataUnit was ""


# sections whose formatting and purpose is completely unclear

# Var44 ? magic-in-action? M. K\"uhbach, could well for AMETEK development only
# EXPECTED_SECTIONS["Var44"].set_section_name("Var44")


# deprecated sections or sections with detected inconsistencies across versions
# Vref vs Voltage branch issue
EXPECTED_SECTIONS["Vref"] = EXPECTED_SECTIONS["Voltage"]
EXPECTED_SECTIONS["Vref"].set_wc_data_unit("V")
EXPECTED_SECTIONS["Vref"].set_accepted_units(["V"])

# pulseDelta vs Delta Pulse issue
# at least in one case a section Delta Pulse appeared
# at least in one case a section Epos ToF appeared

# other comments and issues
# Need to check APSuite version and build number
# "Voltage" section, M. K\"uhbach: expect to have units but example *.apt files
# from flat test do not encode a unit in the "Voltage" section
