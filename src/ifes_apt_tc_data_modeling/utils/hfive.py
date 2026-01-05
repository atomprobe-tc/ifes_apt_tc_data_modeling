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

import numpy as np

from ifes_apt_tc_data_modeling.utils.pint_custom_unit_registry import ureg


def simple_hfive_file(fpw, idx, mion):
    """Write specific content in existent and opened HDF5 file pointed to by fpw with write access."""
    trg = f"/entry1/ion{idx}"
    grp = fpw.create_group(trg)
    grp.attrs["NX_class"] = "NXion"
    dst = fpw.create_dataset(f"{trg}/comment", data=mion.comment)
    # dst = fpw.create_dataset(f"{trg}/color", data=mion.color)
    # dst = fpw.create_dataset(f"{trg}/volume", dtype=np.float32, data=0.)
    # dst.attrs["units"] = f"{ureg.nanometer ** 3}"
    dst = fpw.create_dataset(
        f"{trg}/nuclide_hash",
        dtype=np.uint16,
        data=mion.nuclide_hash,
        chunks=True,
        compression="gzip",
        compression_opts=1,
    )
    dst = fpw.create_dataset(
        f"{trg}/nuclide_list",
        dtype=np.uint16,
        data=mion.nuclide_list,
        chunks=True,
        compression="gzip",
        compression_opts=1,
    )
    dst = fpw.create_dataset(
        f"{trg}/charge_state", dtype=np.int8, data=mion.charge_state
    )
    dst = fpw.create_dataset(f"{trg}/name", data=mion.name)
    dst = fpw.create_dataset(
        trg + "/mass_to_charge_range", dtype=np.float32, data=mion.ranges.magnitude
    )
    dst.attrs["units"] = f"{mion.ranges.units}"
    sub_group_name = f"{trg}/charge_state_analysis"
    sub_group = fpw.create_group(sub_group_name)
    sub_group.attrs["NX_class"] = "NXcharge_state_analysis"
    # config
    dst = fpw.create_dataset(
        f"{sub_group_name}/min_abundance",
        dtype=np.float64,
        data=mion.charge_state_model["min_abundance"],
    )
    # dst = fpw.create_dataset(f"{sub_group_name}/min_abundance_product", dtype=np.float64,
    #                          data=mion.charge_state_model["min_abundance_product"])
    dst = fpw.create_dataset(
        f"{sub_group_name}/min_half_life",
        dtype=np.float64,
        data=mion.charge_state_model["min_half_life"],
    )
    dst.attrs["units"] = f"{ureg.second}"
    dst = fpw.create_dataset(
        f"{sub_group_name}/sacrifice_isotopic_uniqueness",
        dtype=np.uint8,
        data=mion.charge_state_model["sacrifice_isotopic_uniqueness"],
    )
    opt_field_names = [
        "nuclide_hash",
        "charge_state",
        "mass",
        "natural_abundance_product",
        "shortest_half_life",
    ]
    all_opt_available = True
    for opt_field_name in opt_field_names:
        if opt_field_name not in mion.charge_state_model:
            all_opt_available = False
    if all_opt_available:
        dst = fpw.create_dataset(
            f"{sub_group_name}/nuclide_hash",
            dtype=np.uint16,
            data=mion.charge_state_model["nuclide_hash"],
            chunks=True,
            compression="gzip",
            compression_opts=1,
        )
        dst = fpw.create_dataset(
            f"{sub_group_name}/charge_state",
            dtype=np.int8,
            data=mion.charge_state_model["charge_state"],
        )
        dst = fpw.create_dataset(
            f"{sub_group_name}/mass",
            dtype=np.float64,
            data=mion.charge_state_model["mass"],
        )
        dst.attrs["units"] = f"{ureg.dalton}"
        dst = fpw.create_dataset(
            f"{sub_group_name}/natural_abundance_product",
            dtype=np.float64,
            data=mion.charge_state_model["natural_abundance_product"],
        )
        dst = fpw.create_dataset(
            f"{sub_group_name}/shortest_half_life",
            dtype=np.float64,
            data=mion.charge_state_model["shortest_half_life"],
        )
        dst.attrs["units"] = f"{ureg.second}"
