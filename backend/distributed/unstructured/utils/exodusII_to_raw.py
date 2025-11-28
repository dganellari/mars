#!/usr/bin/env python3

import netCDF4
import numpy as np
import sys
import os
import getopt

# import pdb

try:
    geom_t
except NameError:
    geom_t = np.float32
    idx_t = np.int32
    element_idx_t = np.int32


def exodusII_to_raw(input_mesh, output_folder):
    def mkdir(path):
        if not os.path.exists(path):
            os.makedirs(path)

    mkdir(output_folder)
    mkdir(f"{output_folder}/blocks")

    nc = netCDF4.Dataset(input_mesh)

    # -----------------------------------------------------------------------------
    # Coordinates
    # -----------------------------------------------------------------------------
    if "coord" in nc.variables:
        coords = nc.variables["coord"]
    else:
        coords = []
        if "coordx" in nc.variables:
            coordx = nc.variables["coordx"]
            coords.append(coordx)
        if "coordy" in nc.variables:
            coordy = nc.variables["coordy"]
            coords.append(coordy)
        if "coordz" in nc.variables:
            coordz = nc.variables["coordz"]
            coords.append(coordz)

        coords = np.array(coords)

    dims, nnodes = coords.shape

    coordnames = ["x", "y", "z", "t"]

    for i in range(0, dims):
        x = np.array(coords[i, :]).astype(geom_t)
        x.tofile(f"{output_folder}/{coordnames[i]}.raw")

    # -----------------------------------------------------------------------------
    # Time
    # -----------------------------------------------------------------------------
    n_time_steps = 1
    if "time_whole" in nc.variables:
        time_whole = nc.variables["time_whole"]
        n_time_steps = time_whole.shape[0]
        t = np.array(time_whole[:]).astype(np.float32)
        t.tofile(f"{output_folder}/time_whole.raw")

    print(f"n_time_steps = {n_time_steps}")

    # -----------------------------------------------------------------------------
    # Point data (nodal variables)
    # -----------------------------------------------------------------------------
    if "name_nod_var" in nc.variables:
        name_nod_var = nc.variables["name_nod_var"]
        nvars, __ = name_nod_var.shape
        print(f"Point data, nvars = {nvars}")

        point_data_dir = f"{output_folder}/point_data"
        mkdir(point_data_dir)

        nodal_prefix = "vals_nod_var"
        for i in range(0, nvars):
            var_key = f"{nodal_prefix}{i+1}"
            var = nc.variables[var_key]

            var_name = netCDF4.chartostring(name_nod_var[i, :])
            print(f" - {var_name}, dtype {var.dtype}")

            var_path_prefix = f"{point_data_dir}/{var_name}"

            if n_time_steps <= 1:
                path = f"{var_path_prefix}.raw"

                data = np.array(var[:])
                data.tofile(path)
            else:
                size_padding = int(np.ceil(np.log10(n_time_steps)))

                format_string = f"%s.%0.{size_padding}d.raw"

                for t_idx in range(0, n_time_steps):
                    data = np.array(var[t_idx, :])

                    path = format_string % (var_path_prefix, t_idx)
                    data.tofile(path)

    # -----------------------------------------------------------------------------
    # Side-to-node maps
    # -----------------------------------------------------------------------------
    def s2n_quad4():
        # Exodus QUAD sides: 1: (1,2), 2: (2,3), 3: (3,4), 4: (4,1)
        return [
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
        ]

    def s2n_hex8():
        # Exodus HEX faces:
        # 1: 1,2,6,5; 2: 2,3,7,6; 3: 3,4,8,7; 4: 1,5,8,4; 5: 1,4,3,2; 6: 5,6,7,8
        return [
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
            [0, 3, 2, 1],
            [4, 5, 6, 7],
        ]

    def s2n_tet4():
        # Exodus TETRA faces:
        # 1: 1,2,4; 2: 2,3,4; 3: 1,4,3; 4: 1,3,2
        return [
            [0, 1, 3],
            [1, 2, 3],
            [0, 3, 2],
            [0, 2, 1],
        ]

    def s2n_tri3():
        # Exodus TRI edges for 2D; used here as generic 2-node/edge
        return [
            [0, 1],
            [1, 2],
            [2, 0],
        ]

    def s2n_pyramid5():
        # Exodus PYRAMID faces:
        # 1: 1,2,5; 2: 2,3,5; 3: 3,4,5; 4: 4,1,5; 5: 1,4,3,2
        return [
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
            [0, 3, 2, 1],
        ]

    ss_to_nodelist = {}
    ss_to_nodelist["QUAD4"] = s2n_quad4()

    ss_to_nodelist["HEX8"] = s2n_hex8()
    ss_to_nodelist["HEX"] = s2n_hex8()

    ss_to_nodelist["TET4"] = s2n_tet4()
    ss_to_nodelist["TETRA"] = s2n_tet4()
    ss_to_nodelist["tetra"] = s2n_tet4()
    ss_to_nodelist["tetra4"] = s2n_tet4()

    ss_to_nodelist["TRI3"] = s2n_tri3()

    ss_to_nodelist["PYRAMID5"] = s2n_pyramid5()
    ss_to_nodelist["PYRAMID"] = s2n_pyramid5()

    def normalize_elem_type(elem_type):
        if elem_type is None:
            return None
        et = elem_type.upper()
        if et.startswith("HEX"):
            return "HEX8"
        if et.startswith("TET"):
            return "TET4"
        if et.startswith("PYRAMID"):
            return "PYRAMID5"
        if et.startswith("QUAD"):
            return "QUAD4"
        if et.startswith("TRI"):
            return "TRI3"
        return et

    #########################################
    # Elements (now supports mixed types)
    #########################################
    num_elem = nc.dimensions["num_elem"].size
    print(f"num_elem = {num_elem}")

    num_el_blk = nc.dimensions["num_el_blk"].size
    print(f"num_el_blk = {num_el_blk}")

    # Determine max nodes per element over all blocks
    num_nod_per_el_max = 0
    elem_type_per_block = []
    nodes_per_block = []

    for b in range(num_el_blk):
        var_name = f"num_nod_per_el{b+1}"
        if var_name in nc.dimensions:
            num_nod_per_el = nc.dimensions[var_name].size
        else:
            # Fallback: infer from connectivity shape
            connect_b = nc.variables[f"connect{b+1}"]
            _, num_nod_per_el = connect_b.shape

        connect_b = nc.variables[f"connect{b+1}"]
        elem_type_b = getattr(connect_b, "elem_type", None)

        print(f"block {b+1}: elem_type = {elem_type_b}, num_nod_per_el = {num_nod_per_el}")

        elem_type_per_block.append(elem_type_b)
        nodes_per_block.append(num_nod_per_el)
        num_nod_per_el_max = max(num_nod_per_el_max, num_nod_per_el)

    # Global connectivity (padded) and per-element type
    connect = np.zeros((num_elem, num_nod_per_el_max), dtype=idx_t)
    elem_type_arr = np.empty(num_elem, dtype=object)
    nodes_per_elem = np.zeros(num_elem, dtype=np.int16)

    offset = 0
    for b in range(num_el_blk):
        connect_b = nc.variables[f"connect{b+1}"]
        nelements, nnodesxelem = connect_b.shape

        block_begin = offset
        block_end = offset + nelements

        # Optional block name from eb_prop*
        name = None
        try:
            eb_prop_b = nc.variables[f"eb_prop{b+2}"]
            if eb_prop_b is not None:
                name = eb_prop_b.__dict__.get("name", None)
                if name is not None:
                    name = name.lower()
        except KeyError:
            eb_prop_b = None

        if name is not None:
            np.array([block_begin, block_end], dtype=np.int64).tofile(
                f"{output_folder}/blocks/{name}.int64.raw"
            )

        # Fill connectivity (pad remaining entries with 0 -> becomes -1 after shift)
        connect[block_begin:block_end, :nnodesxelem] = connect_b[:].astype(idx_t)
        elem_type_arr[block_begin:block_end] = elem_type_per_block[b]
        nodes_per_elem[block_begin:block_end] = nnodesxelem

        offset += nelements

    # Write global connectivity columns; padded entries become -1
    for i in range(0, num_nod_per_el_max):
        ii = connect[:, i].astype(idx_t) - 1
        ii.tofile(f"{output_folder}/i{i}.raw")

    # Also dump nodes_per_elem so a consumer knows how many nodes each element really has
    nodes_per_elem.tofile(f"{output_folder}/nodes_per_elem.int16.raw")

    # Precompute normalized element types for fast lookup in sidesets
    elem_type_norm = np.empty(num_elem, dtype=object)
    for i in range(num_elem):
        elem_type_norm[i] = normalize_elem_type(elem_type_arr[i])

    #########################################
    # Sidesets (mixed topologies supported)
    #########################################
    num_sidesets = 0

    if "num_side_sets" in nc.dimensions:
        num_sidesets = nc.dimensions["num_side_sets"].size
    else:
        return

    print(f"num_sidesets={num_sidesets}")

    ss_names = nc.variables["ss_names"]

    sideset_dir = f"{output_folder}/sidesets"
    if num_sidesets > 0:
        mkdir(sideset_dir)

    for i in range(0, num_sidesets):
        ssidx = i + 1

        name = netCDF4.chartostring(ss_names[i])

        if name == "":
            name = f"sideset{ssidx}"

        print(f"sideset = {name}")

        key = f"elem_ss{ssidx}"
        e_ss = nc.variables[key]

        key = f"side_ss{ssidx}"
        s_ss = nc.variables[key]

        this_sideset_dir = f"{sideset_dir}/{name}"
        mkdir(this_sideset_dir)

        nsides = len(e_ss[:])

        parent = np.zeros(nsides, dtype=element_idx_t)
        local_face_idx = np.zeros(nsides, dtype=np.int16)

        # First pass: collect nodes per side (variable length) and track max
        side_nodes = []
        max_nodes_per_side = 0

        for n in range(nsides):
            e = e_ss[n] - 1
            s = s_ss[n] - 1

            parent[n] = e
            local_face_idx[n] = s

            etype = elem_type_norm[e]
            if etype not in ss_to_nodelist:
                raise KeyError(f"No side-to-node map for element type '{etype}' in sideset '{name}'")

            s2n_map = ss_to_nodelist[etype]

            if s < 0 or s >= len(s2n_map):
                raise IndexError(f"Side index {s} out of range for element type '{etype}'")

            lnodes = s2n_map[s]

            nodes_this_side = []
            for ln in lnodes:
                node = connect[e, ln] - 1
                nodes_this_side.append(node)

            side_nodes.append(nodes_this_side)
            max_nodes_per_side = max(max_nodes_per_side, len(nodes_this_side))

        # Write local_face_idx and parent
        local_face_idx.tofile(f"{this_sideset_dir}/lfi.int16.raw")
        parent.tofile(f"{this_sideset_dir}/parent.raw")

        # Second pass: build columns name.0.raw, name.1.raw, ...
        # If some sides have fewer nodes (e.g. pyramid tri faces vs quad base),
        # we pad with -1 as a sentinel.
        for d in range(max_nodes_per_side):
            col = []
            for nodes_this_side in side_nodes:
                if d < len(nodes_this_side):
                    col.append(nodes_this_side[d])
                else:
                    col.append(-1)
            path = f"{this_sideset_dir}/{name}.{d}.raw"
            ii = np.array(col, dtype=idx_t)
            ii.tofile(path)


if __name__ == "__main__":

    usage = f"usage: {sys.argv[0]} <input_mesh> <output_folder>"

    if len(sys.argv) < 3:
        print(usage)
        exit()

    input_mesh = sys.argv[1]
    output_folder = sys.argv[2]

    try:
        opts, args = getopt.getopt(sys.argv[3:], "h", ["help"])

    except getopt.GetoptError as err:
        print(err)
        print(usage)
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(usage)
            sys.exit()

    exodusII_to_raw(input_mesh, output_folder)
