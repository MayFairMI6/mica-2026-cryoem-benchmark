#!/usr/bin/env python3
"""Lossless, reference-independent MRC map frame canonicalization.

Requires numpy and mrcfile. Density samples are only permuted and optionally
zero-padded, never interpolated, sharpened, masked, or fitted to a structure.
Output NumPy order is Z,Y,X; MAPC/MAPR/MAPS = 1,2,3; starts = 0; ORIGIN
stores the first voxel's physical XYZ position in angstroms.

MRC2014 has two spatial origin representations. The default strict policy
accepts one nonzero representation, or two agreeing representations. Conflicting
nonzero ORIGIN and NSTART require an explicit, externally justified choice.
See https://www.ccpem.ac.uk/mrc-format/mrc2014/ and /mrc-format-proposals/.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import mrcfile
import numpy as np


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def data_sha256(data):
    # Preserve the original numerical dtype; use little endian canonical bytes.
    dtype = data.dtype.newbyteorder("<")
    return hashlib.sha256(np.asarray(data, dtype=dtype).tobytes(order="C")).hexdigest()


def spatial_header(mrc):
    axes_crs = np.array([int(mrc.header[k]) for k in ("mapc", "mapr", "maps")])
    if sorted(axes_crs.tolist()) != [1, 2, 3]:
        raise ValueError("MAPC/MAPR/MAPS must be a permutation of 1,2,3")
    if mrc.data.ndim != 3 or np.iscomplexobj(mrc.data):
        raise ValueError("Only real-valued 3D density maps are supported")
    angles = np.array([float(mrc.header.cellb[k]) for k in ("alpha", "beta", "gamma")])
    if not np.allclose(angles, 90.0, atol=1e-3):
        raise ValueError("Nonorthogonal cell angles need a crystallographic basis adapter")
    voxel = np.array([float(mrc.voxel_size[k]) for k in ("x", "y", "z")])
    if np.any(~np.isfinite(voxel)) or np.any(voxel <= 0):
        raise ValueError("Positive finite XYZ voxel sizes required")
    starts_crs = np.array([int(mrc.header[k]) for k in ("nxstart", "nystart", "nzstart")])
    starts_xyz = np.zeros(3, dtype=int)
    for axis, start in zip(axes_crs, starts_crs):
        starts_xyz[axis - 1] = start
    origin = np.array([float(mrc.header.origin[k]) for k in ("x", "y", "z")])
    if np.any(~np.isfinite(origin)):
        raise ValueError("Finite ORIGIN fields required")
    return axes_crs, starts_crs, starts_xyz, voxel, origin


def first_voxel_xyz(origin, starts_xyz, voxel, policy="strict"):
    start_origin = starts_xyz * voxel
    have_origin = not np.allclose(origin, 0.0, atol=1e-7)
    have_starts = bool(np.any(starts_xyz))
    if policy == "origin":
        return origin.copy(), "explicit_ORIGIN"
    if policy == "start":
        return start_origin, "explicit_NSTART"
    if policy == "additive":
        return origin + start_origin, "explicit_additive_ORIGIN_plus_NSTART"
    if policy != "strict":
        raise ValueError("Unknown origin policy")
    if have_origin and have_starts:
        if np.allclose(origin, start_origin, atol=2e-4, rtol=0):
            return origin.copy(), "agreeing_ORIGIN_and_NSTART"
        raise ValueError("Conflicting nonzero ORIGIN and NSTART: establish source convention before choosing --origin-policy")
    if have_origin:
        return origin.copy(), "ORIGIN_with_zero_starts"
    return start_origin, "NSTART_with_zero_ORIGIN"


def canonicalize(source, destination, *, cube_min=None, origin_policy="strict", world_tolerance_A=2e-4):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination:
        raise ValueError("Write to a different path; deposited maps are immutable")
    if destination.exists():
        raise FileExistsError(destination)
    started = time.perf_counter()
    with mrcfile.open(source, permissive=False) as native:
        axes, starts_crs, starts_xyz, voxel, origin = spatial_header(native)
        physical_first, policy_used = first_voxel_xyz(origin, starts_xyz, voxel, origin_policy)
        input_shape = list(native.data.shape)
        input_dtype = str(native.data.dtype)
        native_hash = data_sha256(native.data)
        native_physical_axes = (axes[::-1] - 1).tolist()  # native data S,R,C
        transpose_order = [native_physical_axes.index(axis) for axis in (2, 1, 0)]
        canonical = np.ascontiguousarray(np.transpose(native.data, transpose_order))
        canonical_hash = data_sha256(canonical)
        shape = np.array(canonical.shape, dtype=int)
        before = np.zeros(3, dtype=int)
        after = np.zeros(3, dtype=int)
        if cube_min is not None:
            side = max(int(cube_min), int(shape.max()))
            side += side % 2
            before = side // 2 - shape // 2  # matches ModelAngelo's centering
            after = side - shape - before
        data = np.pad(canonical, list(zip(before, after)), mode="constant", constant_values=0)
        output_origin = physical_first - before[::-1] * voxel
        destination.parent.mkdir(parents=True, exist_ok=True)
        with mrcfile.new(destination, overwrite=False) as out:
            out.set_data(data)
            out.voxel_size = tuple(voxel)
            out.header.mapc, out.header.mapr, out.header.maps = 1, 2, 3
            out.nstart = (0, 0, 0)
            out.header.origin = tuple(output_origin)
            out.header.cellb = (90.0, 90.0, 90.0)
            out.header.ispg = 1
            out.header.nversion = 20141
            out.add_label("MICA benchmark: lossless axes/origin canonicalization; see transform JSON")
            out.update_header_stats()
        # Reopen the actual serialized output: float32 header rounding is tested.
        with mrcfile.open(destination, permissive=False) as checked:
            _, _, _, checked_voxel, checked_origin = spatial_header(checked)
            region = checked.data[tuple(slice(int(p), int(p + n)) for p, n in zip(before, shape))]
            same_values = bool(np.array_equal(region, canonical, equal_nan=True))
            region_hash = data_sha256(region)
            points_zyx = np.array(list(itertools.product(*[[0, int(n // 2), int(n - 1)] for n in shape])), dtype=float)
            native_world = physical_first + points_zyx[:, ::-1] * voxel
            output_world = checked_origin + (points_zyx + before)[:, ::-1] * checked_voxel
            world_errors = np.linalg.norm(native_world - output_world, axis=1)
            max_error = float(world_errors.max())
            shape_out = list(checked.data.shape)
            zero_starts = all(int(checked.header[k]) == 0 for k in ("nxstart", "nystart", "nzstart"))
        if not same_values or canonical_hash != region_hash or max_error > world_tolerance_A or not zero_starts:
            raise AssertionError("Density preservation or physical-coordinate equivalence failed")
        inverse = np.argsort(transpose_order)
        restored_native = np.transpose(canonical, inverse)
        roundtrip_hash = data_sha256(restored_native)
        if roundtrip_hash != native_hash:
            raise AssertionError("Inverse axis permutation changed voxel data")
    manifest = {
        "source": str(source), "destination": str(destination),
        "source_sha256": file_sha256(source), "destination_sha256": file_sha256(destination),
        "source_density_sha256": native_hash, "restored_source_density_sha256": roundtrip_hash,
        "canonical_unpadded_density_sha256": canonical_hash, "output_preserved_region_sha256": region_hash,
        "source_dtype": input_dtype, "source_shape_src": input_shape,
        "source_axes_crs": axes.tolist(), "source_starts_crs": starts_crs.tolist(),
        "source_starts_xyz": starts_xyz.tolist(), "source_origin_xyz_A": origin.tolist(),
        "origin_policy_requested": origin_policy, "origin_policy_used": policy_used,
        "source_first_voxel_xyz_A": physical_first.tolist(), "voxel_xyz_A": voxel.tolist(),
        "transpose_src_to_zyx": transpose_order, "output_shape_zyx": shape_out,
        "padding_before_zyx": before.tolist(), "padding_after_zyx": after.tolist(),
        "output_origin_xyz_A_requested": output_origin.tolist(),
        "output_origin_xyz_A_serialized": checked_origin.tolist(),
        "output_voxel_xyz_A_serialized": checked_voxel.tolist(), "output_starts_xyz": [0, 0, 0],
        "reference_structure_used": False, "density_interpolated": False,
        "density_region_bitwise_preserved": canonical_hash == region_hash,
        "inverse_permutation_hash_matches_native": roundtrip_hash == native_hash,
        "world_coordinate_landmarks_tested": len(points_zyx), "max_world_position_error_A": max_error,
        "world_tolerance_A": world_tolerance_A, "elapsed_s": time.perf_counter() - started,
        "predictor_frame_caution": "Some predictors ignore ORIGIN. Use the recorded first-voxel translation to map predictor-local coordinates back to this physical frame; verify predictor behavior before scoring."
    }
    destination.with_suffix(destination.suffix + ".transform.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--cube-min", type=int)
    parser.add_argument("--origin-policy", choices=["strict", "origin", "start", "additive"], default="strict")
    args = parser.parse_args()
    print(json.dumps(canonicalize(args.source, args.destination, cube_min=args.cube_min, origin_policy=args.origin_policy), indent=2))


if __name__ == "__main__":
    main()
