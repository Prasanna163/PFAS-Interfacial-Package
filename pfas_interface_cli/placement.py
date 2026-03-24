from pathlib import Path

from .geometry import (
    centroid,
    norm,
    read_xyz,
    rotate_atoms_about_point,
    rotation_matrix_from_vectors,
    write_xyz,
    xy_center,
    xyz_bounds,
)


HEAD_GROUP_ELEMENTS = {"O", "N", "S", "P"}
ORIENTATION_TARGET_AXES = {
    "perpendicular": (0.0, 0.0, 1.0),
    "parallel": (1.0, 0.0, 0.0),
    "as-is": None,
}


def infer_head_tail_axis(pfas_atoms):
    heavy_atoms = [atom for atom in pfas_atoms if atom["element"] != "H"]
    candidate_atoms = heavy_atoms if heavy_atoms else pfas_atoms
    head_atoms = [atom for atom in candidate_atoms if atom["element"] in HEAD_GROUP_ELEMENTS]
    if not head_atoms:
        head_atoms = candidate_atoms

    head_center = centroid(head_atoms)

    def squared_distance(atom):
        return (
            (atom["x"] - head_center[0]) ** 2
            + (atom["y"] - head_center[1]) ** 2
            + (atom["z"] - head_center[2]) ** 2
        )

    tail_anchors = sorted(candidate_atoms, key=squared_distance, reverse=True)[:4]
    tail_center = centroid(tail_anchors)
    axis = (
        tail_center[0] - head_center[0],
        tail_center[1] - head_center[1],
        tail_center[2] - head_center[2],
    )
    if norm(axis) < 1e-8:
        axis = (0.0, 0.0, 1.0)

    return {
        "axis_vector": axis,
        "head_centroid": head_center,
        "tail_centroid": tail_center,
        "head_atom_indices": [atom["index"] for atom in head_atoms],
        "tail_anchor_indices": [atom["index"] for atom in tail_anchors],
    }


def orient_pfas_atoms(pfas_atoms, orientation):
    if orientation not in ORIENTATION_TARGET_AXES:
        raise ValueError(
            f"Unsupported orientation '{orientation}'. Expected one of: {', '.join(sorted(ORIENTATION_TARGET_AXES))}."
        )

    inference = infer_head_tail_axis(pfas_atoms)
    target_axis = ORIENTATION_TARGET_AXES[orientation]
    if target_axis is None:
        oriented_atoms = [
            {
                "index": atom["index"],
                "element": atom["element"],
                "x": atom["x"],
                "y": atom["y"],
                "z": atom["z"],
            }
            for atom in pfas_atoms
        ]
        rotation_matrix = (
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        )
    else:
        rotation_matrix = rotation_matrix_from_vectors(inference["axis_vector"], target_axis)
        oriented_atoms = rotate_atoms_about_point(pfas_atoms, rotation_matrix, inference["head_centroid"])

    return {
        "atoms": oriented_atoms,
        "orientation": orientation,
        "axis_inference": inference,
        "target_axis": target_axis,
        "rotation_matrix": rotation_matrix,
    }


def merge_pfas_above_slab(
    slab_xyz,
    pfas_xyz,
    output_xyz,
    gap,
    x_shift=0.0,
    y_shift=0.0,
    orientation="as-is",
):
    slab_comment, slab_atoms = read_xyz(slab_xyz)
    pfas_comment, pfas_atoms = read_xyz(pfas_xyz)
    orientation_payload = orient_pfas_atoms(pfas_atoms, orientation)
    oriented_pfas_atoms = orientation_payload["atoms"]

    slab_bounds = xyz_bounds(slab_atoms)
    pfas_bounds = xyz_bounds(oriented_pfas_atoms)
    slab_center_x, slab_center_y = xy_center(slab_atoms)
    if orientation == "perpendicular":
        head_index_set = set(orientation_payload["axis_inference"]["head_atom_indices"])
        ref_atoms = [atom for atom in oriented_pfas_atoms if atom["index"] in head_index_set]
        if not ref_atoms:
            ref_atoms = oriented_pfas_atoms
    else:
        ref_atoms = oriented_pfas_atoms

    pfas_center_x, pfas_center_y = xy_center(ref_atoms)
    pfas_reference_z_min = min(atom["z"] for atom in ref_atoms)

    dx = slab_center_x - pfas_center_x + x_shift
    dy = slab_center_y - pfas_center_y + y_shift
    dz = (slab_bounds["z_max"] + gap) - pfas_reference_z_min

    placed_pfas_atoms = []
    for atom in oriented_pfas_atoms:
        placed_pfas_atoms.append(
            {
                "index": len(placed_pfas_atoms) + 1,
                "element": atom["element"],
                "x": atom["x"] + dx,
                "y": atom["y"] + dy,
                "z": atom["z"] + dz,
            }
        )

    merged_atoms = []
    for atom in slab_atoms:
        merged_atoms.append(
            {
                "index": len(merged_atoms) + 1,
                "element": atom["element"],
                "x": atom["x"],
                "y": atom["y"],
                "z": atom["z"],
            }
        )
    for atom in placed_pfas_atoms:
        merged_atoms.append(
            {
                "index": len(merged_atoms) + 1,
                "element": atom["element"],
                "x": atom["x"],
                "y": atom["y"],
                "z": atom["z"],
            }
        )

    comment = (
        f"Slab: {Path(slab_xyz).name}; PFAS: {Path(pfas_xyz).name}; "
        f"orientation={orientation}; gap={gap:.3f} A; slab_atoms=1-{len(slab_atoms)}; "
        f"pfas_atoms={len(slab_atoms)+1}-{len(merged_atoms)}"
    )
    write_xyz(output_xyz, comment, merged_atoms)

    return {
        "output_xyz": str(output_xyz),
        "slab_atom_range": [1, len(slab_atoms)],
        "pfas_atom_range": [len(slab_atoms) + 1, len(merged_atoms)],
        "translation": {"dx": dx, "dy": dy, "dz": dz},
        "slab_bounds": slab_bounds,
        "pfas_bounds_before": pfas_bounds,
        "pfas_bounds_after": xyz_bounds(placed_pfas_atoms),
        "slab_comment": slab_comment,
        "pfas_comment": pfas_comment,
        "orientation": orientation,
        "placement_reference": {
            "xy_reference": "head_group_bbox_center" if orientation == "perpendicular" else "whole_molecule_bbox_center",
            "z_reference": "head_group_zmin" if orientation == "perpendicular" else "whole_molecule_zmin",
        },
        "orientation_transform": {
            "target_axis": orientation_payload["target_axis"],
            "axis_inference": orientation_payload["axis_inference"],
            "rotation_matrix": orientation_payload["rotation_matrix"],
        },
    }


__all__ = [
    "HEAD_GROUP_ELEMENTS",
    "ORIENTATION_TARGET_AXES",
    "infer_head_tail_axis",
    "merge_pfas_above_slab",
    "orient_pfas_atoms",
]
