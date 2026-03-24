import argparse
import json
import math
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from knf_core.pipeline import KNFPipeline


ENERGY_RE = re.compile(r"energy:\s*([-+]?\d+(?:\.\d+)?)", re.IGNORECASE)
HEAD_GROUP_ELEMENTS = {"O", "N", "S", "P"}
ORIENTATION_TARGET_AXES = {
    "perpendicular": (0.0, 0.0, 1.0),
    "parallel": (1.0, 0.0, 0.0),
    "as-is": None,
}
KNF_METRIC_NAMES = ["SNCI", "SCDI", "SCDI_variance"] + [f"f{i}" for i in range(1, 10)]

# Production defaults (intentionally fixed to keep CLI surface minimal).
DEFAULT_RUN_ROOT = "pipeline_runs"
DEFAULT_ORIENTATION_SEQUENCE = ["perpendicular", "parallel"]
DEFAULT_GFN_LEVEL = 2
DEFAULT_PFAS_MULTIPLICITY = 1
DEFAULT_GAP_ANGSTROM = 2.5
DEFAULT_X_SHIFT_ANGSTROM = 0.0
DEFAULT_Y_SHIFT_ANGSTROM = 0.0
DEFAULT_CONTACT_ELEMENTS = ["O", "F"]
DEFAULT_KEEP_KNF_INTERMEDIATES = False
DEFAULT_KNF_SCDI_VAR_MIN = None
DEFAULT_KNF_SCDI_VAR_MAX = None


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path} is not a valid XYZ file.")

    n_atoms = int(lines[0].strip())
    comment = lines[1].rstrip("\n")
    atom_lines = lines[2 : 2 + n_atoms]
    if len(atom_lines) != n_atoms:
        raise ValueError(f"{path} declares {n_atoms} atoms but contains {len(atom_lines)} atom lines.")

    atoms = []
    for index, line in enumerate(atom_lines, start=1):
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ line {index} in {path}: {line}")
        atoms.append(
            {
                "index": index,
                "element": parts[0],
                "x": float(parts[1]),
                "y": float(parts[2]),
                "z": float(parts[3]),
            }
        )
    return comment, atoms


def write_xyz(path, comment, atoms):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(f"{comment}\n")
        for atom in atoms:
            handle.write(
                f"{atom['element']:2s} {atom['x']: .6f} {atom['y']: .6f} {atom['z']: .6f}\n"
            )


def xyz_bounds(atoms):
    xs = [atom["x"] for atom in atoms]
    ys = [atom["y"] for atom in atoms]
    zs = [atom["z"] for atom in atoms]
    return {
        "x_min": min(xs),
        "x_max": max(xs),
        "y_min": min(ys),
        "y_max": max(ys),
        "z_min": min(zs),
        "z_max": max(zs),
    }


def xy_center(atoms):
    bounds = xyz_bounds(atoms)
    return (
        0.5 * (bounds["x_min"] + bounds["x_max"]),
        0.5 * (bounds["y_min"] + bounds["y_max"]),
    )


def centroid(atoms):
    if not atoms:
        raise ValueError("Cannot compute centroid of an empty atom list.")
    count = float(len(atoms))
    return (
        sum(atom["x"] for atom in atoms) / count,
        sum(atom["y"] for atom in atoms) / count,
        sum(atom["z"] for atom in atoms) / count,
    )


def dot(vec_a, vec_b):
    return vec_a[0] * vec_b[0] + vec_a[1] * vec_b[1] + vec_a[2] * vec_b[2]


def cross(vec_a, vec_b):
    return (
        vec_a[1] * vec_b[2] - vec_a[2] * vec_b[1],
        vec_a[2] * vec_b[0] - vec_a[0] * vec_b[2],
        vec_a[0] * vec_b[1] - vec_a[1] * vec_b[0],
    )


def norm(vec):
    return math.sqrt(dot(vec, vec))


def normalize(vec):
    magnitude = norm(vec)
    if magnitude < 1e-12:
        raise ValueError("Cannot normalize a zero-length vector.")
    return (vec[0] / magnitude, vec[1] / magnitude, vec[2] / magnitude)


def mat_vec_mul(matrix, vec):
    return (
        matrix[0][0] * vec[0] + matrix[0][1] * vec[1] + matrix[0][2] * vec[2],
        matrix[1][0] * vec[0] + matrix[1][1] * vec[1] + matrix[1][2] * vec[2],
        matrix[2][0] * vec[0] + matrix[2][1] * vec[1] + matrix[2][2] * vec[2],
    )


def identity_matrix():
    return (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )


def rotation_matrix_from_vectors(source_vec, target_vec):
    source = normalize(source_vec)
    target = normalize(target_vec)
    cosine = dot(source, target)

    if cosine > 1.0 - 1e-12:
        return identity_matrix()

    if cosine < -1.0 + 1e-12:
        # 180 degree flip: choose any stable axis orthogonal to source.
        trial = (1.0, 0.0, 0.0) if abs(source[0]) < 0.9 else (0.0, 1.0, 0.0)
        axis = normalize(cross(source, trial))
        x, y, z = axis
        # Rotation by pi: R = -I + 2 * (axis * axis^T)
        return (
            (-1.0 + 2.0 * x * x, 2.0 * x * y, 2.0 * x * z),
            (2.0 * y * x, -1.0 + 2.0 * y * y, 2.0 * y * z),
            (2.0 * z * x, 2.0 * z * y, -1.0 + 2.0 * z * z),
        )

    v = cross(source, target)
    s = norm(v)
    vx = (
        (0.0, -v[2], v[1]),
        (v[2], 0.0, -v[0]),
        (-v[1], v[0], 0.0),
    )
    vx2 = (
        (
            vx[0][0] * vx[0][0] + vx[0][1] * vx[1][0] + vx[0][2] * vx[2][0],
            vx[0][0] * vx[0][1] + vx[0][1] * vx[1][1] + vx[0][2] * vx[2][1],
            vx[0][0] * vx[0][2] + vx[0][1] * vx[1][2] + vx[0][2] * vx[2][2],
        ),
        (
            vx[1][0] * vx[0][0] + vx[1][1] * vx[1][0] + vx[1][2] * vx[2][0],
            vx[1][0] * vx[0][1] + vx[1][1] * vx[1][1] + vx[1][2] * vx[2][1],
            vx[1][0] * vx[0][2] + vx[1][1] * vx[1][2] + vx[1][2] * vx[2][2],
        ),
        (
            vx[2][0] * vx[0][0] + vx[2][1] * vx[1][0] + vx[2][2] * vx[2][0],
            vx[2][0] * vx[0][1] + vx[2][1] * vx[1][1] + vx[2][2] * vx[2][1],
            vx[2][0] * vx[0][2] + vx[2][1] * vx[1][2] + vx[2][2] * vx[2][2],
        ),
    )
    factor = (1.0 - cosine) / (s * s)

    eye = identity_matrix()
    return (
        (
            eye[0][0] + vx[0][0] + vx2[0][0] * factor,
            eye[0][1] + vx[0][1] + vx2[0][1] * factor,
            eye[0][2] + vx[0][2] + vx2[0][2] * factor,
        ),
        (
            eye[1][0] + vx[1][0] + vx2[1][0] * factor,
            eye[1][1] + vx[1][1] + vx2[1][1] * factor,
            eye[1][2] + vx[1][2] + vx2[1][2] * factor,
        ),
        (
            eye[2][0] + vx[2][0] + vx2[2][0] * factor,
            eye[2][1] + vx[2][1] + vx2[2][1] * factor,
            eye[2][2] + vx[2][2] + vx2[2][2] * factor,
        ),
    )


def rotate_atoms_about_point(atoms, rotation_matrix, origin):
    rotated = []
    ox, oy, oz = origin
    for atom in atoms:
        local = (atom["x"] - ox, atom["y"] - oy, atom["z"] - oz)
        rx, ry, rz = mat_vec_mul(rotation_matrix, local)
        rotated.append(
            {
                "index": atom["index"],
                "element": atom["element"],
                "x": rx + ox,
                "y": ry + oy,
                "z": rz + oz,
            }
        )
    return rotated


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
        rotation_matrix = identity_matrix()
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


def parse_xtb_energy_from_xyz(path):
    comment, _ = read_xyz(path)
    match = ENERGY_RE.search(comment)
    if not match:
        raise ValueError(f"Could not parse xTB energy from comment line in {path}.")
    return float(match.group(1))


def run_command(command, cwd, log_path):
    with open(log_path, "w", encoding="utf-8", errors="replace") as log_handle:
        result = subprocess.run(
            command,
            cwd=str(cwd),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(command)}")


def sort_water_slab_by_oxygen_z(input_xyz, output_xyz):
    comment, atoms = read_xyz(input_xyz)
    if len(atoms) % 3 != 0:
        raise ValueError("Expected a water slab with atom count divisible by 3.")

    waters = []
    for start in range(0, len(atoms), 3):
        triplet = atoms[start : start + 3]
        elements = [atom["element"] for atom in triplet]
        if elements != ["O", "H", "H"]:
            raise ValueError(
                f"Expected O-H-H ordering in water slab, found {elements} at atoms {start + 1}-{start + 3}."
            )
        waters.append((triplet[0]["z"], triplet))

    waters.sort(key=lambda item: item[0])
    sorted_atoms = []
    for _, triplet in waters:
        sorted_atoms.extend(
            {
                "index": len(sorted_atoms) + 1,
                "element": atom["element"],
                "x": atom["x"],
                "y": atom["y"],
                "z": atom["z"],
            }
            for atom in triplet
        )

    write_xyz(output_xyz, comment, sorted_atoms)
    return len(waters)


def optimize_with_xtb(input_xyz, work_dir, charge=0, multiplicity=1, gfn=2):
    work_dir.mkdir(parents=True, exist_ok=True)
    local_input = work_dir / input_xyz.name
    if input_xyz.resolve() != local_input.resolve():
        shutil.copy2(input_xyz, local_input)

    command = [
        "xtb",
        local_input.name,
        "--opt",
        "tight",
        "--gfn",
        str(gfn),
        "--charge",
        str(charge),
        "--uhf",
        str(multiplicity - 1),
    ]
    log_path = work_dir / "xtb_run.log"
    run_command(command, work_dir, log_path)

    optimized_xyz = work_dir / "xtbopt.xyz"
    if not optimized_xyz.exists():
        raise FileNotFoundError(f"xTB did not produce {optimized_xyz}.")

    return {
        "work_dir": str(work_dir),
        "input_xyz": str(local_input),
        "optimized_xyz": str(optimized_xyz),
        "log": str(log_path),
        "energy_eh": parse_xtb_energy_from_xyz(optimized_xyz),
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


def distance(atom_a, atom_b):
    return math.sqrt(
        (atom_a["x"] - atom_b["x"]) ** 2
        + (atom_a["y"] - atom_b["y"]) ** 2
        + (atom_a["z"] - atom_b["z"]) ** 2
    )


def analyze_interface_geometry(interface_xyz, slab_atom_count, contact_elements):
    _, atoms = read_xyz(interface_xyz)
    slab_atoms = atoms[:slab_atom_count]
    pfas_atoms = atoms[slab_atom_count:]

    slab_bounds = xyz_bounds(slab_atoms)
    pfas_bounds = xyz_bounds(pfas_atoms)
    total_bounds = xyz_bounds(atoms)

    water_h_atoms = [atom for atom in slab_atoms if atom["element"] == "H"]
    pfas_contact_atoms = [atom for atom in pfas_atoms if atom["element"] in contact_elements]
    if not pfas_contact_atoms:
        pfas_contact_atoms = [atom for atom in pfas_atoms if atom["element"] not in {"C", "H"}]

    closest_pair = None
    for pfas_atom in pfas_contact_atoms:
        for water_atom in water_h_atoms:
            current_distance = distance(pfas_atom, water_atom)
            if closest_pair is None or current_distance < closest_pair["distance_angstrom"]:
                closest_pair = {
                    "distance_angstrom": current_distance,
                    "pfas_atom": pfas_atom,
                    "water_atom": water_atom,
                }

    min_gap = pfas_bounds["z_min"] - slab_bounds["z_max"]
    max_gap = pfas_bounds["z_max"] - slab_bounds["z_max"]

    if closest_pair is None:
        contact_classification = "No PFAS contact atoms found for the requested element filter."
    else:
        distance_value = closest_pair["distance_angstrom"]
        if distance_value <= 2.2:
            contact_classification = "Hydrogen-bond-like / strong specific contact"
        elif distance_value <= 3.5:
            contact_classification = "Weak interfacial contact"
        else:
            contact_classification = "Hovering / no direct specific contact"

    if pfas_bounds["z_min"] >= slab_bounds["z_max"]:
        interface_classification = "Surface-localized above slab top"
    elif pfas_bounds["z_max"] <= slab_bounds["z_max"]:
        interface_classification = "Immersed within slab region"
    else:
        interface_classification = "Partially penetrating interface"

    return {
        "total_bounds": total_bounds,
        "slab_bounds": slab_bounds,
        "pfas_bounds": pfas_bounds,
        "slab_thickness_angstrom": slab_bounds["z_max"] - slab_bounds["z_min"],
        "pfas_offset_from_slab_top": {
            "z_min_minus_slab_top": min_gap,
            "z_max_minus_slab_top": max_gap,
        },
        "interface_classification": interface_classification,
        "contact_elements_used": sorted({atom["element"] for atom in pfas_contact_atoms}),
        "closest_contact": closest_pair,
        "contact_classification": contact_classification,
    }


def run_knf_on_existing_optimized_xyz(
    input_xyz,
    output_root,
    charge,
    multiplicity,
    keep_full_files=False,
    scdi_var_min=None,
    scdi_var_max=None,
):
    pipeline = KNFPipeline(
        input_file=str(input_xyz),
        charge=charge,
        spin=multiplicity,
        water=False,
        force=False,
        clean=False,
        debug=False,
        output_root=str(output_root),
        keep_full_files=keep_full_files,
        nci_backend="torch",
        scdi_var_min=scdi_var_min,
        scdi_var_max=scdi_var_max,
    )
    pipeline.setup_directories()
    seeded_xtbopt = Path(pipeline.results_dir) / "xtbopt.xyz"
    shutil.copy2(input_xyz, seeded_xtbopt)
    context = pipeline.run_pre_nci_stage()
    pipeline.run_post_nci_stage(context)

    knf_json = Path(pipeline.results_dir) / "knf.json"
    output_txt = Path(pipeline.results_dir) / "output.txt"
    if not knf_json.exists():
        raise FileNotFoundError(f"KNF output not found: {knf_json}")

    return {
        "results_dir": str(pipeline.results_dir),
        "knf_json": str(knf_json),
        "output_txt": str(output_txt),
        "data": json.loads(knf_json.read_text(encoding="utf-8")),
    }


def extract_knf_metrics(knf_payload):
    vector = knf_payload.get("KNF_vector") or []
    metrics = {
        "SNCI": knf_payload.get("SNCI"),
        "SCDI": knf_payload.get("SCDI"),
        "SCDI_variance": knf_payload.get("SCDI_variance"),
    }
    for index in range(9):
        metrics[f"f{index + 1}"] = vector[index] if index < len(vector) else None
    return metrics


def metric_delta(current, reference):
    if current is None or reference is None:
        return None
    return current - reference


def format_value(value, decimals=6):
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.{decimals}f}"
    return str(value)


def make_table(headers, rows):
    widths = [len(str(header)) for header in headers]
    normalized_rows = []
    for row in rows:
        normalized = [str(cell) for cell in row]
        normalized_rows.append(normalized)
        for index, cell in enumerate(normalized):
            widths[index] = max(widths[index], len(cell))

    def render_row(row):
        return "| " + " | ".join(cell.ljust(widths[index]) for index, cell in enumerate(row)) + " |"

    separator = "|-" + "-|-".join("-" * width for width in widths) + "-|"
    lines = [render_row(headers), separator]
    lines.extend(render_row(row) for row in normalized_rows)
    return lines


def build_text_report(report):
    interface = report["interface_analysis"]
    closest = interface["closest_contact"]
    pfas_metrics = report["knf"]["pfas"]["metrics"]
    interface_metrics = report["knf"]["interface"]["metrics"]
    metric_names = KNF_METRIC_NAMES
    orientation_summary = report.get("orientation_summary") or {}
    orientation_results = report.get("orientations") or {}
    preferred_orientation = orientation_summary.get("preferred_orientation")
    evaluated_orientations = orientation_summary.get("evaluated") or list(orientation_results.keys())

    lines = [
        "PFAS Interface Workflow Summary",
        "================================",
        "",
        f"Run name: {report['run_name']}",
        f"Run directory: {report['run_directory']}",
        "",
        "Inputs",
        "------",
        f"Slab input: {report['inputs']['slab_xyz']}",
        f"PFAS input: {report['inputs']['pfas_xyz']}",
        "",
        "Settings",
        "--------",
        "",
    ]
    lines.extend(
        make_table(
            ["Setting", "Value"],
            [
                ["GFN level", str(report["settings"]["gfn"])],
                ["PFAS/interface charge", str(report["settings"]["pfas_charge"])],
                ["Spin multiplicity", str(report["settings"]["multiplicity"])],
                ["Orientation mode", str(report["settings"].get("orientation_mode", "single"))],
                ["Preferred orientation", preferred_orientation or "n/a"],
                ["Placement gap", f"{report['settings']['placement_gap_angstrom']:.3f} A"],
            ],
        )
    )
    lines.extend(
        [
            "",
            "Atom Ranges",
            "-----------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Group", "Atom Range"],
            [
                [
                    "Slab atoms in interface",
                    f"{report['atom_ranges']['slab_atoms'][0]}-{report['atom_ranges']['slab_atoms'][1]}",
                ],
                [
                    "PFAS atoms in interface",
                    f"{report['atom_ranges']['pfas_atoms'][0]}-{report['atom_ranges']['pfas_atoms'][1]}",
                ],
            ],
        )
    )
    lines.extend(
        [
            "",
            "Energies (Eh)",
            "-------------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Quantity", "Value (Eh)"],
            [
                ["E_slab", f"{report['energies_eh']['slab']:.12f}"],
                ["E_PFAS", f"{report['energies_eh']['pfas']:.12f}"],
                ["E_interface", f"{report['energies_eh']['interface']:.12f}"],
                ["Delta_E_ads", f"{report['energies_eh']['delta_e_ads']:.12f}"],
            ],
        )
    )
    if orientation_results:
        lines.extend(
            [
                "",
                "Orientation Comparison",
                "----------------------",
                "",
            ]
        )
        orientation_rows = []
        ordered_names = [name for name in evaluated_orientations if name in orientation_results]
        for orientation_name in ordered_names:
            payload = orientation_results[orientation_name]
            orientation_rows.append(
                [
                    orientation_name,
                    f"{payload['energies_eh']['interface']:.12f}",
                    f"{payload['energies_eh']['delta_e_ads']:.12f}",
                    payload["interface_analysis"]["contact_classification"],
                ]
            )
        lines.extend(
            make_table(
                [
                    "Orientation",
                    "E_interface (Eh)",
                    "Delta_E_ads (Eh)",
                    "Contact interpretation",
                ],
                orientation_rows,
            )
        )
        anisotropy = orientation_summary.get("delta_delta_e_ads_parallel_minus_perpendicular")
        if anisotropy is not None:
            lines.extend(
                [
                    "",
                    f"Delta_Delta_E_ads (parallel - perpendicular): {anisotropy:.12f} Eh",
                ]
            )
    lines.extend(
        [
            "",
            "Geometry",
            "--------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Metric", "Value"],
            [
                [
                    "Slab Z range",
                    f"{interface['slab_bounds']['z_min']:.4f} to {interface['slab_bounds']['z_max']:.4f}",
                ],
                ["Slab thickness", f"{interface['slab_thickness_angstrom']:.4f} A"],
                [
                    "Interface total Z range",
                    f"{interface['total_bounds']['z_min']:.4f} to {interface['total_bounds']['z_max']:.4f}",
                ],
                [
                    "PFAS Z range",
                    f"{interface['pfas_bounds']['z_min']:.4f} to {interface['pfas_bounds']['z_max']:.4f}",
                ],
                [
                    "PFAS offset from slab top",
                    (
                        f"{interface['pfas_offset_from_slab_top']['z_min_minus_slab_top']:.4f} "
                        f"to {interface['pfas_offset_from_slab_top']['z_max_minus_slab_top']:.4f} A"
                    ),
                ],
                ["Interface classification", interface["interface_classification"]],
            ],
        )
    )
    lines.extend(
        [
            "",
            "Closest PFAS-to-Water Contact",
            "-----------------------------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Metric", "Value"],
            [
                [
                    "Elements checked on PFAS",
                    ", ".join(interface["contact_elements_used"]) if interface["contact_elements_used"] else "none",
                ],
                ["Minimum distance", "n/a" if closest is None else f"{closest['distance_angstrom']:.4f} A"],
                [
                    "Closest pair",
                    (
                        "n/a"
                        if closest is None
                        else (
                            f"PFAS atom {closest['pfas_atom']['index']} ({closest['pfas_atom']['element']}) "
                            f"to water atom {closest['water_atom']['index']} ({closest['water_atom']['element']})"
                        )
                    ),
                ],
                ["Contact interpretation", interface["contact_classification"]],
            ],
        )
    )
    lines.extend(
        [
            "",
            "KNF Metrics",
            "-----------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Metric", "PFAS", "Interface", "Shift (interface - PFAS)"],
            [
                [
                    name,
                    format_value(pfas_metrics.get(name)),
                    format_value(interface_metrics.get(name)),
                    format_value(report["knf"]["delta_interface_minus_pfas"].get(name)),
                ]
                for name in metric_names
            ],
        )
    )
    lines.extend(
        [
            "",
            "Key Output Files",
            "----------------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Artifact", "Path"],
            [
                ["Sorted slab", report["outputs"]["slab_sorted_xyz"]],
                ["Optimized slab", report["outputs"]["slab_optimized_xyz"]],
                ["Optimized PFAS", report["outputs"]["pfas_optimized_xyz"]],
                ["Merged interface start", report["outputs"]["interface_initial_xyz"]],
                ["Optimized interface", report["outputs"]["interface_optimized_xyz"]],
                ["Summary JSON", report["outputs"]["summary_json"]],
                ["Summary TXT", report["outputs"]["summary_txt"]],
                ["KNF PFAS JSON", report["outputs"]["knf_pfas_json"]],
                ["KNF interface JSON", report["outputs"]["knf_interface_json"]],
            ],
        )
    )
    orientation_outputs = report["outputs"].get("orientation_outputs")
    if orientation_outputs:
        lines.extend(
            [
                "",
                "Orientation Output Files",
                "------------------------",
                "",
            ]
        )
        orientation_output_rows = []
        ordered_names = [name for name in evaluated_orientations if name in orientation_outputs]
        for orientation_name in ordered_names:
            payload = orientation_outputs[orientation_name]
            orientation_output_rows.append(
                [
                    orientation_name,
                    payload["interface_initial_xyz"],
                    payload["interface_optimized_xyz"],
                    payload["knf_interface_json"],
                ]
            )
        lines.extend(
            make_table(
                ["Orientation", "Interface start XYZ", "Optimized interface XYZ", "KNF interface JSON"],
                orientation_output_rows,
            )
        )
    return "\n".join(lines) + "\n"


def build_parser():
    parser = argparse.ArgumentParser(
        description="End-to-end PFAS/water-slab xTB + KNF workflow."
    )
    parser.add_argument("slab_xyz", help="Water slab XYZ input")
    parser.add_argument("pfas_xyz", help="PFAS XYZ input")
    parser.add_argument(
        "--pfas-charge",
        type=int,
        default=0,
        help="Charge used for PFAS-only and interface calculations (default: 0).",
    )
    advanced = parser.add_argument_group("Advanced Options (Second Layer)")
    advanced.add_argument(
        "--run-root",
        default=DEFAULT_RUN_ROOT,
        help=f"Directory where timestamped run folders are created (default: {DEFAULT_RUN_ROOT})",
    )
    advanced.add_argument(
        "--run-name",
        default=None,
        help="Optional fixed run folder name. Default: <pfas_stem>_<timestamp>",
    )
    advanced.add_argument(
        "--gap",
        type=float,
        default=DEFAULT_GAP_ANGSTROM,
        help=f"Initial PFAS zmin placement above slab zmax in angstrom (default: {DEFAULT_GAP_ANGSTROM})",
    )
    advanced.add_argument(
        "--x-shift",
        type=float,
        default=DEFAULT_X_SHIFT_ANGSTROM,
        help=f"Additional lateral PFAS x-shift after centering (default: {DEFAULT_X_SHIFT_ANGSTROM})",
    )
    advanced.add_argument(
        "--y-shift",
        type=float,
        default=DEFAULT_Y_SHIFT_ANGSTROM,
        help=f"Additional lateral PFAS y-shift after centering (default: {DEFAULT_Y_SHIFT_ANGSTROM})",
    )
    advanced.add_argument(
        "--orientation-mode",
        choices=["dual", "perpendicular", "parallel"],
        default="dual",
        help=(
            "Interface orientation workflow to run. "
            "'dual' runs both perpendicular and parallel branches, then reports preferred orientation."
        ),
    )
    advanced.add_argument(
        "--multiplicity",
        type=int,
        default=DEFAULT_PFAS_MULTIPLICITY,
        help=f"Spin multiplicity used for PFAS/interface calculations (default: {DEFAULT_PFAS_MULTIPLICITY})",
    )
    advanced.add_argument(
        "--gfn",
        type=int,
        default=DEFAULT_GFN_LEVEL,
        help=f"GFN level passed to xTB (default: {DEFAULT_GFN_LEVEL})",
    )
    advanced.add_argument(
        "--contact-elements",
        default=",".join(DEFAULT_CONTACT_ELEMENTS),
        help=(
            "Comma-separated PFAS elements to test against slab water H for closest contact "
            f"(default: {','.join(DEFAULT_CONTACT_ELEMENTS)})"
        ),
    )
    advanced.add_argument(
        "--keep-knf-intermediates",
        action="store_true",
        default=DEFAULT_KEEP_KNF_INTERMEDIATES,
        help="Keep large KNF intermediate files instead of storage-efficient cleanup.",
    )
    advanced.add_argument(
        "--knf-scdi-var-min",
        type=float,
        default=DEFAULT_KNF_SCDI_VAR_MIN,
        help="Optional fixed KNF SCDI variance minimum for normalized SCDI output.",
    )
    advanced.add_argument(
        "--knf-scdi-var-max",
        type=float,
        default=DEFAULT_KNF_SCDI_VAR_MAX,
        help="Optional fixed KNF SCDI variance maximum for normalized SCDI output.",
    )
    return parser


def run_workflow(args):
    args = argparse.Namespace(**vars(args))

    slab_input = Path(args.slab_xyz).resolve()
    pfas_input = Path(args.pfas_xyz).resolve()
    if not slab_input.exists():
        raise FileNotFoundError(f"Slab input not found: {slab_input}")
    if not pfas_input.exists():
        raise FileNotFoundError(f"PFAS input not found: {pfas_input}")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_name = args.run_name or f"{pfas_input.stem}_{timestamp}"
    run_root = Path(args.run_root).resolve()
    run_dir = run_root / run_name
    if run_dir.exists():
        raise FileExistsError(f"Run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=False)

    if args.orientation_mode == "dual":
        orientation_sequence = list(DEFAULT_ORIENTATION_SEQUENCE)
    else:
        orientation_sequence = [args.orientation_mode]

    sorted_step = run_dir / "01_slab_sort"
    slab_opt_step = run_dir / "02_slab_opt"
    pfas_opt_step = run_dir / "03_pfas_opt"
    orientation_step_map = {}
    next_step_index = 4
    for orientation_name in orientation_sequence:
        root = run_dir / f"{next_step_index:02d}_{orientation_name}"
        orientation_step_map[orientation_name] = {
            "root": root,
            "merge": root / "01_merge",
            "interface_opt": root / "02_interface_opt",
            "knf_interface": root / "03_knf_interface",
        }
        next_step_index += 1
    knf_pfas_step = run_dir / f"{next_step_index:02d}_knf_pfas"

    for step_dir in [sorted_step, slab_opt_step, pfas_opt_step, knf_pfas_step]:
        step_dir.mkdir(parents=True, exist_ok=True)
    for orientation_paths in orientation_step_map.values():
        for step_dir in orientation_paths.values():
            step_dir.mkdir(parents=True, exist_ok=True)

    slab_sorted_xyz = sorted_step / "slab_sorted.xyz"
    sort_water_slab_by_oxygen_z(slab_input, slab_sorted_xyz)
    slab_opt = optimize_with_xtb(
        slab_sorted_xyz,
        slab_opt_step,
        charge=0,
        multiplicity=1,
        gfn=args.gfn,
    )

    pfas_opt = optimize_with_xtb(
        pfas_input,
        pfas_opt_step,
        charge=args.pfas_charge,
        multiplicity=args.multiplicity,
        gfn=args.gfn,
    )

    contact_elements = [part.strip() for part in args.contact_elements.split(",") if part.strip()]
    if not contact_elements:
        raise ValueError("At least one value must be provided to --contact-elements.")

    knf_pfas = run_knf_on_existing_optimized_xyz(
        Path(pfas_opt["optimized_xyz"]),
        knf_pfas_step,
        charge=args.pfas_charge,
        multiplicity=args.multiplicity,
        keep_full_files=args.keep_knf_intermediates,
        scdi_var_min=args.knf_scdi_var_min,
        scdi_var_max=args.knf_scdi_var_max,
    )
    knf_pfas_metrics = extract_knf_metrics(knf_pfas["data"])
    orientation_results = {}
    for orientation_name in orientation_sequence:
        orientation_paths = orientation_step_map[orientation_name]

        interface_initial_xyz = orientation_paths["merge"] / f"{pfas_input.stem}_{orientation_name}_on_slab.xyz"
        merge_info = merge_pfas_above_slab(
            Path(slab_opt["optimized_xyz"]),
            Path(pfas_opt["optimized_xyz"]),
            interface_initial_xyz,
            gap=args.gap,
            x_shift=args.x_shift,
            y_shift=args.y_shift,
            orientation=orientation_name,
        )

        _, slab_atom_end = merge_info["slab_atom_range"]
        interface_opt = optimize_with_xtb(
            interface_initial_xyz,
            orientation_paths["interface_opt"],
            charge=args.pfas_charge,
            multiplicity=args.multiplicity,
            gfn=args.gfn,
        )

        interface_analysis = analyze_interface_geometry(
            Path(interface_opt["optimized_xyz"]),
            slab_atom_count=slab_atom_end,
            contact_elements=contact_elements,
        )

        knf_interface = run_knf_on_existing_optimized_xyz(
            Path(interface_opt["optimized_xyz"]),
            orientation_paths["knf_interface"],
            charge=args.pfas_charge,
            multiplicity=args.multiplicity,
            keep_full_files=args.keep_knf_intermediates,
            scdi_var_min=args.knf_scdi_var_min,
            scdi_var_max=args.knf_scdi_var_max,
        )
        knf_interface_metrics = extract_knf_metrics(knf_interface["data"])
        knf_delta = {
            name: metric_delta(knf_interface_metrics.get(name), knf_pfas_metrics.get(name))
            for name in KNF_METRIC_NAMES
        }

        energies = {
            "slab": slab_opt["energy_eh"],
            "pfas": pfas_opt["energy_eh"],
            "interface": interface_opt["energy_eh"],
        }
        energies["delta_e_ads"] = energies["interface"] - energies["slab"] - energies["pfas"]

        orientation_results[orientation_name] = {
            "atom_ranges": {
                "slab_atoms": merge_info["slab_atom_range"],
                "pfas_atoms": merge_info["pfas_atom_range"],
            },
            "energies_eh": energies,
            "merge": merge_info,
            "interface_analysis": interface_analysis,
            "xtb_runs": {
                "interface": interface_opt,
            },
            "knf": {
                "interface": {
                    "results_dir": knf_interface["results_dir"],
                    "knf_json": knf_interface["knf_json"],
                    "output_txt": knf_interface["output_txt"],
                    "metrics": knf_interface_metrics,
                },
                "delta_interface_minus_pfas": knf_delta,
            },
        }

    preferred_orientation = min(
        orientation_results,
        key=lambda name: orientation_results[name]["energies_eh"]["delta_e_ads"],
    )
    preferred_payload = orientation_results[preferred_orientation]
    anisotropy_parallel_minus_perpendicular = None
    if "parallel" in orientation_results and "perpendicular" in orientation_results:
        anisotropy_parallel_minus_perpendicular = (
            orientation_results["parallel"]["energies_eh"]["delta_e_ads"]
            - orientation_results["perpendicular"]["energies_eh"]["delta_e_ads"]
        )

    report = {
        "run_name": run_name,
        "run_directory": str(run_dir),
        "inputs": {
            "slab_xyz": str(slab_input),
            "pfas_xyz": str(pfas_input),
        },
        "settings": {
            "gfn": args.gfn,
            "placement_gap_angstrom": args.gap,
            "x_shift_angstrom": args.x_shift,
            "y_shift_angstrom": args.y_shift,
            "pfas_charge": args.pfas_charge,
            "multiplicity": args.multiplicity,
            "orientation_mode": args.orientation_mode,
            "contact_elements": contact_elements,
            "knf_scdi_var_min": args.knf_scdi_var_min,
            "knf_scdi_var_max": args.knf_scdi_var_max,
        },
        "orientation_summary": {
            "evaluated": orientation_sequence,
            "preferred_orientation": preferred_orientation,
            "delta_delta_e_ads_parallel_minus_perpendicular": anisotropy_parallel_minus_perpendicular,
        },
        "orientations": orientation_results,
        "atom_ranges": {
            "slab_atoms": preferred_payload["atom_ranges"]["slab_atoms"],
            "pfas_atoms": preferred_payload["atom_ranges"]["pfas_atoms"],
        },
        "energies_eh": preferred_payload["energies_eh"],
        "merge": preferred_payload["merge"],
        "interface_analysis": preferred_payload["interface_analysis"],
        "xtb_runs": {
            "slab": slab_opt,
            "pfas": pfas_opt,
            "interface": preferred_payload["xtb_runs"]["interface"],
        },
        "knf": {
            "pfas": {
                "results_dir": knf_pfas["results_dir"],
                "knf_json": knf_pfas["knf_json"],
                "output_txt": knf_pfas["output_txt"],
                "metrics": knf_pfas_metrics,
            },
            "interface": preferred_payload["knf"]["interface"],
            "delta_interface_minus_pfas": preferred_payload["knf"]["delta_interface_minus_pfas"],
        },
    }

    summary_json = run_dir / "summary.json"
    summary_txt = run_dir / "summary.txt"
    orientation_outputs = {}
    for orientation_name in orientation_sequence:
        payload = orientation_results[orientation_name]
        orientation_outputs[orientation_name] = {
            "interface_initial_xyz": payload["merge"]["output_xyz"],
            "interface_optimized_xyz": payload["xtb_runs"]["interface"]["optimized_xyz"],
            "knf_interface_json": payload["knf"]["interface"]["knf_json"],
            "knf_interface_txt": payload["knf"]["interface"]["output_txt"],
        }
    report["outputs"] = {
        "slab_sorted_xyz": str(slab_sorted_xyz),
        "slab_optimized_xyz": slab_opt["optimized_xyz"],
        "pfas_optimized_xyz": pfas_opt["optimized_xyz"],
        "interface_initial_xyz": preferred_payload["merge"]["output_xyz"],
        "interface_optimized_xyz": preferred_payload["xtb_runs"]["interface"]["optimized_xyz"],
        "summary_json": str(summary_json),
        "summary_txt": str(summary_txt),
        "knf_pfas_json": knf_pfas["knf_json"],
        "knf_interface_json": preferred_payload["knf"]["interface"]["knf_json"],
        "orientation_outputs": orientation_outputs,
    }

    summary_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    summary_txt.write_text(build_text_report(report), encoding="utf-8")

    print(f"Run completed: {run_dir}")
    print(f"Summary JSON: {summary_json}")
    print(f"Summary TXT: {summary_txt}")
    print(f"Delta_E_ads (Eh): {preferred_payload['energies_eh']['delta_e_ads']:.12f}")
    for orientation_name in orientation_sequence:
        delta = orientation_results[orientation_name]["energies_eh"]["delta_e_ads"]
        print(f"{orientation_name} Delta_E_ads (Eh): {delta:.12f}")
    print(f"Preferred orientation: {preferred_orientation}")
    if anisotropy_parallel_minus_perpendicular is not None:
        print(
            "Delta_Delta_E_ads (parallel - perpendicular) (Eh): "
            f"{anisotropy_parallel_minus_perpendicular:.12f}"
        )
    closest = preferred_payload["interface_analysis"]["closest_contact"]
    if closest is not None:
        print(f"Min PFAS contact distance (A): {closest['distance_angstrom']:.4f}")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    run_workflow(args)


if __name__ == "__main__":
    main(sys.argv[1:])
