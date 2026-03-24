import math
from pathlib import Path


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
        trial = (1.0, 0.0, 0.0) if abs(source[0]) < 0.9 else (0.0, 1.0, 0.0)
        axis = normalize(cross(source, trial))
        x, y, z = axis
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


__all__ = [
    "centroid",
    "cross",
    "dot",
    "identity_matrix",
    "mat_vec_mul",
    "norm",
    "normalize",
    "read_xyz",
    "rotate_atoms_about_point",
    "rotation_matrix_from_vectors",
    "sort_water_slab_by_oxygen_z",
    "write_xyz",
    "xy_center",
    "xyz_bounds",
]
