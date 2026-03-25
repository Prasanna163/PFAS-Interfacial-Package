import math
import shutil
from datetime import datetime
from pathlib import Path

from .geometry import sort_water_slab_by_oxygen_z, write_xyz
from .xtb import optimize_with_xtb


DEFAULT_CUSTOM_SLAB_WORKDIR = "slab_library/custom_builds"


def _load_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:  # pragma: no cover - depends on local RDKit install
        raise RuntimeError(
            "RDKit is required for custom slab generation. "
            "Install RDKit, then retry custom slab creation."
        ) from exc
    return Chem, AllChem


def _rotate_z(point, angle_rad):
    x, y, z = point
    cosine = math.cos(angle_rad)
    sine = math.sin(angle_rad)
    return (x * cosine - y * sine, x * sine + y * cosine, z)


def _water_template_atoms():
    Chem, AllChem = _load_rdkit()
    molecule = Chem.AddHs(Chem.MolFromSmiles("O"))
    params = AllChem.ETKDGv3()
    params.randomSeed = 104729
    embed_code = AllChem.EmbedMolecule(molecule, params)
    if embed_code != 0:
        raise RuntimeError("RDKit failed to embed a 3D water conformer.")
    AllChem.UFFOptimizeMolecule(molecule, maxIters=500)

    conformer = molecule.GetConformer()
    oxygen_index = next(atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetSymbol() == "O")
    oxygen_position = conformer.GetAtomPosition(oxygen_index)
    centered = []
    for atom in molecule.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        centered.append(
            {
                "element": atom.GetSymbol(),
                "x": float(position.x - oxygen_position.x),
                "y": float(position.y - oxygen_position.y),
                "z": float(position.z - oxygen_position.z),
            }
        )

    oxygen = [atom for atom in centered if atom["element"] == "O"]
    hydrogens = [atom for atom in centered if atom["element"] == "H"]
    if len(oxygen) != 1 or len(hydrogens) != 2:
        raise RuntimeError("Unexpected RDKit water template topology.")
    return oxygen + hydrogens


def _grid_count(length, spacing):
    return max(1, int(math.ceil(length / spacing)))


def generate_water_slab_atoms(length_x, length_y, length_z, spacing_xy=3.1, spacing_z=2.9):
    if min(length_x, length_y, length_z) <= 0.0:
        raise ValueError("Custom slab dimensions must be positive.")
    if spacing_xy <= 0.0 or spacing_z <= 0.0:
        raise ValueError("Custom slab spacing values must be positive.")

    nx = _grid_count(length_x, spacing_xy)
    ny = _grid_count(length_y, spacing_xy)
    nz = _grid_count(length_z, spacing_z)
    water_template = _water_template_atoms()
    atoms = []
    for iz in range(nz):
        for ix in range(nx):
            for iy in range(ny):
                x_shift = (ix - 0.5 * (nx - 1)) * spacing_xy
                y_shift = (iy - 0.5 * (ny - 1)) * spacing_xy
                z_shift = iz * spacing_z
                angle = 0.5 * math.pi * ((ix + iy + iz) % 4)
                for atom in water_template:
                    rx, ry, rz = _rotate_z((atom["x"], atom["y"], atom["z"]), angle)
                    atoms.append(
                        {
                            "index": len(atoms) + 1,
                            "element": atom["element"],
                            "x": rx + x_shift,
                            "y": ry + y_shift,
                            "z": rz + z_shift,
                        }
                    )
    return atoms, {"nx": nx, "ny": ny, "nz": nz, "water_count": nx * ny * nz}


def write_custom_slab_xyz(
    output_xyz,
    length_x,
    length_y,
    length_z,
    spacing_xy=3.1,
    spacing_z=2.9,
):
    output_xyz = Path(output_xyz)
    output_xyz.parent.mkdir(parents=True, exist_ok=True)
    atoms, grid = generate_water_slab_atoms(
        length_x=length_x,
        length_y=length_y,
        length_z=length_z,
        spacing_xy=spacing_xy,
        spacing_z=spacing_z,
    )
    comment = (
        f"RDKit-generated water slab; dims=({length_x:.3f},{length_y:.3f},{length_z:.3f}) A; "
        f"grid=({grid['nx']},{grid['ny']},{grid['nz']}); spacing=({spacing_xy:.3f},{spacing_z:.3f}) A"
    )
    write_xyz(output_xyz, comment, atoms)
    return {
        "raw_xyz": str(output_xyz),
        "grid": grid,
        "dimensions_angstrom": {"x": length_x, "y": length_y, "z": length_z},
        "spacing_angstrom": {"xy": spacing_xy, "z": spacing_z},
    }


def build_optimize_custom_slab(
    work_root,
    length_x,
    length_y,
    length_z,
    spacing_xy=3.1,
    spacing_z=2.9,
    gfn=2,
    default_slab_path=None,
):
    work_root = Path(work_root).resolve()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = work_root / f"custom_slab_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=False)

    raw_xyz = run_dir / "custom_slab_raw.xyz"
    sorted_xyz = run_dir / "custom_slab_sorted.xyz"
    slab_info = write_custom_slab_xyz(
        raw_xyz,
        length_x=length_x,
        length_y=length_y,
        length_z=length_z,
        spacing_xy=spacing_xy,
        spacing_z=spacing_z,
    )
    sort_water_slab_by_oxygen_z(raw_xyz, sorted_xyz)
    xtb_result = optimize_with_xtb(sorted_xyz, run_dir / "xtb_opt", charge=0, multiplicity=1, gfn=gfn)

    default_target = None
    if default_slab_path is not None:
        default_target = Path(default_slab_path).resolve()
        default_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(Path(xtb_result["optimized_xyz"]).resolve(), default_target)

    return {
        "run_dir": str(run_dir),
        "raw_xyz": slab_info["raw_xyz"],
        "sorted_xyz": str(sorted_xyz),
        "optimized_xyz": xtb_result["optimized_xyz"],
        "xtb": xtb_result,
        "grid": slab_info["grid"],
        "dimensions_angstrom": slab_info["dimensions_angstrom"],
        "spacing_angstrom": slab_info["spacing_angstrom"],
        "default_slab_path": str(default_target) if default_target else None,
    }


__all__ = [
    "DEFAULT_CUSTOM_SLAB_WORKDIR",
    "build_optimize_custom_slab",
    "generate_water_slab_atoms",
    "write_custom_slab_xyz",
]
