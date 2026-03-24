import argparse
from pathlib import Path


def read_xyz(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path} is not a valid XYZ file.")

    n_atoms = int(lines[0].strip())
    comment = lines[1].rstrip("\n")
    atom_lines = lines[2 : 2 + n_atoms]
    if len(atom_lines) != n_atoms:
        raise ValueError(f"{path} declares {n_atoms} atoms but does not contain that many lines.")

    atoms = []
    for line in atom_lines:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed XYZ line in {path}: {line}")
        el = parts[0]
        x, y, z = map(float, parts[1:4])
        atoms.append((el, x, y, z))
    return comment, atoms


def bounds(atoms):
    xs = [atom[1] for atom in atoms]
    ys = [atom[2] for atom in atoms]
    zs = [atom[3] for atom in atoms]
    return (
        min(xs),
        max(xs),
        min(ys),
        max(ys),
        min(zs),
        max(zs),
    )


def center_xy(atoms):
    min_x, max_x, min_y, max_y, _, _ = bounds(atoms)
    return (0.5 * (min_x + max_x), 0.5 * (min_y + max_y))


def translate(atoms, dx, dy, dz):
    return [(el, x + dx, y + dy, z + dz) for el, x, y, z in atoms]


def write_xyz(path, comment, atoms):
    with open(path, "w") as handle:
        handle.write(f"{len(atoms)}\n")
        handle.write(f"{comment}\n")
        for el, x, y, z in atoms:
            handle.write(f"{el:2s} {x: .6f} {y: .6f} {z: .6f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Place a PFAS XYZ above a slab XYZ with a fixed vertical gap."
    )
    parser.add_argument("slab_xyz", help="Relaxed slab XYZ, e.g. xtbopt.xyz")
    parser.add_argument("pfas_xyz", help="PFAS XYZ to place above the slab")
    parser.add_argument("output_xyz", help="Merged slab+PFAS XYZ output")
    parser.add_argument(
        "--gap",
        type=float,
        default=2.5,
        help="Vertical gap between slab zmax and PFAS zmin in angstrom (default: 2.5)",
    )
    parser.add_argument(
        "--x-shift",
        type=float,
        default=0.0,
        help="Additional x shift applied after centering (default: 0.0)",
    )
    parser.add_argument(
        "--y-shift",
        type=float,
        default=0.0,
        help="Additional y shift applied after centering (default: 0.0)",
    )
    args = parser.parse_args()

    slab_comment, slab_atoms = read_xyz(args.slab_xyz)
    pfas_comment, pfas_atoms = read_xyz(args.pfas_xyz)

    slab_min_x, slab_max_x, slab_min_y, slab_max_y, slab_min_z, slab_max_z = bounds(slab_atoms)
    pfas_min_x, pfas_max_x, pfas_min_y, pfas_max_y, pfas_min_z, pfas_max_z = bounds(pfas_atoms)

    slab_cx, slab_cy = center_xy(slab_atoms)
    pfas_cx, pfas_cy = center_xy(pfas_atoms)

    dx = slab_cx - pfas_cx + args.x_shift
    dy = slab_cy - pfas_cy + args.y_shift
    dz = (slab_max_z + args.gap) - pfas_min_z

    placed_pfas = translate(pfas_atoms, dx, dy, dz)
    merged_atoms = slab_atoms + placed_pfas

    comment = (
        f"Slab: {Path(args.slab_xyz).name}; PFAS: {Path(args.pfas_xyz).name}; "
        f"gap={args.gap:.3f} A; slab_atoms=1-{len(slab_atoms)}; "
        f"pfas_atoms={len(slab_atoms)+1}-{len(merged_atoms)}"
    )
    write_xyz(args.output_xyz, comment, merged_atoms)

    _, _, _, _, placed_pfas_min_z, placed_pfas_max_z = bounds(placed_pfas)

    print(f"Output written to: {args.output_xyz}")
    print(f"Slab atoms: 1-{len(slab_atoms)}")
    print(f"PFAS atoms: {len(slab_atoms)+1}-{len(merged_atoms)}")
    print(f"Applied translation: dx={dx:.6f} dy={dy:.6f} dz={dz:.6f}")
    print(
        f"Slab bbox: x=[{slab_min_x:.3f}, {slab_max_x:.3f}] "
        f"y=[{slab_min_y:.3f}, {slab_max_y:.3f}] z=[{slab_min_z:.3f}, {slab_max_z:.3f}]"
    )
    print(
        f"Placed PFAS bbox: x=[{pfas_min_x + dx:.3f}, {pfas_max_x + dx:.3f}] "
        f"y=[{pfas_min_y + dy:.3f}, {pfas_max_y + dy:.3f}] "
        f"z=[{placed_pfas_min_z:.3f}, {placed_pfas_max_z:.3f}]"
    )


if __name__ == "__main__":
    main()
