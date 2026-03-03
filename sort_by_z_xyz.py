import sys

inp = sys.argv[1]
out = sys.argv[2]

lines = open(inp).read().splitlines()
n = int(lines[0].strip())
comment = lines[1].strip()

if n % 3 != 0:
    raise ValueError("Expected a water slab with atom count divisible by 3.")

waters = []
for i in range(0, n, 3):
    atoms = []
    for j in range(3):
        el, x, y, z = lines[2 + i + j].split()[:4]
        x, y, z = float(x), float(y), float(z)
        atoms.append((el, x, y, z))

    if [atom[0] for atom in atoms] != ["O", "H", "H"]:
        raise ValueError(f"Expected O-H-H ordering for water starting at atom {i + 1}.")

    waters.append((atoms[0][3], atoms))

waters.sort(key=lambda t: t[0])

with open(out, "w") as f:
    f.write(f"{n}\n{comment}\n")
    for _, atoms in waters:
        for el, x, y, z in atoms:
            f.write(f"{el:2s} {x: .6f} {y: .6f} {z: .6f}\n")
