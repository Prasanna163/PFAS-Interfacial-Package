# PFAS Interfacial Package

This project packages the PFAS water-slab adsorption workflow as a CLI.

It performs:

- slab sorting by water-layer z order
- frozen-bottom slab optimization with xTB
- PFAS-only optimization
- PFAS placement above the optimized water slab
- frozen-slab interface optimization
- adsorption-energy analysis
- z-range and closest-contact analysis
- KNF analysis on the final optimized PFAS and interface geometries

The repository includes the reusable optimized slab master:

- `master_slab_xtbopt.xyz`
- `master_slab_xtbopt.mol`
- `master_slab_xtbtopo.mol`

## Install

From this folder:

```powershell
python -m pip install -e .
```

## CLI

Primary command:

```powershell
PFAS <input file or directory>
```

Example:

```powershell
PFAS ".\molecules\PFOA.xyz"
```

Directory mode:

```powershell
PFAS ".\molecules"
```

By default the CLI uses `master_slab_xtbopt.xyz` from the current working directory, or falls back to the package copy bundled with this project.

If you want to override the slab:

```powershell
PFAS ".\molecules\PFOA.xyz" --slab ".\another_slab.xyz"
```

Example with charge:

```powershell
PFAS ".\molecules\PFOS.xyz" --pfas-charge -1
```

The older explicit command is also still available:

```powershell
pfas-interface-knf slab.xyz ".\molecules\PFOA.xyz"
```

## Notes

- `xtb` must be available on `PATH`.
- The local KNF installation must be working.
- KNF runs its single-point/NCI stage on the final optimized geometry; the workflow avoids KNF geometry reoptimization.
