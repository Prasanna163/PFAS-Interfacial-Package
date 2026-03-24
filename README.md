# PFAS Interfacial Package

This project packages the PFAS water-slab adsorption workflow as a CLI.

It performs:

- slab sorting by water-layer z order
- slab optimization with xTB
- PFAS-only optimization
- dual-orientation PFAS placement above the optimized water slab (`perpendicular` and `parallel`)
- interface optimization for each orientation branch
- adsorption-energy analysis
- z-range and closest-contact analysis
- KNF analysis on the final optimized PFAS and interface geometries

The repository includes the reusable optimized slab master:

- `references/master_slab_xtbopt.xyz`
- `references/master_slab_xtbopt.mol`
- `references/master_slab_xtbtopo.mol`

Project layout highlights:

- `pfas_interface_cli/` -> installed CLI package and workflow modules:
  `workflow.py`, `geometry.py`, `xtb.py`, `placement.py`, `analysis.py`, `reporting.py`, `knf.py`
- `run_pfas_interface_knf.py` -> compatibility shim that calls `pfas_interface_cli.workflow`
- `references/` -> reusable slab templates
- `tools/` -> helper and analysis scripts
- `examples/` -> sample input files (`slab.xyz`, `slab.mol`)
- `docs/` -> manuscript and documentation assets

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

By default the CLI uses `master_slab_xtbopt.xyz` from the current working directory, or falls back to `references/master_slab_xtbopt.xyz` in this project.

If you want to override the slab:

```powershell
PFAS ".\molecules\PFOA.xyz" --slab ".\another_slab.xyz"
```

Example with charge:

```powershell
PFAS ".\molecules\PFOS.xyz" --pfas-charge -1
```

Advanced options are available as a second layer when needed (for reproducibility sweeps/tuning):

- `--run-root`, `--run-name`
- `--gap`, `--x-shift`, `--y-shift`
- `--orientation-mode` (`dual|perpendicular|parallel`)
- `--multiplicity`, `--gfn`
- `--contact-elements`
- `--keep-knf-intermediates`
- `--knf-scdi-var-min`, `--knf-scdi-var-max`

Example (advanced):

```powershell
PFAS ".\molecules\PFOS.xyz" --pfas-charge -1 --orientation-mode perpendicular --gfn 2 --run-name PFOS_test
```

The older explicit command is also still available:

```powershell
pfas-interface-knf .\examples\slab.xyz ".\molecules\PFOA.xyz"
```

## Notes

- `xtb` must be available on `PATH`.
- The local KNF installation must be working.
- KNF runs its single-point/NCI stage on the final optimized geometry; the workflow avoids KNF geometry reoptimization.
