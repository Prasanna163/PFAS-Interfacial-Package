import argparse
import json
from datetime import datetime
from pathlib import Path

from .analysis import KNF_METRIC_NAMES, analyze_interface_geometry, extract_knf_metrics, metric_delta
from .geometry import sort_water_slab_by_oxygen_z
from .knf import run_knf_on_existing_optimized_xyz
from .placement import merge_pfas_above_slab
from .reporting import build_text_report
from .xtb import optimize_with_xtb


# Production defaults.
DEFAULT_RUN_ROOT = "pipeline_runs"
DEFAULT_ORIENTATION_SEQUENCE = ["perpendicular", "parallel"]
DEFAULT_GFN_LEVEL = 2
DEFAULT_PFAS_MULTIPLICITY = 1
DEFAULT_GAP_ANGSTROM = 1.5
DEFAULT_X_SHIFT_ANGSTROM = 0.0
DEFAULT_Y_SHIFT_ANGSTROM = 0.0
DEFAULT_CONTACT_ELEMENTS = ["O", "F"]
DEFAULT_KEEP_KNF_INTERMEDIATES = False
DEFAULT_KNF_SCDI_VAR_MIN = None
DEFAULT_KNF_SCDI_VAR_MAX = None


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


__all__ = [
    "DEFAULT_CONTACT_ELEMENTS",
    "DEFAULT_GAP_ANGSTROM",
    "DEFAULT_GFN_LEVEL",
    "DEFAULT_KEEP_KNF_INTERMEDIATES",
    "DEFAULT_KNF_SCDI_VAR_MAX",
    "DEFAULT_KNF_SCDI_VAR_MIN",
    "DEFAULT_ORIENTATION_SEQUENCE",
    "DEFAULT_PFAS_MULTIPLICITY",
    "DEFAULT_RUN_ROOT",
    "DEFAULT_X_SHIFT_ANGSTROM",
    "DEFAULT_Y_SHIFT_ANGSTROM",
    "build_parser",
    "main",
    "run_workflow",
]
