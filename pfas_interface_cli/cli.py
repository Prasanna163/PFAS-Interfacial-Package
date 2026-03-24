import argparse
from pathlib import Path

from .workflow import main as workflow_main

ALLOWED_EXTENSIONS = {".xyz", ".mol", ".sdf", ".mol2"}
MASTER_SLAB_FILENAME = "master_slab_xtbopt.xyz"


def _default_slab_candidates():
    project_root = Path(__file__).resolve().parent.parent
    return [
        Path.cwd() / MASTER_SLAB_FILENAME,
        project_root / "references" / MASTER_SLAB_FILENAME,
        project_root / MASTER_SLAB_FILENAME,  # legacy location fallback
    ]


def _default_slab_path():
    for candidate in _default_slab_candidates():
        if candidate.exists():
            return candidate
    return _default_slab_candidates()[0]


def _collect_inputs(input_path, slab_path):
    input_path = Path(input_path).resolve()
    slab_path = Path(slab_path).resolve()

    if input_path.is_file():
        return [input_path]

    if input_path.is_dir():
        files = []
        for candidate in sorted(input_path.iterdir()):
            if not candidate.is_file():
                continue
            if candidate.suffix.lower() not in ALLOWED_EXTENSIONS:
                continue
            if candidate.resolve() == slab_path:
                continue
            if candidate.name.lower().startswith("master_slab"):
                continue
            files.append(candidate.resolve())
        if not files:
            raise FileNotFoundError(f"No supported PFAS input files found in {input_path}.")
        return files

    raise FileNotFoundError(f"Input path not found: {input_path}")


def build_parser():
    parser = argparse.ArgumentParser(
        description="PFAS adsorption workflow CLI using the default master slab."
    )
    parser.add_argument("input_path", help="PFAS input file or directory")
    parser.add_argument(
        "--slab",
        default=str(_default_slab_path()),
        help="Slab XYZ to use. Default: master_slab_xtbopt.xyz in cwd or references folder.",
    )
    parser.add_argument("--pfas-charge", type=int, default=0)
    advanced = parser.add_argument_group("Advanced Options (Second Layer)")
    advanced.add_argument("--run-root", default="pipeline_runs", help="Output run root directory")
    advanced.add_argument("--run-name", default=None, help="Optional run name for single-file mode")
    advanced.add_argument("--gap", type=float, default=2.5)
    advanced.add_argument("--x-shift", type=float, default=0.0)
    advanced.add_argument("--y-shift", type=float, default=0.0)
    advanced.add_argument(
        "--orientation-mode",
        choices=["dual", "perpendicular", "parallel"],
        default="dual",
    )
    advanced.add_argument("--multiplicity", type=int, default=1)
    advanced.add_argument("--gfn", type=int, default=2)
    advanced.add_argument("--contact-elements", default="O,F")
    advanced.add_argument("--keep-knf-intermediates", action="store_true")
    advanced.add_argument("--knf-scdi-var-min", type=float, default=None)
    advanced.add_argument("--knf-scdi-var-max", type=float, default=None)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    slab_path = Path(args.slab).resolve()
    if not slab_path.exists():
        raise FileNotFoundError(
            f"Slab file not found: {slab_path}. "
            "Use --slab or place master_slab_xtbopt.xyz in the working folder or references folder."
        )

    input_files = _collect_inputs(args.input_path, slab_path)
    run_root = Path(args.run_root)
    if len(input_files) > 1 and args.run_name:
        raise ValueError("--run-name is only valid for single-file mode.")

    for input_file in input_files:
        workflow_argv = [
            str(slab_path),
            str(input_file),
            "--run-root",
            str(run_root),
            "--gap",
            str(args.gap),
            "--x-shift",
            str(args.x_shift),
            "--y-shift",
            str(args.y_shift),
            "--orientation-mode",
            args.orientation_mode,
            "--pfas-charge",
            str(args.pfas_charge),
            "--multiplicity",
            str(args.multiplicity),
            "--gfn",
            str(args.gfn),
            "--contact-elements",
            args.contact_elements,
        ]
        if args.keep_knf_intermediates:
            workflow_argv.append("--keep-knf-intermediates")
        if args.knf_scdi_var_min is not None:
            workflow_argv.extend(["--knf-scdi-var-min", str(args.knf_scdi_var_min)])
        if args.knf_scdi_var_max is not None:
            workflow_argv.extend(["--knf-scdi-var-max", str(args.knf_scdi_var_max)])
        if len(input_files) == 1 and args.run_name:
            workflow_argv.extend(["--run-name", args.run_name])
        workflow_main(workflow_argv)


__all__ = ["main", "build_parser"]
