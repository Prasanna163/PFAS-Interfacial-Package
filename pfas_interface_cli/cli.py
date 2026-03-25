import argparse
import re
from pathlib import Path

from .interactive import choose_option, ensure_interactive_terminal, prompt_float, prompt_int, prompt_text
from .slab_builder import DEFAULT_CUSTOM_SLAB_WORKDIR, build_optimize_custom_slab
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


def _slug(text):
    slug = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
    return slug or "slab"


def _discover_preoptimized_slabs():
    project_root = Path(__file__).resolve().parent.parent
    cwd_master = (Path.cwd() / MASTER_SLAB_FILENAME).resolve()
    reference_master = (project_root / "references" / MASTER_SLAB_FILENAME).resolve()
    default_candidate_set = {candidate.resolve() for candidate in _default_slab_candidates() if candidate.exists()}

    search_roots = [
        Path.cwd(),
        project_root / "references",
        Path.cwd() / "slab_library",
        project_root / "slab_library",
    ]

    discovered_paths = []
    for candidate in _default_slab_candidates():
        if candidate.exists():
            discovered_paths.append(candidate.resolve())

    for root in search_roots:
        if not root.exists():
            continue
        for candidate in sorted(root.rglob("*.xyz")):
            resolved = candidate.resolve()
            name = resolved.name.lower()
            if "slab" not in name:
                continue
            if resolved not in default_candidate_set and "xtbopt" not in name and "optimized" not in name:
                continue
            discovered_paths.append(resolved)

    unique_paths = []
    seen = set()
    for path in discovered_paths:
        if path in seen:
            continue
        seen.add(path)
        unique_paths.append(path)

    named = {}
    for path in unique_paths:
        if path == cwd_master:
            base_name = "cwd-master"
        elif path == reference_master:
            base_name = "reference-master"
        else:
            base_name = _slug(path.stem)
        name = base_name
        counter = 2
        while name in named:
            name = f"{base_name}-{counter}"
            counter += 1
        named[name] = path
    return dict(sorted(named.items(), key=lambda item: item[0]))


def _resolve_preoptimized_slab(selector, discovered_slabs):
    candidate_path = Path(selector).expanduser()
    if candidate_path.exists():
        return candidate_path.resolve()
    normalized = selector.strip().lower()
    for name, slab_path in discovered_slabs.items():
        if name.lower() == normalized:
            return slab_path.resolve()
    valid_names = ", ".join(discovered_slabs.keys()) if discovered_slabs else "(none found)"
    raise FileNotFoundError(
        f"Could not resolve preoptimized slab '{selector}'. "
        f"Use --list-preoptimized-slabs to inspect available names. Available: {valid_names}"
    )


def _configure_custom_slab_from_prompts(args):
    args.create_custom_slab = True
    args.preoptimized_slab = None
    args.custom_slab_x = prompt_float(
        "Custom slab X dimension (angstrom)",
        default=args.custom_slab_x,
        min_value=0.1,
    )
    args.custom_slab_y = prompt_float(
        "Custom slab Y dimension (angstrom)",
        default=args.custom_slab_y,
        min_value=0.1,
    )
    args.custom_slab_z = prompt_float(
        "Custom slab Z thickness (angstrom)",
        default=args.custom_slab_z,
        min_value=0.1,
    )
    args.custom_slab_spacing_xy = prompt_float(
        "Water spacing in XY (angstrom)",
        default=args.custom_slab_spacing_xy,
        min_value=0.1,
    )
    args.custom_slab_spacing_z = prompt_float(
        "Layer spacing in Z (angstrom)",
        default=args.custom_slab_spacing_z,
        min_value=0.1,
    )
    args.custom_slab_workdir = prompt_text(
        "Custom slab working directory",
        default=args.custom_slab_workdir,
    )
    set_default_choice = choose_option(
        "Set optimized custom slab as future default?",
        ["Yes", "No"],
    )
    args.no_set_custom_slab_default = set_default_choice == 1


def _interactive_wizard(args):
    ensure_interactive_terminal()
    discovered_slabs = _discover_preoptimized_slabs()
    default_slab = _default_slab_path().resolve()
    top_choice = choose_option(
        "PFAS Interactive Menu",
        [
            "Run PFAS workflow",
            "Create custom slab and exit",
            "Exit",
        ],
        footer=f"Current default slab: {default_slab}",
    )
    if top_choice == 2:
        return None

    if top_choice == 1:
        _configure_custom_slab_from_prompts(args)
        args.custom_slab_only = True
        args.input_path = None
        return args

    args.input_path = prompt_text("PFAS input file or directory", default=args.input_path)
    args.pfas_charge = prompt_int("PFAS charge", default=args.pfas_charge)
    orientation_modes = ["dual", "perpendicular", "parallel"]
    mode_index = choose_option("Select orientation mode", orientation_modes)
    args.orientation_mode = orientation_modes[mode_index]

    slab_choice = choose_option(
        "Choose slab source",
        [
            "Use current default slab",
            "Select preoptimized slab",
            "Build custom slab now",
        ],
    )
    if slab_choice == 0:
        args.preoptimized_slab = None
        args.create_custom_slab = False
        args.custom_slab_only = False
    elif slab_choice == 1:
        if not discovered_slabs:
            raise FileNotFoundError("No preoptimized slabs were found for selection.")
        names = list(discovered_slabs.keys())
        labels = [f"{name} -> {discovered_slabs[name]}" for name in names]
        selected_index = choose_option("Select preoptimized slab", labels)
        args.preoptimized_slab = names[selected_index]
        args.create_custom_slab = False
        args.custom_slab_only = False
    else:
        _configure_custom_slab_from_prompts(args)
        args.custom_slab_only = False

    args.interactive = False
    return args


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
        description=(
            "PFAS adsorption workflow CLI with slab management support "
            "(default slab, preoptimized slabs, and custom RDKit/xTB slab builder)."
        )
    )
    parser.add_argument("input_path", nargs="?", help="PFAS input file or directory")
    parser.add_argument(
        "--slab",
        default=str(_default_slab_path()),
        help="Slab XYZ to use. Default: master_slab_xtbopt.xyz in cwd or references folder.",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Launch an arrow-key interactive UI for slab selection/building and workflow setup.",
    )
    parser.add_argument(
        "--list-preoptimized-slabs",
        action="store_true",
        help="List discovered preoptimized slab names for --preoptimized-slab and exit.",
    )
    parser.add_argument(
        "--preoptimized-slab",
        default=None,
        help=(
            "Select a discovered preoptimized slab by name (see --list-preoptimized-slabs), "
            "or pass a direct XYZ path."
        ),
    )
    parser.add_argument(
        "--create-custom-slab",
        action="store_true",
        help="Build a custom water slab via RDKit, optimize with xTB, and use it for this run.",
    )
    parser.add_argument(
        "--custom-slab-only",
        action="store_true",
        help="Build and optimize the custom slab, set default if selected, then exit.",
    )
    parser.add_argument(
        "--custom-slab-x",
        type=float,
        default=20.0,
        help="Custom slab target X dimension in angstrom (default: 20.0).",
    )
    parser.add_argument(
        "--custom-slab-y",
        type=float,
        default=20.0,
        help="Custom slab target Y dimension in angstrom (default: 20.0).",
    )
    parser.add_argument(
        "--custom-slab-z",
        type=float,
        default=10.0,
        help="Custom slab target Z thickness in angstrom (default: 10.0).",
    )
    parser.add_argument(
        "--custom-slab-spacing-xy",
        type=float,
        default=3.1,
        help="Custom slab water spacing in XY (angstrom, default: 3.1).",
    )
    parser.add_argument(
        "--custom-slab-spacing-z",
        type=float,
        default=2.9,
        help="Custom slab water-layer spacing in Z (angstrom, default: 2.9).",
    )
    parser.add_argument(
        "--custom-slab-workdir",
        default=DEFAULT_CUSTOM_SLAB_WORKDIR,
        help=(
            "Folder where custom slab generation/optimization runs are stored "
            f"(default: {DEFAULT_CUSTOM_SLAB_WORKDIR})."
        ),
    )
    parser.add_argument(
        "--no-set-custom-slab-default",
        action="store_true",
        help="Do not copy optimized custom slab to ./master_slab_xtbopt.xyz.",
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

    if args.interactive:
        args = _interactive_wizard(args)
        if args is None:
            print("Interactive session closed without running the workflow.")
            return

    discovered_slabs = _discover_preoptimized_slabs()
    if args.list_preoptimized_slabs:
        if not discovered_slabs:
            print("No preoptimized slabs found.")
            return
        print("Available preoptimized slabs:")
        for name, slab in discovered_slabs.items():
            print(f"  {name:20s} {slab}")
        return

    if args.preoptimized_slab and args.create_custom_slab:
        raise ValueError("Use either --preoptimized-slab or --create-custom-slab, not both.")
    if args.custom_slab_only and not args.create_custom_slab:
        args.create_custom_slab = True

    slab_path = None
    if args.preoptimized_slab:
        slab_path = _resolve_preoptimized_slab(args.preoptimized_slab, discovered_slabs)
    elif args.create_custom_slab:
        default_target = None
        if not args.no_set_custom_slab_default:
            default_target = Path.cwd() / MASTER_SLAB_FILENAME
        try:
            custom_slab_result = build_optimize_custom_slab(
                work_root=args.custom_slab_workdir,
                length_x=args.custom_slab_x,
                length_y=args.custom_slab_y,
                length_z=args.custom_slab_z,
                spacing_xy=args.custom_slab_spacing_xy,
                spacing_z=args.custom_slab_spacing_z,
                gfn=args.gfn,
                default_slab_path=default_target,
            )
        except Exception as exc:
            raise RuntimeError(
                "Custom slab generation failed. Ensure RDKit is installed and xTB runs correctly."
            ) from exc
        slab_path = Path(
            custom_slab_result["default_slab_path"] or custom_slab_result["optimized_xyz"]
        ).resolve()
        print(f"Custom slab generated: {custom_slab_result['optimized_xyz']}")
        print(f"Custom slab grid (nx,ny,nz): {custom_slab_result['grid']}")
        if custom_slab_result["default_slab_path"]:
            print(f"Updated default slab: {custom_slab_result['default_slab_path']}")
        if args.custom_slab_only:
            return
    else:
        slab_path = Path(args.slab).resolve()

    if not slab_path.exists():
        raise FileNotFoundError(
            f"Slab file not found: {slab_path}. "
            "Use --slab or place master_slab_xtbopt.xyz in the working folder or references folder."
        )

    if not args.input_path:
        raise ValueError("input_path is required unless --custom-slab-only is used.")

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
