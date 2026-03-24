import re
import shutil
import subprocess

from .geometry import read_xyz


ENERGY_RE = re.compile(r"energy:\s*([-+]?\d+(?:\.\d+)?)", re.IGNORECASE)


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


__all__ = [
    "ENERGY_RE",
    "optimize_with_xtb",
    "parse_xtb_energy_from_xyz",
    "run_command",
]
