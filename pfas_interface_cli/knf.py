import json
import shutil
from pathlib import Path

from knf_core.pipeline import KNFPipeline


def run_knf_on_existing_optimized_xyz(
    input_xyz,
    output_root,
    charge,
    multiplicity,
    keep_full_files=False,
    scdi_var_min=None,
    scdi_var_max=None,
):
    pipeline = KNFPipeline(
        input_file=str(input_xyz),
        charge=charge,
        spin=multiplicity,
        water=False,
        force=False,
        clean=False,
        debug=False,
        output_root=str(output_root),
        keep_full_files=keep_full_files,
        nci_backend="torch",
        scdi_var_min=scdi_var_min,
        scdi_var_max=scdi_var_max,
    )
    pipeline.setup_directories()
    seeded_xtbopt = Path(pipeline.results_dir) / "xtbopt.xyz"
    shutil.copy2(input_xyz, seeded_xtbopt)
    context = pipeline.run_pre_nci_stage()
    pipeline.run_post_nci_stage(context)

    knf_json = Path(pipeline.results_dir) / "knf.json"
    output_txt = Path(pipeline.results_dir) / "output.txt"
    if not knf_json.exists():
        raise FileNotFoundError(f"KNF output not found: {knf_json}")

    return {
        "results_dir": str(pipeline.results_dir),
        "knf_json": str(knf_json),
        "output_txt": str(output_txt),
        "data": json.loads(knf_json.read_text(encoding="utf-8")),
    }


__all__ = ["run_knf_on_existing_optimized_xyz"]
