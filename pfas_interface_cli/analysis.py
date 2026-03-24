import math

from .geometry import read_xyz, xyz_bounds


KNF_METRIC_NAMES = ["SNCI", "SCDI", "SCDI_variance"] + [f"f{i}" for i in range(1, 10)]


def distance(atom_a, atom_b):
    return math.sqrt(
        (atom_a["x"] - atom_b["x"]) ** 2
        + (atom_a["y"] - atom_b["y"]) ** 2
        + (atom_a["z"] - atom_b["z"]) ** 2
    )


def analyze_interface_geometry(interface_xyz, slab_atom_count, contact_elements):
    _, atoms = read_xyz(interface_xyz)
    slab_atoms = atoms[:slab_atom_count]
    pfas_atoms = atoms[slab_atom_count:]

    slab_bounds = xyz_bounds(slab_atoms)
    pfas_bounds = xyz_bounds(pfas_atoms)
    total_bounds = xyz_bounds(atoms)

    water_h_atoms = [atom for atom in slab_atoms if atom["element"] == "H"]
    pfas_contact_atoms = [atom for atom in pfas_atoms if atom["element"] in contact_elements]
    if not pfas_contact_atoms:
        pfas_contact_atoms = [atom for atom in pfas_atoms if atom["element"] not in {"C", "H"}]

    closest_pair = None
    for pfas_atom in pfas_contact_atoms:
        for water_atom in water_h_atoms:
            current_distance = distance(pfas_atom, water_atom)
            if closest_pair is None or current_distance < closest_pair["distance_angstrom"]:
                closest_pair = {
                    "distance_angstrom": current_distance,
                    "pfas_atom": pfas_atom,
                    "water_atom": water_atom,
                }

    min_gap = pfas_bounds["z_min"] - slab_bounds["z_max"]
    max_gap = pfas_bounds["z_max"] - slab_bounds["z_max"]

    if closest_pair is None:
        contact_classification = "No PFAS contact atoms found for the requested element filter."
    else:
        distance_value = closest_pair["distance_angstrom"]
        if distance_value <= 2.2:
            contact_classification = "Hydrogen-bond-like / strong specific contact"
        elif distance_value <= 3.5:
            contact_classification = "Weak interfacial contact"
        else:
            contact_classification = "Hovering / no direct specific contact"

    if pfas_bounds["z_min"] >= slab_bounds["z_max"]:
        interface_classification = "Surface-localized above slab top"
    elif pfas_bounds["z_max"] <= slab_bounds["z_max"]:
        interface_classification = "Immersed within slab region"
    else:
        interface_classification = "Partially penetrating interface"

    return {
        "total_bounds": total_bounds,
        "slab_bounds": slab_bounds,
        "pfas_bounds": pfas_bounds,
        "slab_thickness_angstrom": slab_bounds["z_max"] - slab_bounds["z_min"],
        "pfas_offset_from_slab_top": {
            "z_min_minus_slab_top": min_gap,
            "z_max_minus_slab_top": max_gap,
        },
        "interface_classification": interface_classification,
        "contact_elements_used": sorted({atom["element"] for atom in pfas_contact_atoms}),
        "closest_contact": closest_pair,
        "contact_classification": contact_classification,
    }


def extract_knf_metrics(knf_payload):
    vector = knf_payload.get("KNF_vector") or []
    metrics = {
        "SNCI": knf_payload.get("SNCI"),
        "SCDI": knf_payload.get("SCDI"),
        "SCDI_variance": knf_payload.get("SCDI_variance"),
    }
    for index in range(9):
        metrics[f"f{index + 1}"] = vector[index] if index < len(vector) else None
    return metrics


def metric_delta(current, reference):
    if current is None or reference is None:
        return None
    return current - reference


__all__ = [
    "KNF_METRIC_NAMES",
    "analyze_interface_geometry",
    "distance",
    "extract_knf_metrics",
    "metric_delta",
]
