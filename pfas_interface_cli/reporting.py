from .analysis import KNF_METRIC_NAMES


def format_value(value, decimals=6):
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.{decimals}f}"
    return str(value)


def make_table(headers, rows):
    widths = [len(str(header)) for header in headers]
    normalized_rows = []
    for row in rows:
        normalized = [str(cell) for cell in row]
        normalized_rows.append(normalized)
        for index, cell in enumerate(normalized):
            widths[index] = max(widths[index], len(cell))

    def render_row(row):
        return "| " + " | ".join(cell.ljust(widths[index]) for index, cell in enumerate(row)) + " |"

    separator = "|-" + "-|-".join("-" * width for width in widths) + "-|"
    lines = [render_row(headers), separator]
    lines.extend(render_row(row) for row in normalized_rows)
    return lines


def build_text_report(report):
    interface = report["interface_analysis"]
    closest = interface["closest_contact"]
    pfas_metrics = report["knf"]["pfas"]["metrics"]
    interface_metrics = report["knf"]["interface"]["metrics"]
    metric_names = KNF_METRIC_NAMES
    orientation_summary = report.get("orientation_summary") or {}
    orientation_results = report.get("orientations") or {}
    preferred_orientation = orientation_summary.get("preferred_orientation")
    evaluated_orientations = orientation_summary.get("evaluated") or list(orientation_results.keys())

    lines = [
        "PFAS Interface Workflow Summary",
        "================================",
        "",
        f"Run name: {report['run_name']}",
        f"Run directory: {report['run_directory']}",
        "",
        "Inputs",
        "------",
        f"Slab input: {report['inputs']['slab_xyz']}",
        f"PFAS input: {report['inputs']['pfas_xyz']}",
        "",
        "Settings",
        "--------",
        "",
    ]
    lines.extend(
        make_table(
            ["Setting", "Value"],
            [
                ["GFN level", str(report["settings"]["gfn"])],
                ["PFAS/interface charge", str(report["settings"]["pfas_charge"])],
                ["Spin multiplicity", str(report["settings"]["multiplicity"])],
                ["Orientation mode", str(report["settings"].get("orientation_mode", "single"))],
                ["Preferred orientation", preferred_orientation or "n/a"],
                ["Placement gap", f"{report['settings']['placement_gap_angstrom']:.3f} A"],
            ],
        )
    )
    lines.extend(
        [
            "",
            "Atom Ranges",
            "-----------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Group", "Atom Range"],
            [
                [
                    "Slab atoms in interface",
                    f"{report['atom_ranges']['slab_atoms'][0]}-{report['atom_ranges']['slab_atoms'][1]}",
                ],
                [
                    "PFAS atoms in interface",
                    f"{report['atom_ranges']['pfas_atoms'][0]}-{report['atom_ranges']['pfas_atoms'][1]}",
                ],
            ],
        )
    )
    lines.extend(
        [
            "",
            "Energies (Eh)",
            "-------------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Quantity", "Value (Eh)"],
            [
                ["E_slab", f"{report['energies_eh']['slab']:.12f}"],
                ["E_PFAS", f"{report['energies_eh']['pfas']:.12f}"],
                ["E_interface", f"{report['energies_eh']['interface']:.12f}"],
                ["Delta_E_ads", f"{report['energies_eh']['delta_e_ads']:.12f}"],
            ],
        )
    )
    if orientation_results:
        lines.extend(
            [
                "",
                "Orientation Comparison",
                "----------------------",
                "",
            ]
        )
        orientation_rows = []
        ordered_names = [name for name in evaluated_orientations if name in orientation_results]
        for orientation_name in ordered_names:
            payload = orientation_results[orientation_name]
            orientation_rows.append(
                [
                    orientation_name,
                    f"{payload['energies_eh']['interface']:.12f}",
                    f"{payload['energies_eh']['delta_e_ads']:.12f}",
                    payload["interface_analysis"]["contact_classification"],
                ]
            )
        lines.extend(
            make_table(
                [
                    "Orientation",
                    "E_interface (Eh)",
                    "Delta_E_ads (Eh)",
                    "Contact interpretation",
                ],
                orientation_rows,
            )
        )
        anisotropy = orientation_summary.get("delta_delta_e_ads_parallel_minus_perpendicular")
        if anisotropy is not None:
            lines.extend(
                [
                    "",
                    f"Delta_Delta_E_ads (parallel - perpendicular): {anisotropy:.12f} Eh",
                ]
            )
    lines.extend(
        [
            "",
            "Geometry",
            "--------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Metric", "Value"],
            [
                [
                    "Slab Z range",
                    f"{interface['slab_bounds']['z_min']:.4f} to {interface['slab_bounds']['z_max']:.4f}",
                ],
                ["Slab thickness", f"{interface['slab_thickness_angstrom']:.4f} A"],
                [
                    "Interface total Z range",
                    f"{interface['total_bounds']['z_min']:.4f} to {interface['total_bounds']['z_max']:.4f}",
                ],
                [
                    "PFAS Z range",
                    f"{interface['pfas_bounds']['z_min']:.4f} to {interface['pfas_bounds']['z_max']:.4f}",
                ],
                [
                    "PFAS offset from slab top",
                    (
                        f"{interface['pfas_offset_from_slab_top']['z_min_minus_slab_top']:.4f} "
                        f"to {interface['pfas_offset_from_slab_top']['z_max_minus_slab_top']:.4f} A"
                    ),
                ],
                ["Interface classification", interface["interface_classification"]],
            ],
        )
    )
    lines.extend(
        [
            "",
            "Closest PFAS-to-Water Contact",
            "-----------------------------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Metric", "Value"],
            [
                [
                    "Elements checked on PFAS",
                    ", ".join(interface["contact_elements_used"]) if interface["contact_elements_used"] else "none",
                ],
                ["Minimum distance", "n/a" if closest is None else f"{closest['distance_angstrom']:.4f} A"],
                [
                    "Closest pair",
                    (
                        "n/a"
                        if closest is None
                        else (
                            f"PFAS atom {closest['pfas_atom']['index']} ({closest['pfas_atom']['element']}) "
                            f"to water atom {closest['water_atom']['index']} ({closest['water_atom']['element']})"
                        )
                    ),
                ],
                ["Contact interpretation", interface["contact_classification"]],
            ],
        )
    )
    lines.extend(
        [
            "",
            "KNF Metrics",
            "-----------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Metric", "PFAS", "Interface", "Shift (interface - PFAS)"],
            [
                [
                    name,
                    format_value(pfas_metrics.get(name)),
                    format_value(interface_metrics.get(name)),
                    format_value(report["knf"]["delta_interface_minus_pfas"].get(name)),
                ]
                for name in metric_names
            ],
        )
    )
    lines.extend(
        [
            "",
            "Key Output Files",
            "----------------",
            "",
        ]
    )
    lines.extend(
        make_table(
            ["Artifact", "Path"],
            [
                ["Sorted slab", report["outputs"]["slab_sorted_xyz"]],
                ["Optimized slab", report["outputs"]["slab_optimized_xyz"]],
                ["Optimized PFAS", report["outputs"]["pfas_optimized_xyz"]],
                ["Merged interface start", report["outputs"]["interface_initial_xyz"]],
                ["Optimized interface", report["outputs"]["interface_optimized_xyz"]],
                ["Summary JSON", report["outputs"]["summary_json"]],
                ["Summary TXT", report["outputs"]["summary_txt"]],
                ["KNF PFAS JSON", report["outputs"]["knf_pfas_json"]],
                ["KNF interface JSON", report["outputs"]["knf_interface_json"]],
            ],
        )
    )
    orientation_outputs = report["outputs"].get("orientation_outputs")
    if orientation_outputs:
        lines.extend(
            [
                "",
                "Orientation Output Files",
                "------------------------",
                "",
            ]
        )
        orientation_output_rows = []
        ordered_names = [name for name in evaluated_orientations if name in orientation_outputs]
        for orientation_name in ordered_names:
            payload = orientation_outputs[orientation_name]
            orientation_output_rows.append(
                [
                    orientation_name,
                    payload["interface_initial_xyz"],
                    payload["interface_optimized_xyz"],
                    payload["knf_interface_json"],
                ]
            )
        lines.extend(
            make_table(
                ["Orientation", "Interface start XYZ", "Optimized interface XYZ", "KNF interface JSON"],
                orientation_output_rows,
            )
        )
    return "\n".join(lines) + "\n"


__all__ = ["build_text_report", "format_value", "make_table"]
