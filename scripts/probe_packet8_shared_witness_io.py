#!/usr/bin/env python3
"""Load diagnostic shared-drive fixed witnesses for profile-cover probes."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from probe_packet8_adaptive_simplex_atlas import split_cap_table
from probe_packet8_profile_simplex_landscape import Anchor, build_fixed_witness


def load_shared_witness_arrays(
    upgraded_cache: Path,
    report_paths: list[Path],
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Return names, constants, and charges from one cache plus retuned reports."""
    with np.load(upgraded_cache, allow_pickle=False) as cache:
        names = [str(value) for value in cache["names"]]
        constants = np.asarray(cache["constants"], dtype=np.float64).copy()
        charges = np.asarray(cache["charges"], dtype=np.float64).copy()

    if constants.ndim != 1 or charges.shape != (len(constants), 9):
        raise ValueError(f"invalid shared witness cache: {upgraded_cache}")
    if len(names) != len(constants):
        raise ValueError(f"shared witness name count mismatch: {upgraded_cache}")

    if not report_paths:
        return names, constants, charges

    split_caps = split_cap_table()
    extra_constants = []
    extra_charges = []
    for path in report_paths:
        report = json.loads(path.read_text(encoding="utf-8"))
        source = report["anchor"]
        anchor = Anchor(
            str(source["name"]),
            tuple(int(value) for value in report["profile"]),
            float(source["pole"]),
            tuple(float(value) for value in source["fugacities"]),
            float(source["outer_log_bound"]),
        )
        fixed = build_fixed_witness(anchor, split_caps)
        names.append(f"shared_report:{path.stem}")
        extra_constants.append(
            float(fixed.constant_log2)
            - float(report["shared_drive_improvement_log2"])
        )
        extra_charges.append(np.asarray(fixed.linear_charge, dtype=np.float64))

    return (
        names,
        np.concatenate((constants, np.asarray(extra_constants))),
        np.vstack((charges, np.vstack(extra_charges))),
    )
