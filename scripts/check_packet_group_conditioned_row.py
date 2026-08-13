#!/usr/bin/env python3
"""Regression checks for the packet-width conditioned-row evaluator."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from packet_group_drive_stratified import feasible_profile_count
from packet_group_outer_profile import K, N, atom_count
from probe_packet_group_conditioned_row_outer import (
    evaluate_conditioned_outer,
    load_split_spectrum,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WITNESS = ROOT / "out" / "g4_conditioned_row_residual_best140_asymmetric.json"
G8_REGRESSION = ROOT / "scripts" / "packet_group_g8_conditioned_row_regression.json"
DEFAULT_ROW = 132
EXPECTED_G4_OUTER_LOG2 = 198830.2490738372


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"conditioned-row regression: {message}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--witness", type=Path, default=DEFAULT_WITNESS)
    parser.add_argument("--row", type=int, default=DEFAULT_ROW)
    args = parser.parse_args()

    spectrum01 = load_split_spectrum(ROOT / "out" / "ebch85_band01_split_spectrum.csv")
    punctured01 = load_split_spectrum(
        ROOT / "out" / "ebch84_punctured_band01_split_spectrum.csv",
        count_field="pair_count",
        divisor=42,
    )
    spectrum12 = load_split_spectrum(ROOT / "out" / "ebch86_band12_split_spectrum.csv")

    artifact = json.loads(args.witness.read_text(encoding="utf-8"))
    row = artifact["rows"][args.row]
    require(int(row["group_bits"]) == 4, "frozen regression row is not g=4")
    details = row["outer_details"]
    common = (
        row["profile"],
        np.asarray(details["log_variables"], dtype=np.float64),
        float(details["band1_coefficient"]),
        float(details["pair_cauchy_theta"]),
        spectrum01,
        punctured01,
        spectrum12,
        1,
        tuple(float(value) for value in details["band_coefficients"]),
    )
    legacy, _legacy_details = evaluate_conditioned_outer(*common)
    generalized, generalized_details = evaluate_conditioned_outer(*common, group_bits=4)
    require(legacy == generalized, "explicit g=4 path is not bit-identical")
    require(legacy == float(row["outer_log2"]), "stored g=4 witness changed")
    require(legacy == EXPECTED_G4_OUTER_LOG2, "canonical g=4 scalar changed")
    require(generalized_details["group_bits"] == 4, "g=4 detail metadata changed")

    g8_profile = [atom_count(8)] + [0] * 8
    g8_value, g8_details = evaluate_conditioned_outer(
        g8_profile,
        np.zeros(9, dtype=np.float64),
        0.6,
        0.5,
        spectrum01,
        punctured01,
        spectrum12,
        group_bits=8,
    )
    require(math.isclose(g8_value, K, rel_tol=0.0, abs_tol=1e-6), "g=8 mass check failed")
    require(g8_details["group_bits"] == 8, "g=8 detail metadata changed")
    require(
        feasible_profile_count(2, N, 21) == 549757386632,
        "g=2 feasible profile count changed",
    )
    require(
        feasible_profile_count(4, N, 21) == 3148304370308993121588,
        "g=4 feasible profile count changed",
    )
    require(
        feasible_profile_count(8, N, 21)
        == 553169839211945865258921061892182603726,
        "g=8 feasible profile count changed",
    )

    g8_regression = json.loads(G8_REGRESSION.read_text(encoding="utf-8"))
    require(
        g8_regression.get("schema")
        == "permute-conv.packet-group-g8-conditioned-row-regression.v1",
        "g=8 regression schema changed",
    )
    for vector in g8_regression["rows"]:
        actual, _details = evaluate_conditioned_outer(
            vector["profile"],
            np.asarray(vector["log_variables"], dtype=np.float64),
            float(vector["band1_coefficient"]),
            float(vector["pair_cauchy_theta"]),
            spectrum01,
            punctured01,
            spectrum12,
            group_bits=8,
        )
        require(
            math.isclose(
                actual,
                float(vector["outer_log2"]),
                rel_tol=0.0,
                abs_tol=1e-9,
            ),
            f"g=8 frozen vector changed: {vector['name']}",
        )

    print(f"g4_outer_log2={legacy:.13f}")
    print(f"g8_total_mass_log2={g8_value:.13f}")
    print(f"g8_frozen_vectors={len(g8_regression['rows'])}")
    print("packet-group conditioned-row regression: PASS")


if __name__ == "__main__":
    main()
