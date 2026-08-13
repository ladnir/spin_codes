#!/usr/bin/env python3
"""Fast independent regression gate for the g=8 conditioned-row outer.

The binary64 evaluator and the outward evaluator intentionally compute nearby,
but not identical, real expressions.  The outward path treats each
``exp(binary64_log_tilt)`` as an exact rational.  Accordingly, this checker
uses binary64 only as a conservative proximity diagnostic.  The exact outward
replay is frozen separately by a digest of every interval endpoint.

This gate performs no optimization and does not invoke the inner hardener.
"""

from __future__ import annotations

import hashlib
import json
import math
from decimal import Decimal
from pathlib import Path
from typing import Any

import numpy as np

from certify_packet_group_triangle_ledger import (
    conditioned_row_exact_graph_outer_outward,
)
from packet_group_outer_profile import K, atom_count
from packet_group_conditioned_row_high_precision import (
    evaluate_conditioned_row_high_precision,
)
from probe_packet_group_conditioned_row_outer import (
    evaluate_conditioned_outer,
    load_split_spectrum,
)


ROOT = Path(__file__).resolve().parents[1]
REGRESSION = ROOT / "scripts" / "packet_group_g8_conditioned_row_regression.json"
REGRESSION_SHA256 = "c2cf81a5b2f388fc5afd38cec14dbd0913a0dfdaaae438c57eeb9db14c5a936b"
OUTWARD_REPLAY_SHA256 = "3f038cdab9d6ee64867aece01baccab72fc544d69e0a8f3ce1d43d41b15e4de2"
PROXIMITY_BITS = 1e-6

G4_LEADER = {
    "name": "g4_asymmetric_leader_row132",
    "profile": [429359, 24531, 38395, 28039, 3964],
    "log_variables": [
        0.0,
        -2.8922414207312634,
        -2.614894602879822,
        -2.4693193677304146,
        -3.0107049207087875,
    ],
    "band_coefficients": [0.4851090410068871, 0.5148909589931129, 0.999],
    "pair_cauchy_theta": 0.38277584030042644,
    "outer_log2": 198830.2490738372,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"g=8 conditioned-row outward regression: {message}")


def outward_profile_interval(row: dict[str, Any], group_bits: int):
    coefficient = row.get("band_coefficients", row.get("band1_coefficient"))
    constant, charges, report = conditioned_row_exact_graph_outer_outward(
        group_bits,
        row["log_variables"],
        coefficient,
        row["pair_cauchy_theta"],
    )
    require(
        len(row["profile"]) == group_bits + 1 == len(charges),
        f"{row['name']}: outward dimension mismatch",
    )
    value = constant
    for count, charge in zip(row["profile"], charges):
        value = value - charge.times_int(int(count))
    require(report["group_bits"] == group_bits, f"{row['name']}: wrong metadata")
    return value


def binary64_profile_value(
    row: dict[str, Any],
    group_bits: int,
    spectrum01: np.ndarray,
    punctured01: np.ndarray,
    spectrum12: np.ndarray,
) -> float:
    coefficient = row.get("band_coefficients")
    value, details = evaluate_conditioned_outer(
        row["profile"],
        np.asarray(row["log_variables"], dtype=np.float64),
        float(row.get("band1_coefficient", 0.5)),
        float(row["pair_cauchy_theta"]),
        spectrum01,
        punctured01,
        spectrum12,
        band_coefficients=(
            None if coefficient is None else tuple(float(item) for item in coefficient)
        ),
        group_bits=group_bits,
    )
    require(details["group_bits"] == group_bits, f"{row['name']}: wrong binary64 metadata")
    return value


def endpoint_record(name: str, group_bits: int, interval) -> dict[str, Any]:
    return {
        "name": name,
        "group_bits": group_bits,
        "lo": str(interval.lo),
        "hi": str(interval.hi),
    }


def main() -> None:
    raw_fixture = REGRESSION.read_bytes()
    require(
        hashlib.sha256(raw_fixture).hexdigest() == REGRESSION_SHA256,
        "frozen g=8 fixture digest changed",
    )
    fixture = json.loads(raw_fixture)
    require(
        fixture.get("schema")
        == "permute-conv.packet-group-g8-conditioned-row-regression.v1",
        "frozen g=8 fixture schema changed",
    )

    spectrum01 = load_split_spectrum(ROOT / "out" / "ebch85_band01_split_spectrum.csv")
    punctured01 = load_split_spectrum(
        ROOT / "out" / "ebch84_punctured_band01_split_spectrum.csv",
        count_field="pair_count",
        divisor=42,
    )
    spectrum12 = load_split_spectrum(ROOT / "out" / "ebch86_band12_split_spectrum.csv")

    endpoint_records = []
    high_precision_checks = 0
    for row in fixture["rows"]:
        require(len(row["profile"]) == 9, f"{row['name']}: fixture is not g=8")
        outward = outward_profile_interval(row, group_bits=8)
        diagnostic = binary64_profile_value(
            row, 8, spectrum01, punctured01, spectrum12
        )
        require(math.isfinite(diagnostic), f"{row['name']}: non-finite binary64 value")
        require(
            abs(diagnostic - float(row["outer_log2"])) <= PROXIMITY_BITS,
            f"{row['name']}: binary64 fixture proximity failed",
        )
        require(
            abs(diagnostic - float(outward.lo)) <= PROXIMITY_BITS
            and abs(diagnostic - float(outward.hi)) <= PROXIMITY_BITS,
            f"{row['name']}: binary64/outward proximity failed",
        )
        reference, _masses = evaluate_conditioned_row_high_precision(
            8,
            row["profile"],
            row["log_variables"],
            row["band1_coefficient"],
            row["pair_cauchy_theta"],
            ROOT / "out" / "ebch85_band01_split_spectrum.csv",
            ROOT / "out" / "ebch86_band12_split_spectrum.csv",
            ROOT / "scripts" / "ebch128_graph24_spectrum.csv",
        )
        require(
            outward.lo <= reference <= outward.hi,
            f"{row['name']}: independent high-precision value escaped interval",
        )
        high_precision_checks += 1
        endpoint_records.append(endpoint_record(row["name"], 8, outward))

    zero_tilt = {
        "name": "g8_zero_tilt_total_mass",
        "profile": [atom_count(8)] + [0] * 8,
        "log_variables": [0.0] * 9,
        "band1_coefficient": 0.6,
        "pair_cauchy_theta": 0.5,
    }
    mass = outward_profile_interval(zero_tilt, group_bits=8)
    require(mass.lo <= Decimal(K) <= mass.hi, "zero-tilt interval does not contain K")
    endpoint_records.append(endpoint_record(zero_tilt["name"], 8, mass))

    g4_outward = outward_profile_interval(G4_LEADER, group_bits=4)
    g4_binary64 = binary64_profile_value(
        G4_LEADER, 4, spectrum01, punctured01, spectrum12
    )
    require(
        abs(g4_binary64 - G4_LEADER["outer_log2"]) <= PROXIMITY_BITS,
        "canonical g=4 asymmetric binary64 scalar changed",
    )
    require(
        abs(g4_binary64 - float(g4_outward.lo)) <= PROXIMITY_BITS
        and abs(g4_binary64 - float(g4_outward.hi)) <= PROXIMITY_BITS,
        "canonical g=4 asymmetric outward proximity failed",
    )
    endpoint_records.append(endpoint_record(G4_LEADER["name"], 4, g4_outward))

    replay_bytes = json.dumps(
        endpoint_records, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    replay_digest = hashlib.sha256(replay_bytes).hexdigest()
    require(
        replay_digest == OUTWARD_REPLAY_SHA256,
        f"outward endpoint digest changed: {replay_digest}",
    )

    print(f"g8_outward_vectors={len(fixture['rows'])}")
    print(f"independent_high_precision_checks={high_precision_checks}")
    print(f"g8_zero_tilt_mass_interval=[{mass.lo}, {mass.hi}]")
    print(f"g4_asymmetric_interval=[{g4_outward.lo}, {g4_outward.hi}]")
    print(f"outward_endpoint_sha256={replay_digest}")
    print("packet-group g=8 conditioned-row outward regression: PASS")


if __name__ == "__main__":
    main()
