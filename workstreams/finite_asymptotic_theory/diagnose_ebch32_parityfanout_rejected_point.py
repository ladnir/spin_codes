"""Search saved dense-cell witnesses for the rejected lattice point."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_ebch32_ba_k20_three_band_qL import (  # noqa: E402
    N,
    conditioned_spectrum,
)
from diagnose_ebch32_parityfanout_ba_three_band_cover_all_q import (  # noqa: E402
    DISTANCE,
    certify_with_witness,
)


if os.environ.get("SPIN_EBCH_PARITY_FANOUT", "0") != "1":
    raise RuntimeError("set SPIN_EBCH_PARITY_FANOUT=1")
if os.environ.get("SPIN_EBCH_LOWER_WEIGHT", "22") != "24":
    raise RuntimeError("set SPIN_EBCH_LOWER_WEIGHT=24")

INPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_three_band_cover_all_q_d11.json"
POINT = np.asarray((3344.0, 1.0, 21.0), dtype=np.float64)


def box_distance(row: dict[str, object]) -> float:
    vertices = np.asarray(row["vertices_active_counts"], dtype=np.float64)
    lower = np.min(vertices, axis=0)
    upper = np.max(vertices, axis=0)
    displacement = np.maximum(lower - POINT, 0.0) + np.minimum(upper - POINT, 0.0)
    return float(np.dot(displacement, displacement))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=1000)
    args = parser.parse_args()
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    rows = [
        row
        for row in source["tetrahedra"]
        if "reference_type_probabilities" in row
    ]
    rows.sort(key=box_distance)
    spectrum, _ = conditioned_spectrum()
    envelope = load_envelope(DEFAULT_SELECTION, DISTANCE / N)
    point_vertices = np.repeat(POINT[None, :], 4, axis=0)
    for index, witness in enumerate(rows[: args.limit], 1):
        receipt = certify_with_witness(envelope, point_vertices, witness)
        if receipt is not None:
            print(
                json.dumps(
                    {
                        "status": "FOUND_DIAGNOSTIC_WITNESS",
                        "point": POINT.tolist(),
                        "tested_witnesses": index,
                        "source_mode": witness.get("optimizer_mode"),
                        "source_box_distance_squared": box_distance(witness),
                        "receipt": receipt,
                    },
                    indent=2,
                )
            )
            return
        if index % 100 == 0:
            print(f"tested={index}", flush=True)
    print(
        json.dumps(
            {
                "status": "NO_DIAGNOSTIC_WITNESS_FOUND",
                "point": POINT.tolist(),
                "tested_witnesses": min(args.limit, len(rows)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
