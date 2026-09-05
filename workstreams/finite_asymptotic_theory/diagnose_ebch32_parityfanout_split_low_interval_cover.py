"""Gap-free diagnostic refinement of the rejected low-band composition."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope  # noqa: E402
from diagnose_ebch32_ba_k20_three_band_qL import N, conditioned_spectrum  # noqa: E402
from diagnose_ebch32_parityfanout_ba_three_band_q import DISTANCE  # noqa: E402
from diagnose_ebch32_parityfanout_split_low_point import (  # noqa: E402
    BANDS,
    EXTERNAL_COUNTS,
    evaluate_fixed_witness,
    optimize_split,
)


if os.environ.get("SPIN_EBCH_PARITY_FANOUT", "0") != "1":
    raise RuntimeError("set SPIN_EBCH_PARITY_FANOUT=1")
if os.environ.get("SPIN_EBCH_LOWER_WEIGHT", "22") != "24":
    raise RuntimeError("set SPIN_EBCH_LOWER_WEIGHT=24")

OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_split_low_rejected_point_cover_d11.json"
CELL_CONTRIBUTION_MARGIN = 44.0


def test_interval(envelope, lower: int, upper: int, witness):
    checks = [
        evaluate_fixed_witness(envelope, endpoint, witness)
        for endpoint in sorted({lower, upper})
    ]
    length = upper - lower + 1
    minimum_margin = min(float(row["margin_bits"]) for row in checks)
    maximum_raw = max(float(row["raw_reference_bad_log2_upper"]) for row in checks)
    required_margin = CELL_CONTRIBUTION_MARGIN + math.log2(length)
    if maximum_raw >= 0.0 or minimum_margin <= required_margin:
        return None
    return {
        "lower_first_low_count": lower,
        "upper_first_low_count": upper,
        "integer_count": length,
        "minimum_endpoint_margin_bits": minimum_margin,
        "maximum_endpoint_raw_log2": maximum_raw,
        "log2_contribution_upper": math.log2(length) - minimum_margin,
        "endpoint_checks": checks,
        "reference_type_probabilities": witness["reference_type_probabilities"],
        "reference_value_probabilities": witness["reference_value_probabilities"],
        "reference_bit_probability": witness["reference_bit_probability"],
        "surprisal": witness["surprisal"],
        "holder_order": witness["holder_order"],
        "band_log2_renyi_moments": witness["band_log2_renyi_moments"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--maximum-processed", type=int, default=500)
    args = parser.parse_args()
    spectrum, _ = conditioned_spectrum()
    envelope = load_envelope(DEFAULT_SELECTION, DISTANCE / N)
    pending = [(0, EXTERNAL_COUNTS[0], None)]
    accepted = []
    processed = 0
    while pending and processed < args.maximum_processed:
        lower, upper, inherited = pending.pop()
        processed += 1
        receipt = (
            None
            if inherited is None
            else test_interval(envelope, lower, upper, inherited)
        )
        mode = "inherited-witness"
        witness = inherited
        if receipt is None:
            midpoint = (lower + upper) // 2
            witness = optimize_split(envelope, spectrum, midpoint)
            receipt = test_interval(envelope, lower, upper, witness)
            mode = "midpoint-witness"
        if receipt is not None:
            receipt["optimizer_mode"] = mode
            accepted.append(receipt)
        elif lower < upper:
            midpoint = (lower + upper) // 2
            pending.append((midpoint + 1, upper, witness))
            pending.append((lower, midpoint, witness))
        else:
            raise ArithmeticError(f"singleton split count {lower} did not certify")
        print(
            f"processed={processed},accepted={len(accepted)},pending={len(pending)},"
            f"interval=[{lower},{upper}]",
            flush=True,
        )
    logs = [float(row["log2_contribution_upper"]) for row in accepted]
    if logs:
        largest = max(logs)
        aggregate_log2 = largest + math.log2(
            sum(2.0 ** (value - largest) for value in logs)
        )
    else:
        aggregate_log2 = -math.inf
    payload = {
        "schema": "ebch32-parityfanout31x33-b256-split-low-rejected-point-cover-d11-v1",
        "status": (
            "COMPLETE_BINARY64_DIAGNOSTIC"
            if not pending
            else "INCOMPLETE_BINARY64_DIAGNOSTIC"
        ),
        "external_three_band_counts": list(EXTERNAL_COUNTS),
        "split_bands": [list(band) for band in BANDS],
        "covered_first_low_counts": [0, EXTERNAL_COUNTS[0]],
        "processed_intervals": processed,
        "accepted_intervals": len(accepted),
        "pending_intervals": [
            {"lower": lower, "upper": upper}
            for lower, upper, _ in pending
        ],
        "aggregate_log2_upper": aggregate_log2,
        "aggregate_margin_bits": -aggregate_log2,
        "intervals": accepted,
        "limitations": [
            "Witness discovery and displayed sums use nearest binary64 arithmetic.",
            "A complete diagnostic still requires an outward interval verifier.",
            "This receipt refines only the rejected external type (3344,1,21).",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: payload[key] for key in (
        "status",
        "processed_intervals",
        "accepted_intervals",
        "pending_intervals",
        "aggregate_margin_bits",
    )}, indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
