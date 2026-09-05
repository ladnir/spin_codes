#!/usr/bin/env python3
"""Compare symmetric BA conditioning windows at the finite B=240 point."""

from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import (  # noqa: E402
    DEFAULT_SELECTION,
    expected_ba_log_spectrum,
    load_envelope,
)
from diagnose_finite_k20_renyi_dense import (  # noqa: E402
    B,
    DISTANCE,
    EPOCHS,
    L,
    N,
    log_matrix_power_row_sum,
)


Q1_DIAGNOSTIC = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_k20_d11.json"
OUTPUT = WORKSTREAM / "golay_ba3_rm2sub_finite_B240_conditioning_window_diagnostic.json"
TARGET_MARGIN = 40.0


def endpoint_inner_log(envelope) -> tuple[float, float]:
    def objective(log_surprisal: float) -> float:
        surprisal = math.exp(log_surprisal)
        transfer, _ = envelope.transfer(
            candidate_probability=1.0,
            surprisal=surprisal,
        )
        return (
            log_matrix_power_row_sum(transfer, EPOCHS)
            + DISTANCE * surprisal
        )

    result = minimize_scalar(
        objective,
        bounds=(-4.0, 2.0),
        method="bounded",
        options={"xatol": 1e-13},
    )
    return float(result.fun), math.exp(float(result.x))


def main() -> None:
    spectrum = expected_ba_log_spectrum(B)
    q1 = json.loads(Q1_DIAGNOSTIC.read_text(encoding="utf-8"))
    q1_rows = {
        int(row["outer_weight"]): float(row["pointwise_log2_upper"])
        for row in q1["one_active"]["weight_rows"]
    }
    envelope = load_envelope(DEFAULT_SELECTION, DISTANCE / N)
    inner_log, surprisal = endpoint_inner_log(envelope)
    ideal_endpoint_log = L * math.log((1 << (B // 2)) - 1) + inner_log
    central_unconditioned_log = (
        L
        * (
            float(spectrum[B // 2])
            - (
                math.lgamma(B + 1)
                - 2.0 * math.lgamma(B // 2 + 1)
            )
            + B * math.log(2.0)
        )
        + inner_log
    )

    rows = []
    for lower in range(18, 26):
        upper = B - lower
        tail_logs = np.concatenate((spectrum[1:lower], spectrum[upper + 1 :]))
        tail = math.exp(float(logsumexp(tail_logs)))
        good = 1.0 - tail
        included_q1 = np.asarray(
            [q1_rows[weight] * math.log(2.0) for weight in range(lower, upper + 1)]
        )
        conditional_q1_log = float(logsumexp(included_q1)) - math.log(good)
        conditioned_endpoint_log = central_unconditioned_log - L * math.log(good)
        rows.append(
            {
                "permitted_weights": [lower, upper],
                "expected_tail_word_count": tail,
                "good_event_probability_lower": good,
                "conditioning_cost_bits_at_Q_equals_L": -L * math.log2(good),
                "q1_margin_bits": -conditional_q1_log / math.log(2.0),
                "central_shell_Q_equals_L_margin_bits": (
                    -conditioned_endpoint_log / math.log(2.0)
                ),
                "both_display_margins_above_40": (
                    -conditional_q1_log / math.log(2.0) > TARGET_MARGIN
                    and -conditioned_endpoint_log / math.log(2.0) > TARGET_MARGIN
                ),
            }
        )

    payload = {
        "schema": "golay-ba3-rm2sub-finite-b240-conditioning-window-diagnostic-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "target_margin_bits": TARGET_MARGIN,
        },
        "endpoint": {
            "bernoulli_half_surprisal": surprisal,
            "inner_log2_upper": inner_log / math.log(2.0),
            "uniform_random_block_margin_bits": -ideal_endpoint_log / math.log(2.0),
            "unconditioned_BA_central_shell_margin_bits": (
                -central_unconditioned_log / math.log(2.0)
            ),
        },
        "windows": rows,
        "recommended_window": [23, 217],
        "limitations": [
            "All displayed arithmetic and optimization use nearest binary64.",
            "The central-shell calculation is one dense composition, not the complete dense union.",
            "The Q=1 inputs come from the existing binary64 finite diagnostic.",
            "A final claim requires outward receipts for the chosen window and every occupation.",
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
