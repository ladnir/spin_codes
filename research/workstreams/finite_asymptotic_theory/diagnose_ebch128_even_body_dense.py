#!/usr/bin/env python3
"""Finite dense diagnostic for repeated EBCH128 with parity-aware transfer.

The nonzero, non-all-one EBCH words are dominated pointwise by a scaled
uniform-even reference law.  A uniform-even row has 127 independent fair
coordinates; its final coordinate is their parity.  The diagnostic discards
the final transposed region and applies the finite RM2Sub transfer to the
remaining 127 regions.  Discarding output coordinates increases the Chernoff
moment because the tilt variable lies in (0, 1).

This script uses binary64 optimization.  It is not an outward certificate.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize


WORKSTREAM = Path(__file__).resolve().parent
REPOSITORY = WORKSTREAM.parents[1]
SCRIPTS = REPOSITORY / "scripts"
sys.path.insert(0, str(WORKSTREAM))
sys.path.insert(0, str(SCRIPTS))

from analyze_rm2sub_dense_occupation import (  # noqa: E402
    DEFAULT_SELECTION,
    DenseEnvelope,
    association_matrix,
    state_code_weights,
    verify_transpose_columns,
)
from import_wd_spectrum import parse_wd  # noqa: E402


DEFAULT_SPECTRUM = SCRIPTS / "EBCH128_64.wd"
DEFAULT_OUTPUT = WORKSTREAM / "ebch128_even_body_dense_k20_d11_diagnostic.json"


def log_choose(total: int, selected: int) -> float:
    if selected < 0 or selected > total:
        return -math.inf
    return (
        math.lgamma(total + 1)
        - math.lgamma(selected + 1)
        - math.lgamma(total - selected + 1)
    )


def log_binomial_pmf(total: int, selected: int, probability: float) -> float:
    return (
        log_choose(total, selected)
        + selected * math.log(probability)
        + (total - selected) * math.log1p(-probability)
    )


def body_envelope(path: Path, length: int) -> tuple[float, int]:
    """Return log C for A_w/C(n,w) <= C/2^(n-1), excluding 0 and n."""
    rows = parse_wd(path.read_text(encoding="utf-8"))
    best = -math.inf
    maximizing_weight = -1
    for weight, count in rows:
        if weight in (0, length) or count == 0:
            continue
        candidate = (
            math.log(count)
            - log_choose(length, weight)
            + (length - 1) * math.log(2.0)
        )
        if candidate > best:
            best = candidate
            maximizing_weight = weight
    return best, maximizing_weight


def log_transfer_moment(transfer: np.ndarray, epochs: int) -> float:
    eigenvalues = np.linalg.eigvals(transfer)
    radius = float(np.max(np.abs(eigenvalues)))
    if not radius > 0.0:
        return -math.inf
    normalized = transfer / radius
    powered = np.linalg.matrix_power(normalized, epochs)
    moment = float(np.sum(powered[0]))
    if not moment > 0.0:
        raise ArithmeticError("scaled transfer moment is not positive")
    return epochs * math.log(radius) + math.log(moment)


def finite_objective(
    *,
    envelope: DenseEnvelope,
    occupation: int,
    outer_blocks: int,
    body_log_envelope: float,
    retained_regions: int,
    epochs_per_region: int,
    distance: int,
    candidate_probability: float,
    surprisal: float,
    variant: str,
) -> tuple[float, dict[str, float]]:
    if variant == "three_state":
        transfer, details = envelope.transfer(
            candidate_probability=candidate_probability,
            surprisal=surprisal,
        )
    elif variant == "four_state":
        transfer, details = envelope.transfer_small_density(
            candidate_probability=candidate_probability,
            surprisal=surprisal,
        )
    else:
        raise ValueError(f"unknown transfer variant: {variant}")

    epochs = retained_regions * epochs_per_region
    inner_log = log_transfer_moment(transfer, epochs)
    conditioning_log = log_binomial_pmf(
        outer_blocks, occupation, candidate_probability
    )
    total = (
        log_choose(outer_blocks, occupation)
        + occupation * body_log_envelope
        - retained_regions * conditioning_log
        + distance * surprisal
        + inner_log
    )
    return total, {
        "candidate_probability": candidate_probability,
        "surprisal": surprisal,
        "z": math.exp(-surprisal),
        "conditioning_log2": conditioning_log / math.log(2.0),
        "inner_log2": inner_log / math.log(2.0),
        "perron_radius": float(details.get("perron_radius", 0.0)),
    }


def optimize_occupation(
    *,
    envelope: DenseEnvelope,
    occupation: int,
    outer_blocks: int,
    body_log_envelope: float,
    retained_regions: int,
    epochs_per_region: int,
    distance: int,
) -> dict[str, object]:
    alpha = occupation / outer_blocks
    epsilon = 1e-10
    starts = [
        (min(1.0 - epsilon, max(epsilon, scale * alpha)), max(epsilon, scale * alpha))
        for scale in (1.0, 1.2, 1.6, 2.0)
    ]
    starts.extend(((alpha, 0.2), (alpha, 0.8), (0.5, 1.0)))
    variant_rows = []
    for variant in ("three_state", "four_state"):
        best = None
        for start in starts:
            result = minimize(
                lambda point: finite_objective(
                    envelope=envelope,
                    occupation=occupation,
                    outer_blocks=outer_blocks,
                    body_log_envelope=body_log_envelope,
                    retained_regions=retained_regions,
                    epochs_per_region=epochs_per_region,
                    distance=distance,
                    candidate_probability=float(point[0]),
                    surprisal=float(point[1]),
                    variant=variant,
                )[0],
                x0=np.asarray(start, dtype=np.float64),
                bounds=((epsilon, 1.0 - epsilon), (epsilon, 8.0)),
                method="Nelder-Mead",
                options={"maxiter": 600, "xatol": 1e-11, "fatol": 1e-9},
            )
            if best is None or float(result.fun) < float(best.fun):
                best = result
        assert best is not None
        total, details = finite_objective(
            envelope=envelope,
            occupation=occupation,
            outer_blocks=outer_blocks,
            body_log_envelope=body_log_envelope,
            retained_regions=retained_regions,
            epochs_per_region=epochs_per_region,
            distance=distance,
            candidate_probability=float(best.x[0]),
            surprisal=float(best.x[1]),
            variant=variant,
        )
        variant_rows.append(
            {
                "variant": variant,
                "log2_upper": total / math.log(2.0),
                "margin_bits": -total / math.log(2.0),
                "optimizer_success": bool(best.success),
                **details,
            }
        )
    return {
        "occupation": occupation,
        "alpha": alpha,
        **min(variant_rows, key=lambda row: float(row["log2_upper"])),
        "variant_rows": variant_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--message-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-dimension", type=int, default=64)
    parser.add_argument("--outer-length", type=int, default=128)
    parser.add_argument("--relative-distance", type=float, default=0.11)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 12288, 16384),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    outer_blocks = args.message_bits // args.outer_dimension
    output_bits = outer_blocks * args.outer_length
    distance = math.floor(args.relative_distance * output_bits)
    if outer_blocks % 128:
        raise ValueError("128 must divide the number of outer blocks")
    epochs_per_region = outer_blocks // 128
    retained_regions = args.outer_length - 1

    selection = json.loads(args.selection.read_text(encoding="utf-8"))["selected"]
    generator_words = [int(value, 16) for value in selection["A_generator_words_hex"]]
    columns = [int(value, 16) for value in selection["B_columns_hex"]]
    verify_transpose_columns(
        generator_words=generator_words,
        columns=columns,
        output_bits=len(columns),
    )
    weights = state_code_weights(generator_words, len(columns))
    classes, counts, association = association_matrix(weights)
    envelope = DenseEnvelope(
        classes=classes,
        counts=counts,
        association=association,
        state_bits=len(generator_words),
        step_bits=len(columns),
        delta=args.relative_distance,
    )
    body_log, maximizing_weight = body_envelope(args.spectrum, args.outer_length)
    rows = []
    for occupation in args.occupations:
        if not 2 <= occupation <= outer_blocks:
            raise ValueError("occupations must lie in [2, outer_blocks]")
        row = optimize_occupation(
            envelope=envelope,
            occupation=occupation,
            outer_blocks=outer_blocks,
            body_log_envelope=body_log,
            retained_regions=retained_regions,
            epochs_per_region=epochs_per_region,
            distance=distance,
        )
        rows.append(row)
        print(
            f"occupation,{occupation},margin_bits,{row['margin_bits']:.9f},"
            f"variant,{row['variant']}",
            flush=True,
        )

    payload = {
        "schema": "ebch128-even-body-dense-k20-v1",
        "status": "BINARY64_DIAGNOSTIC",
        "parameters": {
            "message_bits": args.message_bits,
            "output_bits": output_bits,
            "distance": distance,
            "relative_distance": args.relative_distance,
            "outer_length": args.outer_length,
            "outer_dimension": args.outer_dimension,
            "outer_blocks": outer_blocks,
            "retained_regions": retained_regions,
            "discarded_parity_regions": 1,
            "epochs_per_region": epochs_per_region,
        },
        "outer_body_envelope": {
            "log2_counting_mass": body_log / math.log(2.0),
            "excess_over_rate_half_even_reference_bits": (
                body_log / math.log(2.0) - args.outer_dimension
            ),
            "maximizing_weight": maximizing_weight,
            "excluded_weights": [0, args.outer_length],
        },
        "occupation_rows": rows,
        "scope": (
            "The body calculation uses a pointwise envelope by the uniform-even "
            "row law. It discards one transposed region that carries the parity "
            "completion. Arithmetic and saddle optimization use nearest binary64. "
            "All-one outer words require a separate tail calculation."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
