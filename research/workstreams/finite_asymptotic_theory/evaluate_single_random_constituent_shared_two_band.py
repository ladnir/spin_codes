#!/usr/bin/env python3
"""Shared-category two-band diagnostic for one random constituent.

The high Bernoulli tail is dominated by the low Bernoulli tail using
monotonicity of the RandomStepConv weight moment.  The two tails therefore
form one defect category with twice the tail majorant.  A two-dimensional
region recurrence retains the defect/central category assignment once across
all transposed regions.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import evaluate_ebch128_randomstepconv_g1 as transfer
from evaluate_single_random_constituent_highprob_bands import (
    band_majorant,
    log_multinomial,
)
from evaluate_single_random_constituent_highprob_renyi import (
    LOG2,
    spectrum_caps,
)


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_OUTPUT = WORKSTREAM / "single_random_constituent_B512_shared_two_band_q17_64_s22.json"


def exact_two_category_regions(
    zero: np.ndarray,
    defect: np.ndarray,
    central: np.ndarray,
    maximum_occupation: int,
    outer_rows: int,
) -> np.ndarray:
    """Average ordered products with fixed defect and central row counts."""

    side = maximum_occupation + 1
    current = np.full((side, side, 2, 2), -math.inf)
    current[0, 0] = transfer.log_identity()
    defects, centrals = np.indices((side, side))
    total = defects + centrals
    candidates = [transfer.log_entries(defect), transfer.log_entries(central)]
    log_zero = transfer.log_entries(zero)

    for completed in range(outer_rows):
        denominator = float(completed + 1)
        updated = np.full_like(current, -math.inf)

        inactive_weight = (completed + 1 - total) / denominator
        valid = (
            (total <= completed)
            & (total <= maximum_occupation)
            & (inactive_weight > 0.0)
        )
        term = transfer.log_matmul(current, log_zero)
        updated[valid] = term[valid] + np.log(inactive_weight[valid])[:, None, None]

        for axis, candidate in enumerate(candidates):
            source = [slice(None), slice(None)]
            target = [slice(None), slice(None)]
            source[axis] = slice(0, maximum_occupation)
            target[axis] = slice(1, maximum_occupation + 1)
            source_tuple = tuple(source)
            target_tuple = tuple(target)
            products = transfer.log_matmul(current[source_tuple], candidate)
            selected = (defects, centrals)[axis][target_tuple]
            source_totals = total[source_tuple]
            valid = (
                (source_totals <= completed)
                & (source_totals < maximum_occupation)
            )
            terms = np.full_like(products, -math.inf)
            terms[valid] = products[valid] + np.log(
                selected[valid] / denominator
            )[:, None, None]
            updated[target_tuple] = np.logaddexp(updated[target_tuple], terms)
        current = updated
    return current


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-occupation", type=int, default=17)
    parser.add_argument("--maximum-occupation", type=int, default=64)
    parser.add_argument("--block-bits", type=int, default=512)
    parser.add_argument("--dimension", type=int, default=256)
    parser.add_argument("--output-bits", type=int, default=1 << 21)
    parser.add_argument("--memory-bits", type=int, default=22)
    parser.add_argument("--distance-numerator", type=int, default=109)
    parser.add_argument("--distance-denominator", type=int, default=1000)
    parser.add_argument("--per-shell-failure-exponent", type=int, default=51)
    parser.add_argument("--low-upper", type=int, default=79)
    parser.add_argument("--central-upper", type=int, default=432)
    parser.add_argument("--grid-min", type=float, default=-8.0)
    parser.add_argument("--grid-max", type=float, default=-3.0)
    parser.add_argument("--grid-step", type=float, default=0.25)
    parser.add_argument("--prior", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if args.output_bits % args.block_bits:
        parser.error("block length must divide output length")
    if not 1 <= args.minimum_occupation <= args.maximum_occupation:
        parser.error("invalid occupation interval")
    block_bits = args.block_bits
    outer_rows = args.output_bits // block_bits
    if args.maximum_occupation > outer_rows:
        parser.error("maximum occupation exceeds the number of outer rows")
    distance_cutoff = (
        args.distance_numerator * args.output_bits
        + args.distance_denominator
        - 1
    ) // args.distance_denominator

    caps, event = spectrum_caps(
        math.ldexp(1.0, -args.per_shell_failure_exponent),
        block_bits=block_bits,
        dimension=args.dimension,
    )
    support = [weight for weight in range(1, block_bits + 1) if int(caps[weight])]
    low_limits = (min(support), args.low_upper)
    central_limits = (args.low_upper + 1, args.central_upper)
    high_limits = (args.central_upper + 1, max(support))
    low = band_majorant(caps, *low_limits, block_bits=block_bits)
    central = band_majorant(caps, *central_limits, block_bits=block_bits)
    high = band_majorant(caps, *high_limits, block_bits=block_bits)
    tail_log_majorant = max(float(low["log_majorant"]), float(high["log_majorant"]))
    defect_log_majorant = math.log(2.0) + tail_log_majorant
    defect_probability = min(
        float(low["value_probability"]),
        1.0 - float(high["value_probability"]),
    )

    compositions = [
        (defects, occupation - defects)
        for occupation in range(args.minimum_occupation, args.maximum_occupation + 1)
        for defects in range(occupation + 1)
    ]
    best = {composition: math.inf for composition in compositions}
    witnesses = {composition: math.nan for composition in compositions}
    if args.prior is not None:
        prior = json.loads(args.prior.read_text(encoding="utf-8"))
        for row in prior["rows"]:
            composition = (int(row["defect_rows"]), int(row["central_rows"]))
            if composition in best:
                best[composition] = float(row["pointwise_log2_upper"]) * LOG2
                witnesses[composition] = float(row["log_surprisal"])
    u_values = np.arange(
        args.grid_min,
        args.grid_max + args.grid_step / 2.0,
        args.grid_step,
    )
    for index, u in enumerate(u_values):
        surprisal = math.exp(float(u))
        z = math.exp(-surprisal)
        zero, active = transfer.step_matrices(z, args.memory_bits)
        defect_matrix = (1.0 - defect_probability) * zero + defect_probability * active
        central_probability = float(central["value_probability"])
        central_matrix = (1.0 - central_probability) * zero + central_probability * active
        regions = exact_two_category_regions(
            zero,
            defect_matrix,
            central_matrix,
            args.maximum_occupation,
            outer_rows,
        )
        for composition in compositions:
            powered = transfer.log_power(regions[composition], block_bits)
            moment = float(np.logaddexp(powered[0, 0], powered[0, 1]))
            defects, centrals = composition
            counts = (outer_rows - defects - centrals, defects, centrals)
            outer = (
                log_multinomial(counts)
                + defects * defect_log_majorant
                + centrals * float(central["log_majorant"])
            )
            value = outer + min(0.0, moment + distance_cutoff * surprisal)
            if value < best[composition]:
                best[composition] = value
                witnesses[composition] = float(u)
        print(f"tilt,{index + 1},{len(u_values)},u,{u:.6f}", flush=True)

    rows = [
        {
            "occupation": sum(composition),
            "defect_rows": composition[0],
            "central_rows": composition[1],
            "pointwise_log2_upper": best[composition] / LOG2,
            "margin_bits": -best[composition] / LOG2,
            "log_surprisal": witnesses[composition],
        }
        for composition in compositions
    ]
    aggregate = float(logsumexp(list(best.values())))
    payload = {
        "schema": "single-random-constituent-shared-two-band-v1",
        "status": "BINARY64_EXACT_SHARED_CATEGORY_DIAGNOSTIC",
        "parameters": {
            "outer_code": f"one uniform binary [{block_bits},{args.dimension}] generator",
            "outer_rows": outer_rows,
            "output_bits": args.output_bits,
            "distance_cutoff": distance_cutoff,
            "memory_bits": args.memory_bits,
            "occupations": [args.minimum_occupation, args.maximum_occupation],
        },
        "spectrum_event": event,
        "bands": {"low": low, "central": central, "high": high},
        "merged_defect": {
            "value_probability": defect_probability,
            "log2_majorant": defect_log_majorant / LOG2,
            "reason": "high-tail moment is dominated by the low-tail moment; summing the two symmetric tail labels costs one bit per defect row",
        },
        "claim": {
            "aggregate_log2_upper": aggregate / LOG2,
            "aggregate_margin_bits": -aggregate / LOG2,
            "worst_composition": max(rows, key=lambda row: float(row["pointwise_log2_upper"])),
        },
        "rows": rows,
        "limitations": [
            "Nearest binary64 arithmetic is not an outward certificate.",
            "Only the stated occupation interval is covered.",
            "The RandomStepConv monotonicity lemma must be included in the final written proof and outward verifier.",
        ],
        "prior_receipt": None if args.prior is None else str(args.prior),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
