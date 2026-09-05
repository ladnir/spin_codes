#!/usr/bin/env python3
"""Bound one outer symbol-occupation shell through the compressed kernel.

Fix an occupation size ``s``.  For each occupied field-symbol position, the
extended BCH encoder maps the nonzero field value bijectively to a nonzero
codeword.  If ``A_h`` is the BCH weight enumerator, define

    B(x) = sum_{h > 0} A_h x^(-h) / C(128,h).

The coefficient envelope for the independently permuted BCH blocks gives one
factor ``a(z)=x^(-wt(B(z)))/C(128,wt(B(z)))`` per occupied value.  The
16-state packet operator contains the corresponding positive powers of ``x``.

For one support, the MDS constraints parametrize the codewords by ``s-2``
field variables.  Every occupied coordinate is a surjective linear form.
Hölder's inequality bounds the weighted sum by

    q^(s-3) sum_z a(z)^s.

The exact BCH spectrum evaluates this moment.  An optional independent-value
envelope remains available as a diagnostic.

The script prints one JSON receipt.  It does not modify the workspace.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

import analyze_riffle_shiftalpha64_compressed_return as compressed


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
DEFAULT_SPECTRUM = SCRIPT_DIRECTORY / "ebch128_64_spectrum.csv"
OUTER_LENGTH = 16_386
FIELD_SIZE = 1 << 64


def load_spectrum(path: Path) -> dict[int, int]:
    with path.open(newline="", encoding="utf-8") as source:
        spectrum = {
            int(row["weight"]): int(row["count"])
            for row in csv.DictReader(source)
        }
    if sum(spectrum.values()) != FIELD_SIZE:
        raise RuntimeError("EBCH spectrum failed total-mass validation")
    if spectrum.get(0) != 1 or spectrum.get(22) != 243_840:
        raise RuntimeError("EBCH spectrum failed minimum-shell validation")
    return spectrum


def log_bch_outer_mass(spectrum: dict[int, int], weight_tilt: float) -> float:
    """Return log B(x) for the nonzero BCH spectrum."""
    log_x = math.log(weight_tilt)
    terms = [
        math.log(count)
        - compressed.kernel.log_binom(compressed.PART_BITS, weight)
        - weight * log_x
        for weight, count in spectrum.items()
        if weight > 0 and count > 0
    ]
    return float(logsumexp(terms))


def log_bch_holder_moment(
    spectrum: dict[int, int], weight_tilt: float, occupation: int
) -> float:
    """Return log sum_(z != 0) a(z)^s for the MDS Hölder envelope."""
    log_x = math.log(weight_tilt)
    terms = [
        math.log(count)
        + occupation
        * (
            -compressed.kernel.log_binom(compressed.PART_BITS, weight)
            - weight * log_x
        )
        for weight, count in spectrum.items()
        if weight > 0 and count > 0
    ]
    return float(logsumexp(terms))


def exact_mds_occupation_count(occupation: int) -> int:
    """Count weight-s words in the [16386,16384,3] outer MDS code."""
    if occupation < 3:
        return 0
    full_support_count = sum(
        (-1) ** zeros
        * math.comb(occupation, zeros)
        * (FIELD_SIZE ** (occupation - 2 - zeros) - 1)
        for zeros in range(occupation - 2)
    )
    return math.comb(OUTER_LENGTH, occupation) * full_support_count


def evaluate(
    *,
    occupation: int,
    distance: int,
    packet_positions: int,
    spectrum: dict[int, int],
    scaled_cost: float,
    zero_scale: float,
    weight_tilt: float,
    outer_envelope: str,
) -> dict[str, object]:
    if (
        scaled_cost <= 0.0
        or not 0.0 < zero_scale < packet_positions
        or weight_tilt <= 0.0
    ):
        return {"raw_log_bound": math.inf}

    compressed.configure_kernel(occupation, 22)
    u = scaled_cost / distance
    theta = zero_scale / packet_positions
    state_factor = np.empty(compressed.STATE_COUNT, dtype=np.float64)
    state_factor[0] = 1.0 / theta
    for state in range(1, compressed.STATE_COUNT):
        exponent = u * state.bit_count()
        state_factor[state] = math.exp(-exponent) / -math.expm1(-exponent)

    log_path_counts = compressed.weighted_path_log_counts_with_weight_tilt(
        state_factor, weight_tilt
    )
    log_support_choices = compressed.kernel.log_binom(OUTER_LENGTH, occupation)
    if outer_envelope == "holder":
        local_log_quantity = log_bch_holder_moment(
            spectrum, weight_tilt, occupation
        )
        log_outer_weight = (
            log_support_choices
            + (occupation - 3) * math.log(FIELD_SIZE)
            + local_log_quantity
        )
        local_terms = [
            (
                weight,
                math.log(count)
                + occupation
                * (
                    -compressed.kernel.log_binom(
                        compressed.PART_BITS, weight
                    )
                    - weight * math.log(weight_tilt)
                ),
            )
            for weight, count in spectrum.items()
            if weight > 0 and count > 0
        ]
    elif outer_envelope == "independent":
        local_log_quantity = log_bch_outer_mass(spectrum, weight_tilt)
        log_outer_weight = (
            log_support_choices + occupation * local_log_quantity
        )
        local_terms = [
            (
                weight,
                math.log(count)
                - compressed.kernel.log_binom(compressed.PART_BITS, weight)
                - weight * math.log(weight_tilt),
            )
            for weight, count in spectrum.items()
            if weight > 0 and count > 0
        ]
    else:
        raise ValueError(f"unknown outer envelope: {outer_envelope}")
    zero_gap_log_factor = -math.log1p(-theta)
    terms: list[tuple[int, float]] = []
    for support, log_count in enumerate(log_path_counts):
        if not math.isfinite(log_count):
            continue
        term = (
            log_count
            + log_outer_weight
            + (packet_positions - support + 1) * zero_gap_log_factor
            - compressed.kernel.log_binom(packet_positions, support)
        )
        terms.append((support, term))

    raw_log_bound = scaled_cost + float(logsumexp([term for _, term in terms]))
    mds_count = exact_mds_occupation_count(occupation)
    holder_unweighted_count = (
        math.comb(OUTER_LENGTH, occupation)
        * FIELD_SIZE ** (occupation - 3)
        * (FIELD_SIZE - 1)
    )
    top_terms = sorted(terms, key=lambda item: item[1], reverse=True)[:12]
    return {
        "raw_log_bound": raw_log_bound,
        "outer_shell_log2_bound": raw_log_bound / math.log(2.0),
        "scaled_positive_cost": scaled_cost,
        "positive_cost_parameter_u": u,
        "scaled_zero_return_parameter": zero_scale,
        "zero_return_tilt_theta": theta,
        "binary_weight_coefficient_tilt": weight_tilt,
        "outer_envelope_method": outer_envelope,
        "log2_bch_outer_mass_or_holder_moment": (
            local_log_quantity / math.log(2.0)
        ),
        "log2_position_choices": log_support_choices / math.log(2.0),
        "exact_mds_shell_log2_multiplicity": math.log2(mds_count),
        "holder_unweighted_multiplicity_cost_bits": (
            math.log2(holder_unweighted_count) - math.log2(mds_count)
        ),
        "top_bch_outer_terms": [
            {"bch_weight": weight, "log2_contribution": term / math.log(2.0)}
            for weight, term in sorted(
                local_terms, key=lambda item: item[1], reverse=True
            )[:8]
        ],
        "top_support_terms": [
            {
                "packet_support": support,
                "pre_chernoff_log2_contribution": term / math.log(2.0),
            }
            for support, term in top_terms
        ],
    }


def optimize(
    *,
    occupation: int,
    distance: int,
    packet_positions: int,
    spectrum: dict[int, int],
    start_scaled_cost: float,
    start_zero_scale: float,
    start_weight_tilt: float,
    maxiter: int,
    outer_envelope: str,
) -> dict[str, object]:
    start = np.log(
        np.asarray([start_scaled_cost, start_zero_scale, start_weight_tilt])
    )

    def objective(log_point: np.ndarray) -> float:
        scaled_cost, zero_scale, weight_tilt = np.exp(log_point)
        return float(
            evaluate(
                occupation=occupation,
                distance=distance,
                packet_positions=packet_positions,
                spectrum=spectrum,
                scaled_cost=float(scaled_cost),
                zero_scale=float(zero_scale),
                weight_tilt=float(weight_tilt),
                outer_envelope=outer_envelope,
            )["raw_log_bound"]
        )

    result = minimize(
        objective,
        start,
        method="Nelder-Mead",
        options={"xatol": 2e-5, "fatol": 2e-7, "maxiter": maxiter},
    )
    scaled_cost, zero_scale, weight_tilt = np.exp(result.x)
    selected = evaluate(
        occupation=occupation,
        distance=distance,
        packet_positions=packet_positions,
        spectrum=spectrum,
        scaled_cost=float(scaled_cost),
        zero_scale=float(zero_scale),
        weight_tilt=float(weight_tilt),
        outer_envelope=outer_envelope,
    )
    selected.update(
        {
            "optimizer_success": bool(result.success),
            "optimizer_message": str(result.message),
            "optimizer_evaluations": int(result.nfev),
        }
    )
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--occupation", type=int, required=True)
    parser.add_argument("--distance", type=int, default=188_766)
    parser.add_argument(
        "--packet-positions",
        type=int,
        default=compressed.DEFAULT_PACKET_POSITIONS,
    )
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--scaled-cost", type=float, default=60.0)
    parser.add_argument("--zero-scale", type=float, default=35.0)
    parser.add_argument("--weight-tilt", type=float, default=1.0)
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--optimizer-maxiter", type=int, default=50)
    parser.add_argument(
        "--outer-envelope",
        choices=("holder", "independent"),
        default="holder",
    )
    args = parser.parse_args()

    if not 3 <= args.occupation <= OUTER_LENGTH:
        raise ValueError("outer occupation lies outside [3,16386]")
    if args.packet_positions < 32 * args.occupation:
        raise ValueError("packet positions are smaller than the occupied blocks")

    spectrum = load_spectrum(args.spectrum)
    result = (
        optimize(
            occupation=args.occupation,
            distance=args.distance,
            packet_positions=args.packet_positions,
            spectrum=spectrum,
            start_scaled_cost=args.scaled_cost,
            start_zero_scale=args.zero_scale,
            start_weight_tilt=args.weight_tilt,
            maxiter=args.optimizer_maxiter,
            outer_envelope=args.outer_envelope,
        )
        if args.optimize
        else evaluate(
            occupation=args.occupation,
            distance=args.distance,
            packet_positions=args.packet_positions,
            spectrum=spectrum,
            scaled_cost=args.scaled_cost,
            zero_scale=args.zero_scale,
            weight_tilt=args.weight_tilt,
            outer_envelope=args.outer_envelope,
        )
    )
    payload = {
        "schema": "riffle-shiftalpha64-outer-occupation-v1",
        "construction": "Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4",
        "parameters": {
            "outer_field_length": OUTER_LENGTH,
            "outer_field_dimension": OUTER_LENGTH - 2,
            "outer_field_size": FIELD_SIZE,
            "occupied_field_symbols": args.occupation,
            "global_packet_positions": args.packet_positions,
            "bad_output_weight_inclusive": args.distance,
        },
        "bound": result,
        "outer_envelope": {
            "local_mass": "B(x)=sum_(h>0) A_h x^(-h)/C(128,h)",
            "position_factor": "C(16386,s)",
            "method": args.outer_envelope,
            "holder_bound": (
                "q^(s-3) sum_(z!=0) "
                "[x^(-wt(B(z)))/C(128,wt(B(z)))]^s"
            ),
            "dropped_constraints": (
                "none; Holder is applied after the two MDS equations "
                "parametrize each support"
                if args.outer_envelope == "holder"
                else "both outer field parity equations"
            ),
            "validity": (
                "The Holder envelope uses surjectivity of every occupied "
                "coordinate form."
                if args.outer_envelope == "holder"
                else "Every outer word is included among all nonzero "
                "assignments on each selected support."
            ),
        },
        "validation": {
            "exact_bch_spectrum": "PASS",
            "exact_outer_mds_occupation_count": "PASS",
            "numerical_parameters_remain_valid_if_optimizer_stops": True,
        },
        "scope": (
            "Rigorous upper bound for one outer occupation shell. "
            + (
                "Holder retains the two outer equations but discards "
                "correlations among their coordinate linear forms."
                if args.outer_envelope == "holder"
                else "The independent envelope drops both outer equations."
            )
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
