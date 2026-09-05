#!/usr/bin/env python3
"""Analyze one ParityFanout draw applied to the shortened BCH250 code.

The computation is diagnostic.  It applies the exact hypergeometric
ParityFanout transition to two source profiles:

1. a random-even profile restricted to the proved weight window [38, 218];
2. the code-independent constant-weight packing envelope plus exact code mass.

The first profile is a model.  For each output shell, the second calculation
maximizes the layered transition over all source spectra with exact mass
2^125-1 and the packing caps.  Different output shells can have different
maximizers.  The calculation does not prove concentration for one sampled
fanout composition.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def fanout_transition(
    length: int,
    input_weight: int,
    source_size: int,
    target_size: int,
) -> list[float]:
    """Return Pr[wt(F_{S,T}(x))=h] for a fixed weight input x."""

    result = [0.0] * (length + 1)
    source_denominator = math.comb(length, source_size)
    for overlap in range(source_size + 1):
        if overlap > input_weight:
            continue
        if source_size - overlap > length - input_weight:
            continue
        source_probability = (
            math.comb(input_weight, overlap)
            * math.comb(length - input_weight, source_size - overlap)
            / source_denominator
        )
        if overlap % 2 == 0:
            result[input_weight] += source_probability
            continue

        remaining_ones = input_weight - overlap
        remaining_length = length - source_size
        target_denominator = math.comb(remaining_length, target_size)
        for target_overlap in range(target_size + 1):
            if target_overlap > remaining_ones:
                continue
            if target_size - target_overlap > remaining_length - remaining_ones:
                continue
            target_probability = (
                math.comb(remaining_ones, target_overlap)
                * math.comb(
                    remaining_length - remaining_ones,
                    target_size - target_overlap,
                )
                / target_denominator
            )
            output_weight = input_weight + target_size - 2 * target_overlap
            result[output_weight] += source_probability * target_probability

    error = abs(sum(result) - 1.0)
    if error > 2e-14:
        raise RuntimeError(
            f"fanout transition at weight {input_weight} has mass error {error}"
        )
    return result


def random_even_window_profile(
    length: int,
    dimension: int,
    minimum_weight: int,
    maximum_weight: int,
) -> list[float]:
    result = [0.0] * (length + 1)
    permitted = [
        weight
        for weight in range(minimum_weight, maximum_weight + 1)
        if weight % 2 == 0
    ]
    denominator = sum(math.comb(length, weight) for weight in permitted)
    total_nonzero = float((1 << dimension) - 1)
    for weight in permitted:
        result[weight] = total_nonzero * math.comb(length, weight) / denominator
    return result


def load_packing_envelope(path: Path, length: int) -> list[float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result = [0.0] * (length + 1)
    for row in payload["spectrum"]:
        weight = int(row["weight"])
        if weight == 0:
            continue
        result[weight] = float(int(row["multiplicity_upper"]))
    return result


def apply_transition(
    spectrum: list[float], transitions: list[list[float]]
) -> list[float]:
    length = len(spectrum) - 1
    result = [0.0] * (length + 1)
    for input_weight, multiplicity in enumerate(spectrum):
        if multiplicity == 0.0:
            continue
        for output_weight, probability in enumerate(transitions[input_weight]):
            result[output_weight] += multiplicity * probability
    return result


def mass_coupled_output_envelope(
    source_caps: list[float],
    layered_transition: np.ndarray,
    total_mass: float,
) -> list[float]:
    """Maximize each output shell subject to caps and exact source mass."""

    active_source_weights = [
        weight for weight, cap in enumerate(source_caps) if cap > 0.0
    ]
    if sum(source_caps) < total_mass:
        raise ValueError("source shell caps do not contain the exact code mass")
    result = [0.0] * layered_transition.shape[1]
    for output_weight in range(layered_transition.shape[1]):
        order = sorted(
            active_source_weights,
            key=lambda weight: float(layered_transition[weight, output_weight]),
            reverse=True,
        )
        remaining = total_mass
        objective = 0.0
        for source_weight in order:
            allocated = min(source_caps[source_weight], remaining)
            objective += allocated * float(
                layered_transition[source_weight, output_weight]
            )
            remaining -= allocated
            if remaining <= 0.0:
                break
        if remaining > total_mass * 2e-15:
            raise RuntimeError("mass-coupled shell optimizer did not fill the code")
        result[output_weight] = objective
    return result


def spectrum_payload(
    candidate: str,
    spectrum: list[float],
    parameters: dict[str, object],
    scope: str,
) -> dict[str, object]:
    rows: list[dict[str, int | float]] = [
        {"weight": 0, "log2_expected_multiplicity": 0.0}
    ]
    for weight, multiplicity in enumerate(spectrum):
        if multiplicity <= 0.0:
            continue
        rows.append(
            {
                "weight": weight,
                "expected_multiplicity": multiplicity,
                "log2_expected_multiplicity": math.log2(multiplicity),
            }
        )
    return {
        "schema": "one-sampled-bch250-parityfanout-spectrum-v1",
        "candidate": candidate,
        "parameters": parameters,
        "spectrum": rows,
        "scope": scope,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--packing-envelope",
        type=Path,
        default=Path(__file__).with_name(
            "shortened_bch250_125_constant_weight_envelope.json"
        ),
    )
    parser.add_argument("--source-size", type=int, default=31)
    parser.add_argument("--target-size", type=int, default=33)
    parser.add_argument("--layers", type=int, default=1)
    parser.add_argument("--dimension", type=int, default=125)
    parser.add_argument(
        "--independent-row-local",
        action="store_true",
        help="label the one-row expectation for independent wrappers across rows",
    )
    parser.add_argument("--model-output", type=Path, required=True)
    parser.add_argument("--envelope-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    args = parser.parse_args()
    if args.layers <= 0:
        parser.error("--layers must be positive")
    if not 1 <= args.dimension <= 125:
        parser.error("--dimension must lie in [1,125]")

    length = 250
    dimension = args.dimension
    minimum_weight = 38
    maximum_weight = 218
    total_nonzero = float((1 << dimension) - 1)

    transitions = [
        fanout_transition(
            length,
            weight,
            args.source_size,
            args.target_size,
        )
        for weight in range(length + 1)
    ]
    modeled_source = random_even_window_profile(
        length,
        dimension,
        minimum_weight,
        maximum_weight,
    )
    envelope_source = load_packing_envelope(args.packing_envelope, length)
    transition_matrix = np.asarray(transitions, dtype=np.float64)
    layered_transition = np.linalg.matrix_power(transition_matrix, args.layers)
    modeled_output = list(
        np.asarray(modeled_source, dtype=np.float64) @ layered_transition
    )
    raw_envelope_output = list(
        np.asarray(envelope_source, dtype=np.float64) @ layered_transition
    )
    envelope_output = mass_coupled_output_envelope(
        envelope_source,
        layered_transition,
        total_nonzero,
    )

    parameters: dict[str, object] = {
        "outer_bits": length,
        "outer_dimension": dimension,
        "source_size": args.source_size,
        "target_size": args.target_size,
        "fanout_layers": args.layers,
        "source_minimum_weight": minimum_weight,
        "source_maximum_weight": maximum_weight,
        "fanout_sampled_once_and_reused": not args.independent_row_local,
        "fanout_wrappers_independent_across_rows": args.independent_row_local,
    }
    model_payload = spectrum_payload(
        "ShortenedBCH250x125-ParityFanout31x33 random-even model",
        modeled_output,
        parameters,
        (
            "Expected spectrum after one ParityFanout draw under a modeled "
            "random-even source profile restricted to weights 38 through 218. "
            "The fanout transition is exact in form. The source spectrum and "
            "binary64 arithmetic are diagnostic. This is not a concentration "
            "statement for one sampled and reused fanout map."
        ),
    )
    envelope_payload = spectrum_payload(
        "ShortenedBCH250x125-ParityFanout31x33 packing envelope",
        envelope_output,
        parameters,
        (
            "Expected per-shell upper bounds obtained by maximizing the layered "
            "fanout transition over every nonnegative source spectrum with exact "
            "total mass 2^125-1 and the source packing caps. Each output shell "
            "uses a separate continuous linear-program optimum, so the resulting "
            "output bounds are not jointly realizable in general. They do not "
            "prove concentration for one sampled and reused fanout composition."
        ),
    )

    comparison_rows = []
    for weight in range(1, length + 1):
        modeled = modeled_output[weight]
        upper = envelope_output[weight]
        if modeled <= 0.0 and upper <= 0.0:
            continue
        gap = math.inf if modeled <= 0.0 else math.log2(upper / modeled)
        comparison_rows.append(
            {
                "weight": weight,
                "modeled_log2_multiplicity": (
                    -math.inf if modeled <= 0.0 else math.log2(modeled)
                ),
                "envelope_log2_multiplicity": (
                    -math.inf if upper <= 0.0 else math.log2(upper)
                ),
                "envelope_minus_model_bits": gap,
            }
        )
    finite_rows = [
        row
        for row in comparison_rows
        if math.isfinite(float(row["envelope_minus_model_bits"]))
    ]
    dominant_gap = max(
        finite_rows,
        key=lambda row: float(row["envelope_minus_model_bits"]),
    )
    summary = {
        "schema": "one-sampled-bch250-parityfanout-comparison-v1",
        "parameters": parameters,
        "checks": {
            "modeled_source_log2_mass": math.log2(sum(modeled_source)),
            "modeled_output_log2_mass": math.log2(sum(modeled_output)),
            "raw_envelope_source_log2_mass": math.log2(sum(envelope_source)),
            "raw_envelope_output_log2_mass": math.log2(
                sum(raw_envelope_output)
            ),
            "mass_coupled_envelope_output_log2_mass": math.log2(
                sum(envelope_output)
            ),
        },
        "dominant_pointwise_gap": dominant_gap,
        "comparison_rows": comparison_rows,
        "scope": (
            "Diagnostic comparison of a modeled spectrum and a rigorous but "
            "very loose pointwise source envelope. It proves neither a good-"
            "spectrum event nor a finite distance statement."
        ),
    }

    for path, payload in (
        (args.model_output, model_payload),
        (args.envelope_output, envelope_payload),
        (args.summary_output, summary),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote,{path}")
    print(
        "dominant_gap_bits,"
        f"{dominant_gap['envelope_minus_model_bits']:.9f},"
        "weight,"
        f"{dominant_gap['weight']}"
    )


if __name__ == "__main__":
    main()
