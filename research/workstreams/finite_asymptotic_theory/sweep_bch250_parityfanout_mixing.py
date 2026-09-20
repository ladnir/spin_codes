#!/usr/bin/env python3
"""Sweep layered parity-fanout shapes for mass-coupled spectral mixing.

This is a nearest-binary64 diagnostic.  Every source spectrum in the sweep
has exact nonzero mass 2^125-1 and obeys the constant-weight packing caps for
the shortened BCH250 code.  For each output shell or tail band, a continuous
linear program maximizes the layered expected fanout mass over that source
polytope.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from analyze_one_sampled_bch250_parityfanout import (
    fanout_transition,
    load_packing_envelope,
    mass_coupled_output_envelope,
)


def maximize_linear_form(
    source_caps: list[float], coefficients: np.ndarray, total_mass: float
) -> float:
    weights = [weight for weight, cap in enumerate(source_caps) if cap > 0.0]
    weights.sort(key=lambda weight: float(coefficients[weight]), reverse=True)
    remaining = total_mass
    result = 0.0
    for weight in weights:
        allocated = min(source_caps[weight], remaining)
        result += allocated * float(coefficients[weight])
        remaining -= allocated
        if remaining <= 0.0:
            break
    if remaining > total_mass * 2e-15:
        raise RuntimeError("linear optimizer did not fill the exact code mass")
    return result


def bernoulli_envelope_bits(spectrum: list[float], endpoint_width: int) -> tuple[float, int]:
    length = len(spectrum) - 1
    best = -math.inf
    best_weight = -1
    for weight in range(1, length):
        if min(weight, length - weight) <= endpoint_width:
            continue
        multiplicity = spectrum[weight]
        if multiplicity <= 0.0:
            continue
        candidate = (
            math.log2(multiplicity)
            - math.log2(math.comb(length, weight))
            + length
        )
        if candidate > best:
            best = candidate
            best_weight = weight
    return best, best_weight


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--packing-envelope",
        type=Path,
        default=Path(__file__).with_name(
            "shortened_bch250_125_constant_weight_envelope.json"
        ),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    length = 250
    dimension = 125
    outer_blocks = 8448
    output_bits = 2_112_000
    random_dense_margin_bits = 8.399869558095574e-5 * output_bits
    allowed_excess_per_row = random_dense_margin_bits / outer_blocks
    total_mass = float((1 << dimension) - 1)
    source_caps = load_packing_envelope(args.packing_envelope, length)

    shapes = [(15, 17), (31, 33), (47, 49), (63, 65), (79, 81),
              (95, 97), (111, 113), (125, 125)]
    layer_counts = (1, 2, 4, 8, 16, 32, 64, 128)
    endpoint_widths = (
        16, 17, 18, 19, 20, 24, 28, 32, 34, 36, 38, 40, 60, 80, 90, 100
    )
    rows = []

    for source_size, target_size in shapes:
        transitions = np.asarray(
            [
                fanout_transition(length, weight, source_size, target_size)
                for weight in range(length + 1)
            ],
            dtype=np.float64,
        )
        layered = np.eye(length + 1, dtype=np.float64)
        previous_layers = 0
        for layers in layer_counts:
            for _ in range(layers - previous_layers):
                layered = layered @ transitions
            previous_layers = layers
            shell_envelope = mass_coupled_output_envelope(
                source_caps, layered, total_mass
            )
            endpoint_rows = []
            for endpoint_width in endpoint_widths:
                central_bits, central_weight = bernoulli_envelope_bits(
                    shell_envelope, endpoint_width
                )
                endpoint_probability = np.asarray(
                    [
                        sum(
                            layered[source_weight, output_weight]
                            for output_weight in range(length + 1)
                            if min(output_weight, length - output_weight)
                            <= endpoint_width
                        )
                        for source_weight in range(length + 1)
                    ],
                    dtype=np.float64,
                )
                maximum_endpoint_mass = maximize_linear_form(
                    source_caps, endpoint_probability, total_mass
                )
                endpoint_rows.append(
                    {
                        "endpoint_width": endpoint_width,
                        "central_envelope_bits": central_bits,
                        "central_excess_bits_per_row": central_bits - dimension,
                        "central_dominant_weight": central_weight,
                        "maximum_endpoint_mass_log2": math.log2(
                            maximum_endpoint_mass
                        ),
                        "maximum_endpoint_fraction_log2": math.log2(
                            maximum_endpoint_mass / total_mass
                        ),
                        "central_excess_fits_all_active_margin": (
                            central_bits - dimension <= allowed_excess_per_row
                        ),
                    }
                )
            rows.append(
                {
                    "source_size": source_size,
                    "target_size": target_size,
                    "layers": layers,
                    "xor_count_proxy_per_row": (
                        layers * (source_size + target_size - 1)
                    ),
                    "endpoint_rows": endpoint_rows,
                }
            )

    payload = {
        "schema": "bch250-parityfanout-mixing-sweep-v1",
        "parameters": {
            "outer_bits": length,
            "outer_dimension": dimension,
            "outer_blocks": outer_blocks,
            "output_bits": output_bits,
            "random_dense_margin_bits": random_dense_margin_bits,
            "allowed_central_excess_bits_per_row": allowed_excess_per_row,
            "packing_envelope": str(args.packing_envelope),
        },
        "rows": rows,
        "scope": (
            "Nearest-binary64 expected-spectrum diagnostic. The source shell "
            "polytope uses exact total mass and exact integer packing caps. Each "
            "output shell and endpoint band is maximized separately. The receipt "
            "does not prove concentration for one sampled and reused fanout "
            "composition, and its XOR count is only a scalar proxy."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    feasible = []
    for row in rows:
        for endpoint in row["endpoint_rows"]:
            if endpoint["central_excess_fits_all_active_margin"]:
                feasible.append(
                    (
                        row["xor_count_proxy_per_row"],
                        row["source_size"],
                        row["target_size"],
                        row["layers"],
                        endpoint["endpoint_width"],
                        endpoint["central_excess_bits_per_row"],
                    )
                )
    print(f"wrote,{args.output}")
    print(f"feasible_count,{len(feasible)}")
    for item in sorted(feasible)[:20]:
        print("feasible," + ",".join(str(value) for value in item))


if __name__ == "__main__":
    main()
