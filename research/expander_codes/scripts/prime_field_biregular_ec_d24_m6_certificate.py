#!/usr/bin/env python3
"""Generate and verify the degree-24, memory-6 prime-field EC certificate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from ea_certificate import decimal_marker
from prime_field_biregular_ec_certificate import (
    full_support_markers_float,
    result_json,
)
from prime_field_biregular_ec_diagnostic import (
    candidate_parameters,
    singleton_constraint_balanced_ec_block_logbound,
)
from prime_field_biregular_ec_singleton_certificate import (
    EXACT_MARKERS as D26_M4_EXACT_MARKERS,
    SCHEMA,
    verify_certificate,
)
from prime_field_ea_diagnostic import gv_distance


PRIME = 2**127 - 1
DEGREE = 24
MEMORY = 6
K, N, REGION_LENGTH = candidate_parameters(2**21, DEGREE)
CUTOFF = math.floor(N * gv_distance(PRIME, 0.5))

EXACT_MARKERS: dict[tuple[int, int], dict[str, tuple[str, str]]] = {
    band: markers
    for band, markers in D26_M4_EXACT_MARKERS.items()
    if band[1] <= 192
}
EXACT_MARKERS.update({
    (193, 224): {
        "structural": ("0.9949515834398468", "0.0002005170754088813"),
        "field": ("0.9949515836172307", "0.00020145645463620567"),
    },
    (225, 256): {
        "structural": ("0.9942452684235772", "0.00022856329832220094"),
        "field": ("0.9942452670597987", "0.0002296147938680022"),
    },
    (257, 288): {
        "structural": ("0.9935470274615466", "0.00025642489918217864"),
        "field": ("0.9935470308392499", "0.0002573651977353216"),
    },
    (289, 320): {
        "structural": ("0.9928558638851366", "0.0002838222498269103"),
        "field": ("0.9928558267322602", "0.0002850059513677085"),
    },
    (321, 360): {
        "structural": ("0.992087245687715", "0.0003142039588116197"),
        "field": ("0.9920872490080946", "0.00031527046291951526"),
    },
    (361, 400): {
        "structural": ("0.9912397101061964", "0.00034798943340894935"),
        "field": ("0.9912397245291356", "0.0003488266315715851"),
    },
    (401, 500): {
        "structural": ("0.9897909175821032", "0.0004050392378688905"),
        "field": ("0.9897908823408765", "0.00040677752461040453"),
    },
    (501, 600): {
        "structural": ("0.9877296241331784", "0.00048750712617144865"),
        "field": ("0.9877296252749657", "0.0004879741162990873"),
    },
    (601, 700): {
        "structural": ("0.9855286934580273", "0.0005723752099763288"),
        "field": ("0.9855329029918858", "0.0005732364301116299"),
    },
})
EXACT_BANDS = tuple(EXACT_MARKERS)


def support_blocks(k: int) -> tuple[tuple[int, int], ...]:
    """Partition supports 701 through ``k - 1`` into convexity blocks."""
    blocks = [
        (701, 800),
        (801, 1_000),
        (1_001, 1_200),
        (1_201, 1_400),
        (1_401, 1_600),
        (1_601, 3_200),
        (3_201, 6_400),
        (6_401, 12_800),
        (12_801, 25_600),
        (25_601, 51_200),
        (51_201, 102_400),
        (102_401, 204_800),
        (204_801, 409_600),
        (409_601, 819_200),
    ]
    tail_widths = [2**power for power in range(16, 2, -1)] + [3]
    reserved = sum(tail_widths)
    start = blocks[-1][1] + 1
    first_width = (k - 1) - start + 1 - reserved
    if first_width <= 0:
        raise ValueError("support range is too short for the fixed partition")
    blocks.append((start, start + first_width - 1))
    start = blocks[-1][1] + 1
    for width in tail_widths:
        blocks.append((start, start + width - 1))
        start += width
    if start != k:
        raise AssertionError("support partition does not end at k - 1")
    return tuple(blocks)


BLOCKS = support_blocks(K)


def generate_certificate(*, target_bits: int, precision_bits: int) -> dict[str, Any]:
    exact_bands = []
    for lo, hi in EXACT_BANDS:
        selected = EXACT_MARKERS[(lo, hi)]
        exact_bands.append({
            "lo": lo,
            "hi": hi,
            "markers": {
                "structural": list(selected["structural"]),
                "field": list(selected["field"]),
            },
        })
        print(f"selected exact band [{lo},{hi}]", flush=True)

    blocks = []
    for lo, hi in BLOCKS:
        bound = singleton_constraint_balanced_ec_block_logbound(
            prime=PRIME,
            k=K,
            n=N,
            cutoff=CUTOFF,
            region_count=DEGREE,
            memory=MEMORY,
            support_start=lo,
            support_limit=hi,
        )
        blocks.append({
            "lo": lo,
            "hi": hi,
            "markers": {
                "structural": [
                    decimal_marker(value) for value in bound.structural_markers
                ],
                "field": [
                    decimal_marker(value) for value in bound.field_markers
                ],
            },
            "diagnostic_log2": bound.total_log2,
        })
        print(
            f"optimized block [{lo},{hi}]: {bound.total_log2:.6f} bits",
            flush=True,
        )

    return {
        "schema": SCHEMA,
        "parameters": {
            "prime": str(PRIME),
            "k": K,
            "n": N,
            "cutoff": CUTOFF,
            "region_count": DEGREE,
            "memory": MEMORY,
            "target_bits": target_bits,
        },
        "verification": {"precision_bits": precision_bits},
        "exact_bands": exact_bands,
        "blocks": blocks,
        "full_support": {
            "r": K,
            "markers": full_support_markers_float(
                prime=PRIME,
                n=N,
                cutoff=CUTOFF,
                memory=MEMORY,
                message_weight=K,
            ),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument("--target-bits", type=int, default=20)
    generate.add_argument("--precision-bits", type=int, default=256)
    verify = subparsers.add_parser("verify")
    verify.add_argument("certificate", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        certificate = generate_certificate(
            target_bits=args.target_bits,
            precision_bits=args.precision_bits,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(certificate, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.output}")
        return
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    result = verify_certificate(certificate)
    print(json.dumps(result_json(result), indent=2))
    if not result.success:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
