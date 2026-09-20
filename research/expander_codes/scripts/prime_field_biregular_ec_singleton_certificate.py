#!/usr/bin/env python3
"""Generate and verify singleton-trace finite-field EC certificates."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from flint import arb, ctx

from ea_certificate import decimal_marker
from prime_field_biregular_ec_certificate import (
    VerificationResult,
    full_support_markers_float,
    full_support_term_arb,
    result_json,
    singleton_block_term_arb,
    singleton_exact_band_term_arb,
    validate_coverage,
)
from prime_field_biregular_ec_diagnostic import (
    candidate_parameters,
    singleton_constraint_balanced_ec_block_logbound,
)
from prime_field_ea_diagnostic import gv_distance


SCHEMA = "finite-field-biregular-ec-singleton-trace-v1"
LEGACY_SCHEMA = "prime-field-biregular-ec-singleton-trace-v1"
FIELD_ORDER = 2**127 - 1
DEGREE = 26
MEMORY = 4
K, N, REGION_LENGTH = candidate_parameters(2**21, DEGREE)
CUTOFF = math.floor(N * gv_distance(FIELD_ORDER, 0.5))

EXACT_MARKERS: dict[tuple[int, int], dict[str, tuple[str, str]]] = {
    (1, 8): {
        "structural": ("0.9991447773967677", "5.839835298989263e-22"),
        "field": ("0.9991692212810199", "2.5166797721728796e-21"),
    },
    (9, 16): {
        "structural": ("0.9967622332862647", "2.4826055705733867e-55"),
        "field": ("0.9975785283404316", "5.8774717541116e-39"),
    },
    (17, 32): {
        "structural": ("0.9985604162688718", "3.6970904119057333e-13"),
        "field": ("0.998592495268847", "9.185341959987983e-13"),
    },
    (33, 64): {
        "structural": ("0.9985527284526642", "9.420467108948146e-6"),
        "field": ("0.9985564845711785", "1.0295943023022584e-5"),
    },
    (65, 96): {
        "structural": ("0.9979801845055196", "6.855802703817471e-5"),
        "field": ("0.9979785655536838", "6.945969394666614e-5"),
    },
    (97, 128): {
        "structural": ("0.9972801872091132", "9.252994173875877e-5"),
        "field": ("0.997278476458764", "9.36589629959553e-5"),
    },
    (129, 160): {
        "structural": ("0.9963851659450441", "1.229659570575988e-4"),
        "field": ("0.9963794932070447", "1.2426853116162396e-4"),
    },
    (161, 192): {
        "structural": ("0.9954851021587344", "1.5392000341824723e-4"),
        "field": ("0.995488652002496", "1.5441616026343332e-4"),
    },
    (193, 256): {
        "structural": ("0.994", "0.00015"),
        "field": ("0.994", "0.00015"),
    },
    (257, 400): {
        "structural": ("0.992", "0.00025"),
        "field": ("0.992", "0.00025"),
    },
    (401, 800): {
        "structural": ("0.9865", "0.0005"),
        "field": ("0.9865", "0.0005"),
    },
}

BLOCKS = (
    (801, 1_600),
    (1_601, 3_200),
    (3_201, 6_400),
    (6_401, 12_800),
    (12_801, 25_600),
    (25_601, 51_200),
    (51_201, 102_400),
    (102_401, 204_800),
    (204_801, 409_600),
    (409_601, 819_200),
    (819_201, 917_499),
    (917_500, 983_035),
    (983_036, 1_015_803),
    (1_015_804, 1_032_187),
    (1_032_188, 1_040_379),
    (1_040_380, 1_044_475),
    (1_044_476, 1_046_523),
    (1_046_524, 1_047_547),
    (1_047_548, 1_048_059),
    (1_048_060, 1_048_315),
    (1_048_316, 1_048_443),
    (1_048_444, 1_048_507),
    (1_048_508, 1_048_539),
    (1_048_540, 1_048_555),
    (1_048_556, 1_048_563),
    (1_048_564, 1_048_566),
)


def generate_certificate(
    *, target_bits: int, precision_bits: int,
    field_order: int = FIELD_ORDER, cutoff: int | None = None,
) -> dict[str, Any]:
    if cutoff is None:
        cutoff = math.floor(N * gv_distance(field_order, 0.5))
    minimum_field_marker = math.nextafter(
        1.0 / (field_order - 1), math.inf
    )
    exact_bands = [
        {
            "lo": lo,
            "hi": hi,
            "markers": {
                name: [
                    values[0],
                    decimal_marker(max(float(values[1]), minimum_field_marker))
                    if name == "field" else values[1],
                ]
                for name, values in markers.items()
            },
        }
        for (lo, hi), markers in EXACT_MARKERS.items()
    ]
    blocks = []
    for lo, hi in BLOCKS:
        bound = singleton_constraint_balanced_ec_block_logbound(
            prime=field_order,
            k=K,
            n=N,
            cutoff=cutoff,
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
        print(f"optimized block [{lo},{hi}]", flush=True)
    return {
        "schema": SCHEMA,
        "parameters": {
            "field_order": str(field_order),
            "k": K,
            "n": N,
            "cutoff": cutoff,
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
                prime=field_order,
                n=N,
                cutoff=cutoff,
                memory=MEMORY,
                message_weight=K,
            ),
        },
    }


def verify_certificate(certificate: dict[str, Any]) -> VerificationResult:
    if certificate.get("schema") not in {SCHEMA, LEGACY_SCHEMA}:
        raise ValueError("unsupported certificate schema")
    validate_coverage(certificate)
    parameters = certificate["parameters"]
    field_order = int(parameters.get("field_order", parameters.get("prime")))
    k, n = int(parameters["k"]), int(parameters["n"])
    cutoff = int(parameters["cutoff"])
    region_count = int(parameters["region_count"])
    memory = int(parameters["memory"])
    target_bits = int(parameters["target_bits"])
    ctx.prec = int(certificate["verification"]["precision_bits"])

    exact_total = arb(0)
    largest_exact_band, largest_exact_term = (0, 0), arb(0)
    for item in certificate["exact_bands"]:
        lo, hi = int(item["lo"]), int(item["hi"])
        term = singleton_exact_band_term_arb(
            prime=field_order,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=region_count,
            memory=memory,
            lo=lo,
            hi=hi,
            markers=item["markers"],
        )
        exact_total += term
        if largest_exact_band == (0, 0) or term > largest_exact_term:
            largest_exact_band, largest_exact_term = (lo, hi), term
        print(f"verified exact band [{lo},{hi}]", flush=True)

    block_total = arb(0)
    largest_block, largest_block_term = (0, 0), arb(0)
    for item in certificate["blocks"]:
        lo, hi = int(item["lo"]), int(item["hi"])
        term = singleton_block_term_arb(
            prime=field_order,
            k=k,
            n=n,
            cutoff=cutoff,
            region_count=region_count,
            memory=memory,
            lo=lo,
            hi=hi,
            markers=item["markers"],
        )
        block_total += term
        if largest_block == (0, 0) or term > largest_block_term:
            largest_block, largest_block_term = (lo, hi), term
        print(f"verified block [{lo},{hi}]", flush=True)

    full = certificate["full_support"]
    full_term = full_support_term_arb(
        prime=field_order,
        n=n,
        cutoff=cutoff,
        memory=memory,
        message_weight=k,
        markers=full["markers"],
    )
    total = exact_total + block_total + full_term
    return VerificationResult(
        success=bool(total < arb(2) ** (-target_bits)),
        total_bound=total,
        security_bits=-total.log() / arb(2).log(),
        exact_bound=exact_total,
        block_bound=block_total,
        full_support_bound=full_term,
        largest_exact_band=largest_exact_band,
        largest_exact_term=largest_exact_term,
        largest_block=largest_block,
        largest_block_term=largest_block_term,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--output", type=Path, required=True)
    generate.add_argument("--target-bits", type=int, default=20)
    generate.add_argument("--precision-bits", type=int, default=256)
    generate.add_argument("--field-order", type=int, default=FIELD_ORDER)
    generate.add_argument("--cutoff", type=int)
    rebind = subparsers.add_parser(
        "rebind",
        help="reuse valid positive markers at another field order",
    )
    rebind.add_argument("certificate", type=Path)
    rebind.add_argument("--output", type=Path, required=True)
    rebind.add_argument("--field-order", type=int, required=True)
    rebind.add_argument(
        "--cutoff", type=int,
        help="output-weight cutoff (default: the field-order GV cutoff)",
    )
    verify = subparsers.add_parser("verify")
    verify.add_argument("certificate", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        certificate = generate_certificate(
            target_bits=args.target_bits,
            precision_bits=args.precision_bits,
            field_order=args.field_order,
            cutoff=args.cutoff,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(certificate, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.output}")
        return
    if args.command == "rebind":
        certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
        parameters = certificate["parameters"]
        certificate["schema"] = SCHEMA
        parameters["field_order"] = str(args.field_order)
        parameters.pop("prime", None)
        parameters["cutoff"] = args.cutoff
        if parameters["cutoff"] is None:
            parameters["cutoff"] = math.floor(
                int(parameters["n"]) * gv_distance(args.field_order, 0.5)
            )
        minimum_field_marker = math.nextafter(
            1.0 / (args.field_order - 1), math.inf
        )
        for item in certificate["exact_bands"]:
            markers = item["markers"]["field"]
            markers[1] = decimal_marker(
                max(float(markers[1]), minimum_field_marker)
            )
        for item in certificate["blocks"]:
            item.pop("diagnostic_log2", None)
            markers = item["markers"]["field"]
            markers[1] = decimal_marker(
                max(float(markers[1]), minimum_field_marker)
            )
        full_markers = certificate["full_support"]["markers"]["field"]
        full_markers[1] = decimal_marker(
            max(float(full_markers[1]), minimum_field_marker)
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
