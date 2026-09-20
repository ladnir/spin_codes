#!/usr/bin/env python3
"""Audit containment and class accounting in the XBCH outward certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from flint import arb, ctx


COMPONENTS = (
    "exactly_one_tail_no_body",
    "exactly_two_tail_no_body",
    "three_or_more_tail_no_body",
    "body_no_tail",
    "body_and_tail",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outward", type=Path, required=True)
    parser.add_argument("--nearest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ctx.prec = 256
    outward = json.loads(args.outward.read_text(encoding="utf-8"))
    nearest = json.loads(args.nearest.read_text(encoding="utf-8"))
    if set(outward["components"]) != set(COMPONENTS):
        raise AssertionError("outward component partition changed")
    if set(nearest["components"]) != set(COMPONENTS):
        raise AssertionError("nearest component partition changed")

    rows = []
    for name in COMPONENTS:
        lower = float(outward["components"][name]["margin_bits_lower"])
        diagnostic = float(nearest["components"][name]["margin_bits"])
        loss = diagnostic - lower
        if loss < -1e-12:
            raise AssertionError(f"{name}: outward bound is not conservative")
        if loss > 1e-3:
            raise AssertionError(f"{name}: unexpected outward loss {loss}")
        rows.append(
            {
                "component": name,
                "outward_margin_bits_lower": lower,
                "nearest_margin_bits": diagnostic,
                "outward_loss_bits": loss,
            }
        )

    logs = [
        arb(str(outward["components"][name]["log2_upper"]))
        for name in COMPONENTS
    ]
    pivot = max(float(value) for value in logs)
    recomputed = arb(pivot) + sum(
        arb(2) ** (value - pivot) for value in logs
    ).log() / arb(2).log()
    recorded = float(outward["aggregate_log2_upper"])
    if recorded + 1e-12 < float(recomputed.upper()):
        raise AssertionError("recorded aggregate is below its component sum")
    margin = float(outward["aggregate_margin_bits_lower"])
    if margin <= 40.0:
        raise AssertionError("outward certificate does not clear 40 bits")
    if int(outward["parameters"]["arb_precision_bits"]) < 256:
        raise AssertionError("retained receipt was not generated at 256 bits")

    payload = {
        "schema": "riffle-xbch-outward-certificate-audit-v1",
        "status": "PASS",
        "aggregate_margin_bits_lower": margin,
        "reserve_over_40_bits": margin - 40.0,
        "component_containment": rows,
        "aggregate_recomputed_log2_upper": float(recomputed.upper()),
        "inputs": {
            "outward": str(args.outward),
            "outward_sha256": sha256(args.outward),
            "nearest": str(args.nearest),
            "nearest_sha256": sha256(args.nearest),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print("status,PASS")
    print(f"aggregate_margin_bits_lower,{margin:.12f}")
    print(f"reserve_over_40_bits,{margin - 40.0:.12f}")
    print(f"wrote,{args.output}")


if __name__ == "__main__":
    main()
