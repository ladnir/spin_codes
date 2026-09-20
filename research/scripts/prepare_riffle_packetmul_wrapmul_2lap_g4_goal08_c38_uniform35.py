#!/usr/bin/env python3
"""Pack 35 disjoint information sets for each requested C38 zero-anchor shortening."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from analyze_riffle_packetmul_wrapmul_2lap_g4_lifted_windows import (
    build_apply,
    observation_rows,
)
from analyze_systematic_group_kernel import systematic_state_columns
from prepare_riffle_packetmul_wrapmul_2lap_g4_goal08_c38_information_sets import (
    DIMENSION,
    LENGTH,
    NODES,
    pack,
)
from probe_riffle_packetmul_wrapmul_2lap_g4_goal05_zero_anchor import (
    MASK64,
    delete_node,
    kernel_basis,
    xor_selected,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
PACKER_SOURCE = (
    ROOT / "scripts" / "prepare_riffle_packetmul_wrapmul_2lap_g4_goal08_c38_information_sets.py"
)
SET_COUNT = 35


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor-start", type=int, required=True)
    parser.add_argument("--anchor-end", type=int, required=True)
    args = parser.parse_args()
    if not 0 <= args.anchor_start <= args.anchor_end < NODES:
        raise SystemExit("C38 uniform packing: require 0 <= start <= end < 38")

    parity = build_apply(systematic_state_columns())
    full_rows = observation_rows(NODES, parity)
    for anchor in range(args.anchor_start, args.anchor_end + 1):
        exact_path = (
            CANDIDATE / "receipts" / f"goal08_c38_anchor_{anchor:02d}_information_sets.json"
        )
        start = time.perf_counter()
        if exact_path.exists():
            exact = json.loads(exact_path.read_text())
            if exact["maximum_disjoint_information_sets"] < SET_COUNT:
                raise RuntimeError(f"C38 uniform packing: anchor {anchor} exact result is too small")
            sets = exact["information_sets"][:SET_COUNT]
            origin = {
                "method": "SELECTED_FROM_EXACT_MAXIMUM_PACKING",
                "exact_receipt_sha256": digest(exact_path),
            }
        else:
            basis = kernel_basis(
                tuple((row >> (64 * anchor)) & MASK64 for row in full_rows),
                64,
            )
            if len(basis) != DIMENSION:
                raise RuntimeError(f"C38 uniform packing: anchor {anchor} dimension changed")
            generators = tuple(
                delete_node(xor_selected(full_rows, state), anchor)
                for state in basis
            )
            columns = [
                sum(
                    ((generator >> coordinate) & 1) << bit
                    for bit, generator in enumerate(generators)
                )
                for coordinate in range(LENGTH)
            ]
            result = pack(SET_COUNT, columns)
            if not result["success"]:
                raise RuntimeError(
                    f"C38 uniform packing: anchor {anchor} packs only "
                    f"{result['maximum_union_size']}/{SET_COUNT * DIMENSION} elements"
                )
            sets = result["sets"]
            origin = {
                "method": "EXACT_LINEAR_MATROID_UNION_FEASIBILITY",
                "maximum_union_size": result["maximum_union_size"],
            }
        elapsed = time.perf_counter() - start

        flat = [coordinate for row in sets for coordinate in row]
        if len(sets) != SET_COUNT or any(len(row) != DIMENSION for row in sets):
            raise RuntimeError(f"C38 uniform packing: anchor {anchor} has incomplete sets")
        if len(flat) != len(set(flat)):
            raise RuntimeError(f"C38 uniform packing: anchor {anchor} sets overlap")

        output = (
            CANDIDATE / "receipts" / f"goal08_c38_anchor_{anchor:02d}_uniform35.json"
        )
        payload = {
            "schema": "riffle-packetmul-wrapmul-2lap-g4-goal08-c38-uniform35-v1",
            "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
            "evidence_label": "EXACT_DISJOINT_INFORMATION_SET_FEASIBILITY",
            "source_sha256": digest(Path(__file__).resolve()),
            "packer_source_sha256": digest(PACKER_SOURCE),
            "anchor_node": anchor,
            "window_nodes": NODES,
            "shortened_code": {"length": LENGTH, "dimension": DIMENSION},
            "information_set_count": SET_COUNT,
            "information_sets": sets,
            "origin": origin,
            "elapsed_seconds": elapsed,
            "result": "THIRTY_FIVE_DISJOINT_INFORMATION_SETS",
            "scope_limitation": (
                "This receipt proves feasibility of 35 disjoint information sets. "
                "It does not claim that 35 is the maximum for this anchor."
            ),
        }
        output.write_text(json.dumps(payload, indent=2) + "\n")
        print(
            f"anchor={anchor} sets={SET_COUNT} method={origin['method']} "
            f"seconds={elapsed:.6f}"
        )
    print("status=C38_UNIFORM35_BATCH_COMPLETE")


if __name__ == "__main__":
    main()
