#!/usr/bin/env python3
"""Greedily close uniform Hasse-edge witness-interval gaps.

Every uniform-subset vertex must already be covered.  At each iteration this
diagnostic finds the widest first gap in the union of fixed-witness safe
intervals over all relevant Hasse edges, tunes broad and support-local
witnesses at the integralized midpoint, persists them, and repeats.

The interval geometry is continuous and sample-independent, but all tuning,
root finding, and comparisons remain binary64 rather than outward certified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from probe_packet8_adaptive_simplex_atlas import (
    load_witness_cache,
    make_anchor,
    save_witness_cache,
    split_cap_table,
    write_anchors,
)
from probe_packet8_profile_edge_interval_cover import edge_cover
from probe_packet8_profile_ordered_chambers import (
    add_full_bijection,
    uniform_subset,
)
from probe_packet8_profile_simplex_landscape import (
    ANCHORS,
    PROFILE_COUNT_LOG2,
    build_fixed_spectrum_witness,
    build_fixed_witness,
    integral_profile,
    load_anchor_file,
)


def all_edges(family: str) -> list[tuple[int, int]]:
    result = []
    for subset in range(2, 1 << 9):
        complement = ((1 << 9) - 1) ^ subset
        if family == "hasse":
            additions = [
                1 << packet_class
                for packet_class in range(9)
                if complement >> packet_class & 1
            ]
        else:
            additions = []
            addition = complement
            while addition:
                additions.append(addition)
                addition = (addition - 1) & complement
        result.extend((subset, subset | addition) for addition in additions)
    return result


def edge_gaps(witnesses, target: float, edges: list[tuple[int, int]]):
    rows = add_full_bijection(list(witnesses))
    constants = np.asarray([row.constant_log2 for row in rows])
    charges = np.vstack([row.linear_charge for row in rows])
    gaps = []
    for subset, larger in edges:
        left = uniform_subset(subset)
        is_covered, gap_left, gap_right, _leader = edge_cover(
            left,
            uniform_subset(larger),
            constants,
            charges,
            target,
        )
        if is_covered:
            continue
        point = left + 0.5 * (gap_left + gap_right) * (
            uniform_subset(larger) - left
        )
        gaps.append(
            (
                gap_right - gap_left,
                subset,
                larger,
                gap_left,
                gap_right,
                integral_profile(point).astype(int).tolist(),
            )
        )
    gaps.sort(reverse=True)
    return gaps


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument(
        "--edge-family", choices=("hasse", "comparable"), default="hasse"
    )
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--residual-cache", type=Path)
    parser.add_argument("--poles", default="0.1,0.15,0.2,0.3")
    parser.add_argument("--inactive-floor", type=float, default=0.2)
    parser.add_argument("--sharp-floor", type=float, default=1e-30)
    parser.add_argument("--outer-bounds", default="3,4,6,10,40")
    parser.add_argument("--anchor-target", type=float, default=-1000.0)
    parser.add_argument("--tuning-passes", type=int, default=1)
    parser.add_argument("--coordinate-iterations", type=int, default=6)
    parser.add_argument("--witness-iterations", type=int, default=32)
    args = parser.parse_args()
    if args.iterations <= 0 or args.batch_size <= 0:
        raise SystemExit("edge interval atlas: invalid iteration or batch size")
    args.poles = [float(value) for value in args.poles.split(",")]
    args.outer_bounds = [float(value) for value in args.outer_bounds.split(",")]

    extras = list(load_anchor_file(args.atlas))
    anchors = ANCHORS + tuple(extras)
    cache_path = args.atlas.with_suffix(".witnesses.npz")
    witnesses = load_witness_cache(cache_path, anchors)
    if witnesses is None:
        raise SystemExit("edge interval atlas: matching witness cache missing")
    split_caps = split_cap_table()
    target = -40.0 - PROFILE_COUNT_LOG2

    print(f"packet-8 greedy {args.edge_family}-edge interval atlas", flush=True)
    print(
        f"atlas={args.atlas} initial_witnesses={len(witnesses)} "
        f"target={target:.12f}",
        flush=True,
    )
    source_path = Path(__file__).with_name(
        "probe_packet8_profile_edge_interval_cover.py"
    )
    source_digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
    witness_names = [row.name for row in witnesses]
    edges = all_edges(args.edge_family)
    total = len(edges)
    active_edges = edges
    if args.residual_cache is not None and args.residual_cache.exists():
        with args.residual_cache.open(encoding="utf-8") as handle:
            cached = json.load(handle)
        cached_names = list(cached.get("witness_names", ()))
        if cached.get("edge_family") != args.edge_family:
            raise SystemExit("edge interval atlas: residual family mismatch")
        if cached.get("edge_cover_source_sha256") != source_digest:
            raise SystemExit("edge interval atlas: residual source mismatch")
        if witness_names[: len(cached_names)] != cached_names:
            raise SystemExit("edge interval atlas: residual witness prefix mismatch")
        active_edges = [tuple(map(int, row)) for row in cached["edges"]]
        print(
            f"residual_cache=HIT active_edges={len(active_edges)}/{total} "
            f"cached_witnesses={len(cached_names)}",
            flush=True,
        )
    for iteration in range(args.iterations):
        gaps = edge_gaps(witnesses, target, active_edges)
        covered = total - len(gaps)
        if not gaps:
            print(f"all_edges_covered=PASS total={total}", flush=True)
            break
        batch = gaps[: args.batch_size]
        new_anchor_rows = []
        for batch_index, (
            width,
            subset,
            larger,
            gap_left,
            gap_right,
            profile,
        ) in enumerate(batch):
            label = f"edge_{subset:03x}_{larger:03x}_{len(extras):03d}"
            broad = make_anchor(
                profile,
                label + "_broad",
                split_caps,
                args,
                args.inactive_floor,
            )
            sharp = make_anchor(
                profile,
                label + "_sharp",
                split_caps,
                args,
                args.sharp_floor,
            )
            new_anchor_rows.extend((broad[0], sharp[0]))
            print(
                f"iteration={iteration} batch={batch_index} "
                f"covered_before={covered}/{total} "
                f"edge={subset:03x}->{larger:03x} "
                f"gap={gap_left:.12f},{gap_right:.12f} width={width:.12f} "
                f"broad_self={broad[1]:.9f} sharp_self={sharp[1]:.9f} "
                f"profile=" + ",".join(map(str, profile)),
                flush=True,
            )
        for anchor in new_anchor_rows:
            witnesses.extend(
                (
                    build_fixed_witness(anchor, split_caps),
                    build_fixed_spectrum_witness(anchor, split_caps),
                )
            )
        extras.extend(new_anchor_rows)
        write_anchors(args.atlas, extras)
        save_witness_cache(cache_path, ANCHORS + tuple(extras), witnesses)
        active_edges = [(row[1], row[2]) for row in gaps]

    gaps = edge_gaps(witnesses, target, active_edges)
    covered = total - len(gaps)
    if args.residual_cache is not None:
        args.residual_cache.parent.mkdir(parents=True, exist_ok=True)
        with args.residual_cache.open("w", encoding="utf-8") as handle:
            json.dump(
                {
                    "edge_family": args.edge_family,
                    "edge_cover_source_sha256": source_digest,
                    "witness_names": [row.name for row in witnesses],
                    "edges": [[row[1], row[2]] for row in gaps],
                },
                handle,
                indent=2,
            )
            handle.write("\n")
    print(f"final_covered_edges={covered}/{total}")
    print(f"final_uncovered_edges={total-covered}")
    if gaps:
        print(f"final_widest_gap_width={gaps[0][0]:.12f}")
    print("status=DIAGNOSTIC_BINARY64_EDGE_INTERVAL_ATLAS")


if __name__ == "__main__":
    main()
