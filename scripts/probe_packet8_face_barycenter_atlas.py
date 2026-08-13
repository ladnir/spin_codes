#!/usr/bin/env python3
"""Greedily close all nested-prefix triangle barycenters.

This is a higher-dimensional diagnostic after vertex and edge exploration.  It
enumerates every strict nested triple of uniform-subset vertices, evaluates
all 198580 triangle barycenters, then repeatedly tunes broad and support-local
witnesses at the worst uncovered barycenter.  New fixed witnesses are applied
to every barycenter and persisted in the shared atlas/cache.

Barycenter closure is sample-independent for those points but is not a full
triangle or chamber certificate.  Arithmetic remains binary64.
"""

from __future__ import annotations

import argparse
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
from probe_packet8_nested_face_landscape import nested_triples
from probe_packet8_profile_edge_interval_cover import orbit_value
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
    normalization_log2,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--poles", default="0.1,0.15,0.2,0.3")
    parser.add_argument("--inactive-floor", type=float, default=0.2)
    parser.add_argument("--sharp-floor", type=float, default=1e-30)
    parser.add_argument("--outer-bounds", default="3,4,6,10,40")
    parser.add_argument("--anchor-target", type=float, default=-50000.0)
    parser.add_argument("--tuning-passes", type=int, default=1)
    parser.add_argument("--coordinate-iterations", type=int, default=6)
    parser.add_argument("--witness-iterations", type=int, default=32)
    args = parser.parse_args()
    if args.iterations <= 0:
        raise SystemExit("face barycenter atlas: invalid iteration count")
    args.poles = [float(value) for value in args.poles.split(",")]
    args.outer_bounds = [float(value) for value in args.outer_bounds.split(",")]

    extras = list(load_anchor_file(args.atlas))
    anchors = ANCHORS + tuple(extras)
    cache_path = args.atlas.with_suffix(".witnesses.npz")
    witnesses = load_witness_cache(cache_path, anchors)
    if witnesses is None:
        raise SystemExit("face barycenter atlas: matching witness cache missing")
    split_caps = split_cap_table()
    target = -40.0 - PROFILE_COUNT_LOG2

    triples = nested_triples()
    points = {subset: uniform_subset(subset) for subset in range(2, 1 << 9)}
    profiles = np.vstack(
        [
            (points[left] + points[middle] + points[right]) / 3.0
            for left, middle, right in triples
        ]
    )
    normalizations = normalization_log2(profiles)

    rows = add_full_bijection(list(witnesses))
    values = np.full(len(profiles), np.inf)
    for start in range(0, len(rows), 128):
        batch = rows[start : start + 128]
        constants = np.asarray([row.constant_log2 for row in batch])
        charges = np.vstack([row.linear_charge for row in batch])
        values = np.minimum(
            values,
            np.min(
                constants[None, :]
                - profiles @ charges.T
                - normalizations[:, None],
                axis=1,
            ),
        )
    orbit = np.asarray(
        [orbit_value(0.0, point, np.zeros(9), 0.0) for point in profiles]
    )
    values = np.minimum(values, orbit)

    print("packet-8 greedy nested-face barycenter atlas", flush=True)
    print(
        f"atlas={args.atlas} barycenters={len(profiles)} "
        f"initial_witnesses={len(witnesses)} target={target:.12f} "
        f"initial_covered={int(np.sum(values <= target))}/{len(values)}",
        flush=True,
    )
    for iteration in range(args.iterations):
        uncovered = np.flatnonzero(values > target)
        if not len(uncovered):
            print("all_barycenters_covered=PASS", flush=True)
            break
        selected = int(uncovered[np.argmax(values[uncovered])])
        left, middle, right = triples[selected]
        profile = integral_profile(profiles[selected]).astype(int).tolist()
        before = float(values[selected])
        label = f"face_{left:03x}_{middle:03x}_{right:03x}_{len(extras):03d}"
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
        new_anchors = (broad[0], sharp[0])
        new_witnesses = []
        for anchor in new_anchors:
            new_witnesses.extend(
                (
                    build_fixed_witness(anchor, split_caps),
                    build_fixed_spectrum_witness(anchor, split_caps),
                )
            )
        old_values = values
        constants = np.asarray([row.constant_log2 for row in new_witnesses])
        charges = np.vstack([row.linear_charge for row in new_witnesses])
        candidates = (
            constants[None, :]
            - profiles @ charges.T
            - normalizations[:, None]
        )
        values = np.minimum(values, np.min(candidates, axis=1))
        witnesses.extend(new_witnesses)
        extras.extend(new_anchors)
        write_anchors(args.atlas, extras)
        save_witness_cache(cache_path, ANCHORS + tuple(extras), witnesses)
        print(
            f"iteration={iteration} face={left:03x}<{middle:03x}<{right:03x} "
            f"selected_before={before:.9f} broad_self={broad[1]:.9f} "
            f"sharp_self={sharp[1]:.9f} "
            f"newly_covered={int(np.sum((old_values > target) & (values <= target)))} "
            f"covered={int(np.sum(values <= target))}/{len(values)} "
            f"worst_after={float(np.max(values)):.9f} profile="
            + ",".join(map(str, profile)),
            flush=True,
        )
    print(f"final_covered={int(np.sum(values <= target))}/{len(values)}")
    print(f"final_uncovered={int(np.sum(values > target))}")
    print("status=DIAGNOSTIC_BINARY64_FACE_BARYCENTER_ATLAS")


if __name__ == "__main__":
    main()
