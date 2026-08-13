#!/usr/bin/env python3
"""Greedily close the 511 uniform-subset chamber vertices.

The ordered-chamber decomposition uses one uniform profile for every nonempty
subset of the nine packet-weight classes.  This diagnostic repeatedly selects
the currently worst uncovered subset vertex, tunes broad and sharp witnesses
at the nearest exact integer profile, appends them to the adaptive atlas, and
updates the matching fixed-witness cache.

The resulting rows are binary64 atlas-discovery witnesses, not outward
certificates.
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
from probe_packet8_profile_ordered_chambers import (
    add_full_bijection,
    uniform_subset,
)
from probe_packet8_profile_simplex_landscape import (
    ANCHORS,
    M,
    PROFILE_COUNT_LOG2,
    build_fixed_spectrum_witness,
    build_fixed_witness,
    load_anchor_file,
    normalization_log2,
)


def exact_subset_profile(subset: int) -> list[int]:
    active = [index for index in range(9) if subset >> index & 1]
    quotient, remainder = divmod(M, len(active))
    result = [0] * 9
    for position, packet_class in enumerate(active):
        result[packet_class] = quotient + (position < remainder)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--poles", default="0.1,0.15,0.2,0.3")
    parser.add_argument("--inactive-floor", type=float, default=0.2)
    parser.add_argument("--sharp-floor", type=float, default=1e-30)
    parser.add_argument("--outer-bounds", default="3,4,6,10,40")
    parser.add_argument("--anchor-target", type=float, default=-1000.0)
    parser.add_argument("--tuning-passes", type=int, default=1)
    parser.add_argument("--coordinate-iterations", type=int, default=6)
    parser.add_argument("--witness-iterations", type=int, default=32)
    args = parser.parse_args()
    args.poles = [float(value) for value in args.poles.split(",")]
    args.outer_bounds = [float(value) for value in args.outer_bounds.split(",")]
    if args.iterations <= 0:
        raise SystemExit("uniform vertex atlas: invalid iteration count")

    extras = list(load_anchor_file(args.atlas))
    anchors = ANCHORS + tuple(extras)
    cache_path = args.atlas.with_suffix(".witnesses.npz")
    witnesses = load_witness_cache(cache_path, anchors)
    if witnesses is None:
        raise SystemExit("uniform vertex atlas: matching witness cache missing")
    split_caps = split_cap_table()
    target = -40.0 - PROFILE_COUNT_LOG2
    points = np.vstack([uniform_subset(subset) for subset in range(1, 1 << 9)])
    normalizations = normalization_log2(points)

    def values_for(current_witnesses):
        rows = add_full_bijection(list(current_witnesses))
        constants = np.asarray([row.constant_log2 for row in rows])
        charges = np.vstack([row.linear_charge for row in rows])
        return np.min(
            constants[None, :] - points @ charges.T - normalizations[:, None],
            axis=1,
        )

    values = values_for(witnesses)
    print("packet-8 uniform-subset vertex atlas")
    print(
        f"atlas={args.atlas} initial_witnesses={len(witnesses)} "
        f"target={target:.12f} initial_covered={int(np.sum(values <= target))}/511",
        flush=True,
    )
    for iteration in range(args.iterations):
        uncovered = np.flatnonzero(values > target)
        if not len(uncovered):
            print("uniform_vertices_covered=PASS", flush=True)
            break
        selected_index = int(uncovered[np.argmax(values[uncovered])])
        subset = selected_index + 1
        profile = exact_subset_profile(subset)
        before = float(values[selected_index])

        broad = make_anchor(
            profile,
            f"uniform_{subset:03x}_broad",
            split_caps,
            args,
            args.inactive_floor,
        )
        sharp = make_anchor(
            profile,
            f"uniform_{subset:03x}_sharp",
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
        witnesses.extend(new_witnesses)
        extras.extend(new_anchors)
        write_anchors(args.atlas, extras)
        save_witness_cache(
            cache_path, ANCHORS + tuple(extras), witnesses
        )
        old_values = values
        values = values_for(witnesses)
        newly_covered = int(
            np.sum((old_values > target) & (values <= target))
        )
        print(
            f"iteration={iteration} subset={subset:03x} selected_before={before:.9f} "
            f"broad_self={broad[1]:.9f} sharp_self={sharp[1]:.9f} "
            f"newly_covered={newly_covered} "
            f"covered={int(np.sum(values <= target))}/511 "
            f"worst_after={float(np.max(values)):.9f}",
            flush=True,
        )
    print(f"final_covered={int(np.sum(values <= target))}/511")
    print(f"final_uncovered={int(np.sum(values > target))}")
    print("status=DIAGNOSTIC_BINARY64_UNIFORM_VERTEX_ATLAS")


if __name__ == "__main__":
    main()
