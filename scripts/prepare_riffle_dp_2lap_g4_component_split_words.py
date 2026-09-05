#!/usr/bin/env python3
"""Materialize component-split generator words for the exact pair engine."""

from __future__ import annotations

import argparse
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import accumulate, binary_rank, build_apply
from audit_riffle_dp_2lap_g4_component_endpoint_dimension30 import (
    observations,
    support_dimension,
)
from probe_riffle_dp_g4_component_mixing import transpose_columns
from solve_riffle_dp_2lap_g4_component_split_sat import (
    balanced_split,
    component_basis,
    generator_words,
)


def limbs(word: int) -> tuple[int, ...]:
    return tuple((word >> (64 * index)) & ((1 << 64) - 1) for index in range(6))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support", type=lambda value: int(value, 0), required=True)
    parser.add_argument("--packet-value", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    transpose_step = build_apply(
        transpose_columns(tuple(step(1 << bit) for bit in range(64)))
    )
    left_mask, right_mask = balanced_split(args.support)
    left_basis = component_basis(transpose_step, left_mask)
    right_basis = component_basis(transpose_step, right_mask)
    if binary_rank(left_basis + right_basis) != support_dimension(args.support):
        raise RuntimeError("component-split words: basis mismatch")
    states = observations(step, args.packet_value, 24)
    left_words = generator_words(states, left_basis)
    right_words = generator_words(states, right_basis)

    lines = [
        "RIFFLE_DP_2LAP_G4_COMPONENT_SPLIT_WORDS_V1",
        (
            f"{args.support:x} {args.packet_value} {left_mask:x} {right_mask:x} "
            f"{len(left_words)} {len(right_words)}"
        ),
    ]
    for side, words in (("L", left_words), ("R", right_words)):
        for word in words:
            lines.append(side + " " + " ".join(f"{limb:016x}" for limb in limbs(word)))
    args.output.write_text("\n".join(lines) + "\n")
    print(
        f"support={hex(args.support)} split={hex(left_mask)}+{hex(right_mask)} "
        f"dimensions={len(left_words)}+{len(right_words)}"
    )
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
