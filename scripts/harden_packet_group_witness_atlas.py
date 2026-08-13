#!/usr/bin/env python3
"""Outward-harden every witness row in one packet-group atlas."""

from __future__ import annotations

import argparse
from pathlib import Path

import certify_packet_group_g4_anchor_mesh as g4
import certify_packet_group_triangle_ledger as verifier


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--group", type=int, required=True)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--checkpoint-dir", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error("--workers must lie in 1..8")
    if args.group == 4:
        g4.configure_generic_hardener()
    else:
        verifier.GROUP_BITS = args.group
        verifier.BLOCK_ATOMS = 64 // args.group
    witnesses = verifier.load_witnesses([args.atlas])
    names = [f"{args.atlas.name}:{index}" for index in range(len(verifier._rows_from_artifact(args.atlas)))]
    verifier.harden_selected_witnesses(
        names,
        witnesses,
        args.iterations,
        args.workers,
        args.checkpoint_dir,
        group_bits=args.group,
    )
    print(f"outward_hardened_witnesses={len(names)}")


if __name__ == "__main__":
    main()
