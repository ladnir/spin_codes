#!/usr/bin/env python3
"""Tune shared-drive anchors at uniform profiles of residual support masks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from probe_packet8_drive_residual_anchor import retune_shared_report
from probe_packet8_profile_simplex_landscape import M


def uniform_profile(mask: int) -> list[int]:
    active = [weight for weight in range(9) if mask >> weight & 1]
    if not active:
        raise ValueError("empty support mask")
    quotient, remainder = divmod(M, len(active))
    profile = [0] * 9
    for rank, weight in enumerate(active):
        profile[weight] = quotient + (rank < remainder)
    return profile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--count", type=int, default=10)
    parser.add_argument("--skip-mask", action="append", default=[])
    parser.add_argument("--coordinate-iterations", type=int, default=4)
    parser.add_argument("--poles", default="0.1,0.2,0.3")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    census = json.loads(args.census.read_text(encoding="utf-8"))
    skipped = {int(value, 16) for value in args.skip_mask}
    masks = [
        int(row["mask"], 16)
        for row in census["uncovered_support_masks"]
        if int(row["mask"], 16) not in skipped
    ][: args.count]
    poles = [float(value) for value in args.poles.split(",")]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for rank, mask in enumerate(masks):
        profile = uniform_profile(mask)
        report = retune_shared_report(
            profile,
            f"support_{mask:03x}_uniform",
            poles=poles,
            coordinate_iterations=args.coordinate_iterations,
        )
        path = args.output_dir / f"support_{mask:03x}_uniform.json"
        path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(
            f"rank={rank} mask={mask:03x} "
            f"combined={report['shared_drive_combined_log2']:.9f} "
            f"improvement={report['shared_drive_improvement_log2']:.9f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
