#!/usr/bin/env python3
"""Interval certificate for the ultra-late full-split placement cap.

For the ultra-late branch the first active inner block lies in the final
``late_blocks`` blocks.  A placement-only union bound gives

    A_h * binom(b*late_blocks, h) / binom(N, h).

For a fixed outer Cauchy pole z,

    A_h <= (W_loc(z)^outer_blocks - 1) z^{-h}.

The resulting exponent is concave in h, since the adjacent ratio is
``((m-h)/(N-h))/z`` with ``m=b*late_blocks``.  Therefore the interval maximum
is attained at an endpoint or at the unique critical point.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

from block_outer_upgrade_probe import load_local_spectrum, local_spectrum_log2, log2_sub_one, z_grid
from dense_largek_eval import log2_binom, log2add


@dataclass(frozen=True)
class Interval:
    start: int
    stop: int
    z: float | None = None

    @property
    def label(self) -> str:
        return f"{self.start}--{self.stop}"


def parse_intervals(text: str) -> list[Interval]:
    out: list[Interval] = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        fields = part.split(":")
        if len(fields) not in (2, 3):
            raise ValueError(f"bad interval {part!r}; expected h0:h1[:z]")
        start, stop = int(fields[0]), int(fields[1])
        if start > stop:
            raise ValueError(f"bad interval {part!r}: start > stop")
        z = None if len(fields) == 2 else float(fields[2])
        out.append(Interval(start, stop, z))
    return out


def fixed_z_outer_log(
    *,
    blocks: int,
    spectrum: list[tuple[int, int]],
    z: float,
) -> float:
    return log2_sub_one(blocks * local_spectrum_log2(spectrum, z))


def row_exponent(*, N: int, m: int, h: int, z: float, outer_log: float) -> float:
    if h < 0 or h > m or h > N:
        return float("-inf")
    return outer_log - h * math.log2(z) + log2_binom(m, h) - log2_binom(N, h)


def interval_bound(
    *,
    interval: Interval,
    N: int,
    m: int,
    outer_logs: list[tuple[float, float]],
) -> tuple[float, dict[str, float | int]]:
    lo = interval.start
    hi = min(interval.stop, m, N)
    if lo > hi:
        return float("-inf"), {"peak_h": -1, "z": float("nan"), "row_log2": float("-inf")}

    best = float("inf")
    best_z = float("nan")
    best_h = -1
    best_row = float("-inf")
    for z, outer_log in outer_logs:
        if outer_log == float("-inf"):
            continue
        critical = (m - z * N) / (1.0 - z)
        candidates = {lo, hi, math.floor(critical), math.ceil(critical)}
        worst = float("-inf")
        worst_h = -1
        for h in candidates:
            if lo <= h <= hi:
                value = row_exponent(N=N, m=m, h=h, z=z, outer_log=outer_log)
                if value > worst:
                    worst = value
                    worst_h = h
        value = math.log2(hi - lo + 1) + worst
        if value < best:
            best = value
            best_z = z
            best_h = worst_h
            best_row = worst
    return best, {"peak_h": best_h, "z": best_z, "row_log2": best_row}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-spectrum-csv", type=Path, default=Path(__file__).with_name("rm512_256_spectrum.csv"))
    parser.add_argument("--outer-blocks", type=int, default=4096)
    parser.add_argument("--N", type=int, default=2**21)
    parser.add_argument("--block-bits", type=int, default=64)
    parser.add_argument("--late-blocks", type=int, default=5949)
    parser.add_argument(
        "--intervals",
        required=True,
        help="Semicolon-separated h0:h1 or h0:h1:z rows. If z is omitted, the row optimizes over --z-grid.",
    )
    parser.add_argument("--z-min", type=float, default=1e-8)
    parser.add_argument("--z-max", type=float, default=0.999)
    parser.add_argument("--z-count", type=int, default=220)
    args = parser.parse_args()

    intervals = parse_intervals(args.intervals)
    spectrum = load_local_spectrum(str(args.local_spectrum_csv))
    default_outer_logs = [
        (z, fixed_z_outer_log(blocks=args.outer_blocks, spectrum=spectrum, z=z))
        for z in z_grid(args.z_min, args.z_max, args.z_count)
    ]
    m = args.block_bits * args.late_blocks

    total = float("-inf")
    print("Full-split ultra-late placement interval certificate")
    print(f"N,{args.N}")
    print(f"m_late,{m}")
    print(f"late_blocks,{args.late_blocks}")
    print("h_min,h_max,z,total_log2,peak_h,row_log2")
    for interval in intervals:
        outer_logs = default_outer_logs
        if interval.z is not None:
            outer_logs = [
                (
                    interval.z,
                    fixed_z_outer_log(blocks=args.outer_blocks, spectrum=spectrum, z=interval.z),
                )
            ]
        value, meta = interval_bound(interval=interval, N=args.N, m=m, outer_logs=outer_logs)
        total = log2add(total, value)
        print(
            f"{interval.start},{interval.stop},{float(meta['z']):.17g},"
            f"{value:.6f},{int(meta['peak_h'])},{float(meta['row_log2']):.6f}"
        )
    print(f"total_log2,{total:.6f}")
    print(f"margin_bits,{-total:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
