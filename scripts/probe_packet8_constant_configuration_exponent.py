#!/usr/bin/env python3
"""Configuration-model exponent for the constant-packet two-band code.

There are ``n0=42m/8`` band-zero packet variables and ``n1=43m/8`` band-one
packet variables, each repeated to eight block sockets.  Every one of the
``m`` block nodes enforces membership in the exact EBCH projection code with
split enumerator ``W(x,y)`` from ``analyze_ebch85_split_spectrum.py``.

For packet-variable densities ``rho0,rho1``, the ordinary two-edge-type
configuration ensemble has first-moment exponent per block

    -7(42/8 H(rho0) + 43/8 H(rho1))
      + inf_{x,y>0} log2 W(x,y)
          - 42 rho0 log2 x - 43 rho1 log2 y.

This probe evaluates that exponent on diagonal and two-dimensional grids.  It
is a diagnostic for the structured Riffle lift: the actual per-column lane
partitions and affine tile incidence are not yet replaced by the uniform
configuration ensemble in a proof.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


N0 = 42
N1 = 43
VARIABLES0_PER_BLOCK = N0 / 8
VARIABLES1_PER_BLOCK = N1 / 8


def load_spectrum(path: Path):
    cells = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cells.append(
                (
                    int(row["band0_weight"]),
                    int(row["band1_weight"]),
                    int(row["count"]),
                )
            )
    if sum(count for _a, _b, count in cells) != 1 << 64:
        raise SystemExit("constant configuration: local spectrum mass mismatch")
    weights0 = np.array([a for a, _b, _count in cells], dtype=np.float64)
    weights1 = np.array([b for _a, b, _count in cells], dtype=np.float64)
    logs = np.array([math.log(count) for _a, _b, count in cells])
    return weights0, weights1, logs


def entropy(value: float) -> float:
    if value <= 0.0 or value >= 1.0:
        return 0.0
    return -value * math.log2(value) - (1.0 - value) * math.log2(1.0 - value)


class Exponent:
    def __init__(self, spectrum_path: Path):
        self.weights0, self.weights1, self.logs = load_spectrum(spectrum_path)

    def local_entropy(self, rho0: float, rho1: float):
        if rho0 == 0.0 and rho1 == 0.0:
            return 0.0, np.array([-40.0, -40.0])
        if rho0 == 1.0 and rho1 == 1.0:
            return 0.0, np.array([40.0, 40.0])

        def objective(point):
            log_w = self.logs + self.weights0 * point[0] + self.weights1 * point[1]
            return (
                logsumexp(log_w)
                - N0 * rho0 * point[0]
                - N1 * rho1 * point[1]
            )

        def logit(value: float) -> float:
            clipped = min(max(value, 1e-8), 1.0 - 1e-8)
            return math.log(clipped / (1.0 - clipped))

        starts = (
            np.array([logit(rho0), logit(rho1)]),
            np.zeros(2),
        )
        candidates = []
        for start in starts:
            result = minimize(
                objective,
                start,
                method="L-BFGS-B",
                bounds=[(-40.0, 40.0), (-40.0, 40.0)],
                options={"maxiter": 1000, "ftol": 1e-14, "gtol": 1e-10},
            )
            candidates.append((float(result.fun), result.x))
        value, point = min(candidates, key=lambda row: row[0])
        return value / math.log(2.0), point

    def ensemble(self, rho0: float, rho1: float):
        local, point = self.local_entropy(rho0, rho1)
        exponent = local - 7.0 * (
            VARIABLES0_PER_BLOCK * entropy(rho0)
            + VARIABLES1_PER_BLOCK * entropy(rho1)
        )
        return exponent, local, point


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spectrum",
        type=Path,
        default=Path(__file__).resolve().parent.parent
        / "out"
        / "ebch85_band01_split_spectrum.csv",
    )
    parser.add_argument("--diagonal-step", type=float, default=0.01)
    parser.add_argument("--grid-step", type=float, default=0.05)
    args = parser.parse_args()
    if not 0.0 < args.diagonal_step <= 0.25 or not 0.0 < args.grid_step <= 0.25:
        raise SystemExit("constant configuration: invalid grid step")
    evaluator = Exponent(args.spectrum)

    diagonal = []
    steps = round(1.0 / args.diagonal_step)
    for index in range(steps + 1):
        rho = index / steps
        exponent, local, point = evaluator.ensemble(rho, rho)
        diagonal.append((exponent, rho, local, point))
    interior_diagonal = [row for row in diagonal if 0.0 < row[1] < 1.0]
    worst_diagonal = max(interior_diagonal)

    grid = []
    grid_steps = round(1.0 / args.grid_step)
    for first in range(grid_steps + 1):
        rho0 = first / grid_steps
        for second in range(grid_steps + 1):
            rho1 = second / grid_steps
            if (rho0, rho1) in ((0.0, 0.0), (1.0, 1.0)):
                continue
            exponent, local, point = evaluator.ensemble(rho0, rho1)
            grid.append((exponent, rho0, rho1, local, point))
    worst_grid = max(grid)

    print("constant-packet two-band configuration exponent probe")
    print(
        f"diagonal_worst_exponent_per_block={worst_diagonal[0]:.12f} "
        f"rho={worst_diagonal[1]:.6f}"
    )
    print(
        f"grid_worst_exponent_per_block={worst_grid[0]:.12f} "
        f"rho0={worst_grid[1]:.6f} rho1={worst_grid[2]:.6f}"
    )
    for rho in (0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99):
        exponent, local, point = evaluator.ensemble(rho, rho)
        print(
            f"diagonal_rho={rho:.6f} exponent_per_block={exponent:.12f} "
            f"local_entropy={local:.12f} log_x={point[0]:.8f} log_y={point[1]:.8f}"
        )
    print("model=UNIFORM_TWO_EDGE_TYPE_CONFIGURATION_ENSEMBLE")
    print("status=DIAGNOSTIC_BINARY64_CONFIGURATION_EXPONENT")


if __name__ == "__main__":
    main()
