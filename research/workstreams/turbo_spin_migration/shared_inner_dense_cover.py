"""Binary64 vertex cover for the scalable route with Turbo's s22 inner.

The exponent is convex in (alpha, alpha*x) for a fixed witness on each
affine outer segment, so four negative vertices cover a box analytically.
Floating-point evaluations and the discovered witnesses are NOT a proof.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
import math
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from research.workstreams.turbo_spin_migration.shared_inner_dense_screen import (
    OUTER, entropy_kl, optimize, transfer_rate,
)


def segments():
    payload = json.loads(OUTER.read_text(encoding="utf-8"))
    assert payload["status"] == "proved"
    left = []
    for record in payload["left_segments"]:
        lo, hi = map(lambda value: float(Fraction(value)), record["weight"])
        slope = float(Fraction(record["slope"]))
        intercept = float(Fraction(record["intercept"]))
        left.append((lo, hi, slope, intercept))
    center_lo = left[-1][1]
    center = [(center_lo, 1 - center_lo, 0.,
               float(Fraction(payload["central_constant_upper"])))]
    right = [(1 - hi, 1 - lo, -slope, slope + intercept)
             for lo, hi, slope, intercept in reversed(left)]
    return left + center + right


@dataclass(frozen=True)
class Box:
    alo: float
    ahi: float
    xlo: float
    xhi: float
    slope: float
    intercept: float
    depth: int = 0

    def vertices(self):
        return [(a, x) for a in (self.alo, self.ahi)
                for x in (self.xlo, self.xhi)]

    def split(self):
        if (self.ahi - self.alo) / (self.ahi + self.alo) > 4 * (self.xhi - self.xlo):
            middle = (self.alo + self.ahi) / 2
            return [Box(self.alo, middle, self.xlo, self.xhi,
                        self.slope, self.intercept, self.depth + 1),
                    Box(middle, self.ahi, self.xlo, self.xhi,
                        self.slope, self.intercept, self.depth + 1)]
        middle = (self.xlo + self.xhi) / 2
        return [Box(self.alo, self.ahi, self.xlo, middle,
                    self.slope, self.intercept, self.depth + 1),
                Box(self.alo, self.ahi, middle, self.xhi,
                    self.slope, self.intercept, self.depth + 1)]


def vertex_exponent(a, x, slope, intercept, witness):
    p, y, lam = witness
    return (a * (slope * x + intercept) + entropy_kl(a, p)
            + a * entropy_kl(x, y) + transfer_rate(p * y, lam) + .11 * lam)


def upper_for_witness(box: Box, witness):
    return max(vertex_exponent(a, x, box.slope, box.intercept, witness)
               for a, x in box.vertices())


def upper_for_bank(box: Box, bank):
    p, y, lam, rate = bank
    values = []
    for a, x in box.vertices():
        kl_a = a * (math.log(a) - np.log(p))
        if a < 1:
            kl_a += (1 - a) * (math.log1p(-a) - np.log1p(-p))
        kl_x = x * (math.log(x) - np.log(y)) + (1 - x) * (
            math.log1p(-x) - np.log1p(-y))
        values.append(a * (box.slope * x + box.intercept)
                      + kl_a + a * kl_x + rate + .11 * lam)
    maxima = np.maximum.reduce(values)
    index = int(np.argmin(maxima))
    return float(maxima[index]), index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-boxes", type=int, default=10_000)
    parser.add_argument("--max-depth", type=int, default=32)
    parser.add_argument("--margin", type=float, default=1e-7)
    args = parser.parse_args()
    witnesses = []
    for a in (1e-4, .001, .01, .1, .5, 1.):
        for x in (.104, .16, .25, .4, .48, .5, .52, .6, .75, .84, .896):
            row = optimize(a, x, .11)
            witnesses.append((row["p"], row["y"], row["lam"]))
    for a in (.9, .99, .999, 1.):
        for x in (.485, .49, .494, .498, .502, .506, .51, .515):
            row = optimize(a, x, .11)
            witnesses.append((row["p"], row["y"], row["lam"]))
    bank = tuple(np.asarray(values) for values in zip(*[
        (*w, transfer_rate(w[0] * w[1], w[2])) for w in witnesses]))
    pending = [Box(1e-4, 1., lo, hi, slope, intercept)
               for lo, hi, slope, intercept in segments()]
    accepted = processed = max_depth = 0
    max_accepted = -math.inf
    tightest = None
    worst = None
    while pending and processed < args.max_boxes:
        box = pending.pop()
        processed += 1
        max_depth = max(max_depth, box.depth)
        best = upper_for_bank(box, bank)
        if best[0] < -args.margin:
            accepted += 1
            if best[0] > max_accepted:
                max_accepted = best[0]
                tightest = {"alpha": [box.alo, box.ahi],
                            "x": [box.xlo, box.xhi], "witness_index": best[1]}
            continue
        if box.depth < args.max_depth:
            pending.extend(box.split())
            continue
        midpoint = ((box.alo + box.ahi) / 2, (box.xlo + box.xhi) / 2)
        proposal = optimize(*midpoint, .11)
        witness = (proposal["p"], proposal["y"], proposal["lam"])
        value = upper_for_witness(box, witness)
        if value < -args.margin:
            accepted += 1
            if value > max_accepted:
                max_accepted = value
                tightest = {"alpha": [box.alo, box.ahi],
                            "x": [box.xlo, box.xhi], "witness_index": "optimized"}
            bank = tuple(np.append(component, item) for component, item in
                         zip(bank, (*witness, transfer_rate(witness[0] * witness[1],
                                                            witness[2]))))
            continue
        worst = {"alpha": [box.alo, box.ahi],
                 "x": [box.xlo, box.xhi], "upper": min(best[0], value),
                 "midpoint_optimized": proposal["exponent"]}
        break
        if processed % 1000 == 0:
            print(json.dumps({"processed": processed, "pending": len(pending),
                              "max_accepted": max_accepted}), flush=True)
    print(json.dumps({"status": "covered_binary64" if not pending and not worst else "open",
                      "processed": processed, "accepted": accepted,
                      "pending": len(pending), "max_depth": max_depth,
                      "max_accepted": max_accepted, "unresolved": worst,
                      "tightest": tightest,
                      "next": ({"alpha": [pending[-1].alo, pending[-1].ahi],
                                "x": [pending[-1].xlo, pending[-1].xhi]}
                               if pending else None)}))


if __name__ == "__main__":
    main()
