"""Outward-rounded six-variable Q2 exponent branch check.

Uses the selected shared-outer per-packet majorant and the continuum
two-packet region kernel. Finite-step collision control is separate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import heapq
import json
import math
from pathlib import Path
import sys

from mpmath import iv

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from research.workstreams.paper_architecture.certificates.single_sampled_ba_rm2sub.linear_time_audit.verify_outer_interval import (
    HALF, ONE, ZERO, add_upper, as_interval, g_upper, interval,
    negative_entropy_upper, p_upper, upper_float,
)
from research.workstreams.turbo_spin_migration.verify_packet_q1_interval import (
    DELTA, P0, SCHEDULE, M, outer_witness, packet_upper, perron_2x2,
)
from research.workstreams.turbo_spin_migration.q2_inner_exponent_screen import (
    screen as q2_screen,
)


FLOOR = Fraction(19, 100)


@lru_cache(maxsize=100_000)
def inner_witness(rho1: Fraction, rho2: Fraction):
    result = q2_screen(float(rho1), float(rho2), float(DELTA),
                       float(SCHEDULE))
    return tuple(Fraction(format(result[key], ".17g"))
                 for key in ("u", "v", "gamma"))


def inner_upper(rho1, rho2, y1, y2, gamma) -> float:
    yi, zi, gi = map(as_interval, (y1, y2, gamma))
    E = iv.exp(-gi)
    J = (1-E)/gi
    H = 2*(gi-1+E)/(gi**2)
    G = 2*(1-(1+gi)*E)/(gi**2)
    r = iv.mpf(1)/M
    a = 1+yi*zi*r*H
    b = (yi+zi)*J+yi*zi*(1-r)*G
    c = (yi+zi)*r*J+yi*zi*r*(1-r)*G
    d = E+(yi+zi)*(1-r)*E+yi*zi*(r*G+(1-r)**2*E)
    root = perron_2x2(a, b, c, d)
    p = interval(*rho1)
    q = interval(*rho2)
    value = (iv.log(root)-p*iv.log(yi)-q*iv.log(zi)
             +gi*as_interval(DELTA)/as_interval(P0))/4
    return add_upper(
        upper_float(value),
        negative_entropy_upper(*rho1)/4,
        negative_entropy_upper(*rho2)/4,
    )


def outer_upper(alpha, beta, rho) -> float:
    beta_mid = sum(beta)*HALF
    rho_mid = sum(rho)*HALF
    x, u = outer_witness(beta_mid, rho_mid)
    return add_upper(
        g_upper(*alpha), p_upper(*alpha, *beta),
        packet_upper(*beta, *rho, x, u),
        negative_entropy_upper(*beta),
    )


class Box:
    def __init__(self, bounds):
        self.bounds = tuple(bounds)

    def feasible(self) -> bool:
        for i in (0, 3):
            a, b = self.bounds[i], self.bounds[i+1]
            if b[1] < a[0]*HALF or b[0] > ONE-a[0]*HALF:
                return False
        return True

    def upper(self) -> float:
        if not self.feasible():
            return -math.inf
        a1, b1, r1, a2, b2, r2 = self.bounds
        y1, y2, gamma = inner_witness(sum(r1)*HALF, sum(r2)*HALF)
        return add_upper(
            upper_float(2*iv.log(iv.mpf(2))/as_interval(SCHEDULE)),
            outer_upper(a1, b1, r1),
            outer_upper(a2, b2, r2),
            inner_upper(r1, r2, y1, y2, gamma),
        )

    def split(self):
        widths = [hi-lo for lo, hi in self.bounds]
        dim = max(range(6), key=lambda i: widths[i])
        lo, hi = self.bounds[dim]
        mid = (lo+hi)*HALF
        left, right = list(self.bounds), list(self.bounds)
        left[dim], right[dim] = (lo, mid), (mid, hi)
        return Box(left), Box(right)

    def as_json(self):
        return [[str(lo), str(hi)] for lo, hi in self.bounds]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-boxes", type=int, default=100_000)
    parser.add_argument("--max-depth", type=int, default=80)
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()
    root = Box(((ZERO, ONE), (ZERO, ONE), (FLOOR, ONE),
                (ZERO, ONE), (ZERO, ONE), (FLOOR, ONE)))
    heap = [(-root.upper(), 0, 0, root)]
    serial = processed = accepted = deepest = 0
    max_accepted = -math.inf
    unresolved = None
    while heap and processed < args.max_boxes:
        neg_upper, depth, _, box = heapq.heappop(heap)
        upper = -neg_upper
        processed += 1
        deepest = max(deepest, depth)
        if upper < 0:
            accepted += 1
            max_accepted = max(max_accepted, upper)
            continue
        if depth >= args.max_depth:
            unresolved = {"upper": upper, "depth": depth,
                          "box": box.as_json()}
            break
        for child in box.split():
            value = child.upper()
            if value < 0:
                accepted += 1
                max_accepted = max(max_accepted, value)
            else:
                serial += 1
                heapq.heappush(heap, (-value, depth+1, serial, child))
        if args.progress_every and processed % args.progress_every == 0:
            print(json.dumps({"processed": processed, "queued": len(heap),
                              "next_upper": -heap[0][0] if heap else None}),
                  flush=True)
    if unresolved is None and heap:
        neg_upper, depth, _, box = heap[0]
        unresolved = {"upper": -neg_upper, "depth": depth,
                      "box": box.as_json()}
    closed = not heap and unresolved is None
    print(json.dumps({"result": "closed" if closed else "open",
                      "processed": processed, "accepted": accepted,
                      "deepest": deepest,
                      "max_accepted_upper": max_accepted,
                      "unresolved": unresolved}, sort_keys=True))
    if not closed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
