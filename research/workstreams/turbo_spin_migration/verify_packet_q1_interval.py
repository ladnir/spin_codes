"""Outward-rounded branch check for Turbo's one-active-block exponent.

The dense-source domain uses alpha>=8/50000 and packet fraction>=19/100.
Sparse source messages require the separate simple count bound described
in Q1_PACKET_TRANSFER.md. This verifier checks the continuum exponent,
not the finite-step matrix correction.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import heapq
import json
import math
from pathlib import Path
import sys

from mpmath import iv
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from research.workstreams.paper_architecture.certificates.single_sampled_ba_rm2sub.linear_time_audit.verify_outer_interval import (
    HALF, ONE, ZERO, add_upper, as_interval, g_upper, interval,
    negative_entropy_upper, p_upper, upper_float,
)
from research.workstreams.turbo_spin_migration.packet_q1_exponent_screen import (
    inner_exponent,
)
from research.workstreams.turbo_spin_migration.packet_tail_exponent import (
    packet_transfer_radius,
)


RHO_FLOOR = Fraction(19, 100)
ALPHA_FLOOR = Fraction(1, 6250)
DELTA = Fraction(11, 100)
SCHEDULE = Fraction(39, 4)
M = (1 << 22) - 1
P0 = Fraction(1 << 21, M)


@lru_cache(maxsize=100_000)
def outer_witness(beta: Fraction, rho: Fraction) -> tuple[Fraction, Fraction]:
    b, r = float(beta), float(rho)

    def value(point):
        log_x, log_u = point
        return (math.log(packet_transfer_radius(math.exp(log_x), math.exp(log_u)))/4
                -b*log_x-r*log_u/4)

    fit = minimize(value, np.asarray([0.0, 0.0]), method="L-BFGS-B",
                   bounds=[(-30, 30), (-30, 30)])
    return (Fraction(format(math.exp(fit.x[0]), ".17g")),
            Fraction(format(math.exp(fit.x[1]), ".17g")))


@lru_cache(maxsize=10_000)
def inner_witness(rho: Fraction) -> tuple[Fraction, Fraction]:
    _, log_y, gamma = inner_exponent(float(rho), float(DELTA))
    return (Fraction(format(math.exp(log_y), ".17g")),
            Fraction(format(gamma, ".17g")))


def perron_2x2(a, b, c, d):
    return (a+d+iv.sqrt((a-d)**2+4*b*c))/2


def packet_upper(beta_lo: Fraction, beta_hi: Fraction,
                 rho_lo: Fraction, rho_hi: Fraction,
                 x: Fraction, u: Fraction) -> float:
    xi, ui = as_interval(x), as_interval(u)
    x2, x4 = xi**2, xi**4
    root = perron_2x2(
        1+ui*(6*x2+x4),
        4*ui*xi*(1+x2),
        xi+ui*(3*xi+4*xi**3),
        ui*(1+6*x2+x4),
    )
    beta = interval(beta_lo, beta_hi)
    rho = interval(rho_lo, rho_hi)
    return upper_float(iv.log(root)/4 - beta*iv.log(xi)
                       -rho*iv.log(ui)/4)


def inner_upper(rho_lo: Fraction, rho_hi: Fraction,
                y: Fraction, gamma: Fraction) -> float:
    yi, gi = as_interval(y), as_interval(gamma)
    E = iv.exp(-gi)
    J = (1-E)/gi
    root = perron_2x2(
        iv.mpf(1), yi*J, yi*J/M,
        E*(1+yi*(1-iv.mpf(1)/M)),
    )
    rho = interval(rho_lo, rho_hi)
    value = (iv.log(root)-rho*iv.log(yi)
             +gi*as_interval(DELTA)/as_interval(P0))/4
    return add_upper(upper_float(value),
                     negative_entropy_upper(rho_lo, rho_hi)/4)


@dataclass(frozen=True)
class Box:
    alpha_lo: Fraction
    alpha_hi: Fraction
    beta_lo: Fraction
    beta_hi: Fraction
    rho_lo: Fraction
    rho_hi: Fraction

    def feasible(self) -> bool:
        return (self.beta_hi >= self.alpha_lo*HALF
                and self.beta_lo <= ONE-self.alpha_lo*HALF)

    def upper(self) -> float:
        if not self.feasible():
            return -math.inf
        beta_mid = (self.beta_lo+self.beta_hi)*HALF
        rho_mid = (self.rho_lo+self.rho_hi)*HALF
        x, u = outer_witness(beta_mid, rho_mid)
        y, gamma = inner_witness(rho_mid)
        return add_upper(
            upper_float(iv.log(as_interval(Fraction(2)))
                        /as_interval(SCHEDULE)),
            g_upper(self.alpha_lo, self.alpha_hi),
            p_upper(self.alpha_lo, self.alpha_hi,
                    self.beta_lo, self.beta_hi),
            packet_upper(self.beta_lo, self.beta_hi,
                         self.rho_lo, self.rho_hi, x, u),
            negative_entropy_upper(self.beta_lo, self.beta_hi),
            inner_upper(self.rho_lo, self.rho_hi, y, gamma),
        )

    def split(self) -> tuple[Box, Box]:
        widths = (self.alpha_hi-self.alpha_lo,
                  self.beta_hi-self.beta_lo,
                  self.rho_hi-self.rho_lo)
        dim = max(range(3), key=lambda i: widths[i])
        if dim == 0:
            mid = (self.alpha_lo+self.alpha_hi)*HALF
            return (Box(self.alpha_lo, mid, self.beta_lo, self.beta_hi,
                        self.rho_lo, self.rho_hi),
                    Box(mid, self.alpha_hi, self.beta_lo, self.beta_hi,
                        self.rho_lo, self.rho_hi))
        if dim == 1:
            mid = (self.beta_lo+self.beta_hi)*HALF
            return (Box(self.alpha_lo, self.alpha_hi, self.beta_lo, mid,
                        self.rho_lo, self.rho_hi),
                    Box(self.alpha_lo, self.alpha_hi, mid, self.beta_hi,
                        self.rho_lo, self.rho_hi))
        mid = (self.rho_lo+self.rho_hi)*HALF
        return (Box(self.alpha_lo, self.alpha_hi, self.beta_lo, self.beta_hi,
                    self.rho_lo, mid),
                Box(self.alpha_lo, self.alpha_hi, self.beta_lo, self.beta_hi,
                    mid, self.rho_hi))

    def as_json(self):
        return {"alpha": [str(self.alpha_lo), str(self.alpha_hi)],
                "beta": [str(self.beta_lo), str(self.beta_hi)],
                "rho": [str(self.rho_lo), str(self.rho_hi)]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-boxes", type=int, default=100_000)
    parser.add_argument("--max-depth", type=int, default=80)
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()
    root = Box(ALPHA_FLOOR, ONE, ZERO, ONE, RHO_FLOOR, ONE)
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
