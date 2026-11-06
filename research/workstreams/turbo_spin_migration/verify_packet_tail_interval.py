"""Outward-rounded 2D branch check for the 16% BAA packet-tail exponent.

Reuses the paper's Golay and first-accumulator interval bounds. A fixed
rational (x,u) witness is selected for each box and evaluated with mpmath
interval arithmetic. Prints the first unresolved box if a limit is hit.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
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
from research.workstreams.turbo_spin_migration.packet_tail_exponent import (
    final_packet_tail_exponent,
)


RHO = Fraction(4, 25)
ALPHA_LOWER = Fraction(1, 625)  # 8 eta with eta=1/5000


def rational_witness(beta: Fraction) -> tuple[Fraction, Fraction]:
    _, log_x, log_u = final_packet_tail_exponent(float(beta), float(RHO))
    x = Fraction(format(math.exp(log_x), ".17g"))
    u = Fraction(format(math.exp(log_u), ".17g"))
    if not (x > 0 and 0 < u <= 1):
        raise ValueError("optimizer produced invalid packet-tail witness")
    return x, u


def packet_term_upper(beta_lo: Fraction, beta_hi: Fraction,
                      x: Fraction, u: Fraction) -> float:
    xi, ui = as_interval(x), as_interval(u)
    x2, x4 = xi**2, xi**4
    a = 1 + ui * (6*x2 + x4)
    b = 4*ui*xi*(1+x2)
    c = xi + ui*(3*xi + 4*xi**3)
    d = ui*(1+6*x2+x4)
    root = (a+d+iv.sqrt((a-d)**2+4*b*c))/2
    beta = interval(beta_lo, beta_hi)
    objective = (iv.log(root)/4 - beta*iv.log(xi)
                 - as_interval(RHO)*iv.log(ui)/4)
    return upper_float(objective)


@dataclass(frozen=True)
class Box:
    alpha_lo: Fraction
    alpha_hi: Fraction
    beta_lo: Fraction
    beta_hi: Fraction

    def feasible(self) -> bool:
        return (self.beta_hi >= self.alpha_lo*HALF
                and self.beta_lo <= ONE-self.alpha_lo*HALF)

    def upper(self) -> float:
        if not self.feasible():
            return -math.inf
        beta_mid = (self.beta_lo+self.beta_hi)*HALF
        x, u = rational_witness(beta_mid)
        return add_upper(
            g_upper(self.alpha_lo, self.alpha_hi),
            p_upper(self.alpha_lo, self.alpha_hi,
                    self.beta_lo, self.beta_hi),
            packet_term_upper(self.beta_lo, self.beta_hi, x, u),
            negative_entropy_upper(self.beta_lo, self.beta_hi),
        )

    def split(self) -> tuple[Box, Box]:
        if self.alpha_hi-self.alpha_lo >= self.beta_hi-self.beta_lo:
            mid = (self.alpha_lo+self.alpha_hi)*HALF
            return (Box(self.alpha_lo, mid, self.beta_lo, self.beta_hi),
                    Box(mid, self.alpha_hi, self.beta_lo, self.beta_hi))
        mid = (self.beta_lo+self.beta_hi)*HALF
        return (Box(self.alpha_lo, self.alpha_hi, self.beta_lo, mid),
                Box(self.alpha_lo, self.alpha_hi, mid, self.beta_hi))

    def as_json(self) -> dict[str, list[str]]:
        return {"alpha": [str(self.alpha_lo), str(self.alpha_hi)],
                "beta": [str(self.beta_lo), str(self.beta_hi)]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-boxes", type=int, default=100_000)
    parser.add_argument("--max-depth", type=int, default=80)
    parser.add_argument("--progress-every", type=int, default=10_000)
    args = parser.parse_args()

    root = Box(ALPHA_LOWER, ONE, ZERO, ONE)
    heap = [(-root.upper(), 0, 0, root)]
    serial = processed = accepted = 0
    deepest = 0
    largest_accepted = -math.inf
    unresolved = None
    while heap and processed < args.max_boxes:
        neg_upper, depth, _, box = heapq.heappop(heap)
        upper = -neg_upper
        processed += 1
        deepest = max(deepest, depth)
        if upper < 0:
            accepted += 1
            largest_accepted = max(largest_accepted, upper)
            continue
        if depth >= args.max_depth:
            unresolved = {"upper": upper, "depth": depth,
                          "box": box.as_json()}
            break
        for child in box.split():
            value = child.upper()
            if value < 0:
                accepted += 1
                largest_accepted = max(largest_accepted, value)
            else:
                serial += 1
                heapq.heappush(heap, (-value, depth+1, serial, child))
        if args.progress_every and processed % args.progress_every == 0:
            print(json.dumps({"processed": processed,
                              "queued": len(heap),
                              "next_upper": -heap[0][0] if heap else None}),
                  flush=True)
    if unresolved is None and heap:
        upper, depth, _, box = heap[0]
        unresolved = {"upper": -upper, "depth": depth,
                      "box": box.as_json()}
    print(json.dumps({"result": "closed" if not heap and unresolved is None else "open",
                      "processed": processed, "accepted": accepted,
                      "deepest": deepest,
                      "max_accepted_upper": largest_accepted,
                      "unresolved": unresolved}, sort_keys=True))
    if heap or unresolved is not None:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
