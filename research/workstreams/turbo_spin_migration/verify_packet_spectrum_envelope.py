"""Certified interval upper envelope for the expected BAA packet spectrum.

For each packet-fraction band, branch over Golay and first-accumulator
weights. Every leaf upper-bounds the asymptotic expected number of local
nonzero words in that band, normalized by the outer length b. The JSON
output is a reusable certificate input for occupation-two bounds.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import heapq
import json
import math

from verify_packet_q2_interval import outer_upper
from research.workstreams.paper_architecture.certificates.single_sampled_ba_rm2sub.linear_time_audit.verify_outer_interval import (
    HALF, ONE, ZERO,
)


def region_upper(rho_lo: Fraction, rho_hi: Fraction,
                 max_boxes: int) -> tuple[float, int, tuple]:
    """Return a valid bound even if the requested refinement budget ends."""
    root = ((ZERO, ONE), (ZERO, ONE))

    def feasible(bounds):
        alpha, beta = bounds
        return beta[1] >= alpha[0]*HALF and beta[0] <= ONE-alpha[0]*HALF

    def bound(bounds):
        if not feasible(bounds):
            return -math.inf
        return outer_upper(bounds[0], bounds[1], (rho_lo, rho_hi))

    heap = [(-bound(root), 0, root)]
    serial = processed = 0
    while heap and processed < max_boxes:
        neg_value, _, bounds = heapq.heappop(heap)
        processed += 1
        alpha, beta = bounds
        dimension = 0 if alpha[1]-alpha[0] >= beta[1]-beta[0] else 1
        lo, hi = bounds[dimension]
        mid = (lo+hi)*HALF
        for part in ((lo, mid), (mid, hi)):
            child = list(bounds)
            child[dimension] = part
            child = tuple(child)
            value = bound(child)
            if value != -math.inf:
                serial += 1
                heapq.heappush(heap, (-value, serial, child))
    # All live leaves partition the feasible domain, including their
    # boundaries. The largest leaf bound therefore bounds the supremum.
    return (-heap[0][0] if heap else -math.inf), processed, (heap[0][2] if heap else ())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bands", type=int, default=32)
    parser.add_argument("--max-boxes-per-band", type=int, default=1000)
    parser.add_argument("--rho-floor", type=Fraction,
                        default=Fraction(19, 100))
    parser.add_argument("--output", default="")
    args = parser.parse_args()
    if args.bands < 1 or args.max_boxes_per_band < 1:
        raise ValueError("positive band and box counts required")
    records = []
    for i in range(args.bands):
        lo = args.rho_floor+(ONE-args.rho_floor)*i/args.bands
        hi = args.rho_floor+(ONE-args.rho_floor)*(i+1)/args.bands
        upper, processed, top = region_upper(lo, hi, args.max_boxes_per_band)
        records.append({"rho_lo": str(lo), "rho_hi": str(hi),
                        "upper_nats_per_bit": upper,
                        "boxes_processed": processed,
                        "top_box": [[str(x) for x in pair] for pair in top]})
        print(json.dumps({"band": i, **records[-1]}), flush=True)
    result = {"rho_floor": str(args.rho_floor), "bands": records,
              "scope": "outward-rounded asymptotic expected packet-spectrum upper envelope"}
    if args.output:
        with open(args.output, "w", encoding="utf-8") as out:
            json.dump(result, out, indent=2)
            out.write("\n")


if __name__ == "__main__":
    main()
