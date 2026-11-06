"""Exploratory Q1 screen with a separate Chernoff tilt per packet count.

The finite transfer, when selected, is a valid analytic upper envelope
for ideal setup. The printed binary64 values are not interval certificates.
"""

import argparse
from math import ceil, exp, log, log2

from packet_baa_q1_screen import (
    LIVE_BIT_DENSITY, accumulator_packet_counts, expected_packet_spectrum,
    one_block_moment,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", type=int, default=108)
    parser.add_argument("--rho", type=float, default=0.16)
    parser.add_argument("--delta", type=float, default=0.11)
    parser.add_argument("--input-bits", type=int, default=1 << 20)
    parser.add_argument("--outer-blocks", type=int, default=0)
    parser.add_argument("--grid-steps", type=int, default=240)
    args = parser.parse_args()
    if (args.regions < 1 or (4*args.regions) % 24
            or args.grid_steps < 2):
        raise ValueError("invalid regions or grid")
    if args.outer_blocks and args.outer_blocks % 16:
        raise ValueError("finite outer blocks must be a multiple of 16")

    R = args.regions
    b = 4*R
    L = args.outer_blocks or args.input_bits / (b/2)
    spectrum = expected_packet_spectrum(R, accumulator_packet_counts(R))
    floor = ceil(args.rho*R)
    best = [float("inf")] * (R+1)
    selected_gamma = [0.0] * (R+1)
    for i in range(args.grid_steps+1):
        gamma = 0.025 + i * (8 - 0.025) / args.grid_steps
        moment = one_block_moment(R, gamma, args.outer_blocks or None)
        theta = (gamma/(4*L*LIVE_BIT_DENSITY)
                 if args.outer_blocks else 0.0)
        prefactor = exp(gamma*args.delta*R/LIVE_BIT_DENSITY)
        for v in range(floor, R+1):
            value = prefactor * moment[v] * exp(128*theta*v)
            if value < best[v]:
                best[v] = value
                selected_gamma[v] = gamma
    charges = [(spectrum[v]*best[v], v) for v in range(floor, R+1)]
    total = sum(value for value, _ in charges)
    charges.sort(reverse=True)
    print({"b": b, "L": L, "packet_floor": floor,
           "log2_q1_upper": log2(L*total),
           "log2_per_block_charge": log2(total),
           "schedule_margin_nats_per_bit":
               (log(total)/b + log(2)/(39/4)),
           "dominant_packet_counts":
               [(v, round(value/total, 5), round(selected_gamma[v], 4))
                for value, v in charges[:8]],
           "scope": "binary64 selected-tilt screen, not a certificate"})


if __name__ == "__main__":
    main()
