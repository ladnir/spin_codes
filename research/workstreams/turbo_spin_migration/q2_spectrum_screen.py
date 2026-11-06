"""Binary64 occupation-two exponent screen using exact finite BAA spectrum.

The b-normalized spectrum is a finite-size proxy for the asymptotic
packet exponent. This is not a distance certificate.
"""

import argparse
from math import ceil, log

from packet_baa_q1_screen import accumulator_packet_counts, expected_packet_spectrum
from q2_inner_exponent_screen import screen


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", type=int, default=108)
    parser.add_argument("--rho-floor", type=float, default=0.19)
    parser.add_argument("--step", type=int, default=8)
    args = parser.parse_args()
    R = args.regions
    if R < 1 or 4*R % 24 or args.step < 1:
        raise ValueError("invalid regions or step")
    b = 4*R
    spectrum = expected_packet_spectrum(R, accumulator_packet_counts(R))
    floor = ceil(args.rho_floor*R)
    selected = sorted(set([floor, floor+1, floor+2, R]
                          + list(range(floor+args.step, R, args.step))))
    outcomes = []
    for i, v1 in enumerate(selected):
        for v2 in selected[i:]:
            if not (spectrum[v1] and spectrum[v2]):
                continue
            inner = screen(v1/R, v2/R, 0.11, 39/4)
            total = (inner["schedule_plus_inner_nats_per_bit"]
                     +(log(spectrum[v1])+log(spectrum[v2]))/b)
            outcomes.append((total, v1, v2,
                             inner["schedule_plus_inner_nats_per_bit"],
                             (log(spectrum[v1])+log(spectrum[v2]))/b))
    outcomes.sort(reverse=True)
    print({"b": b, "pairs_checked": len(outcomes),
           "worst_pairs": outcomes[:12],
           "scope": "binary64 finite-spectrum proxy; not a certificate"})


if __name__ == "__main__":
    main()
