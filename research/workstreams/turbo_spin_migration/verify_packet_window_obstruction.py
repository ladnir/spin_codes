"""Outward interval check for the explicit packet-window obstruction.

The combinatorial proof is in PACKET_WINDOW_OBSTRUCTION.md. This script
only certifies its strict numerical entropy inequality and displays the
critical window fraction for the paper's growth constant.
"""

from fractions import Fraction
from pathlib import Path
import sys

from mpmath import iv

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from research.workstreams.paper_architecture.certificates.single_sampled_ba_rm2sub.linear_time_audit.verify_outer_interval import (
    as_interval, lower_float, upper_float,
)


def main():
    iv.dps = 80
    a = as_interval(Fraction(13, 25))
    f = as_interval(Fraction(1, 5))
    divergence = a*iv.log(a/f)+(1-a)*iv.log((1-a)/(1-f))
    route_entropy = as_interval(Fraction(16, 39))*iv.log(iv.mpf(2))
    assert upper_float(divergence) < lower_float(route_entropy)
    c = as_interval(Fraction(39, 4))
    critical_window = (1-iv.sqrt(1-iv.power(2, -8/c)))/2
    target_window = as_interval(Fraction(11, 50))
    target_divergence = (iv.log(iv.mpf(1)/(4*target_window*(1-target_window))))/2
    required_c = 4*iv.log(iv.mpf(2))/target_divergence
    two_bit_window = (1-iv.sqrt(1-iv.power(2, -4/c)))/2
    print({"D(13/25 || 1/5) upper": upper_float(divergence),
           "(16/39) ln 2 lower": lower_float(route_entropy),
           "critical_window_enclosure": [lower_float(critical_window),
                                          upper_float(critical_window)],
           "distance_upper_enclosure": [lower_float(critical_window/2),
                                        upper_float(critical_window/2)],
           "two_bit_distance_upper_enclosure": [lower_float(two_bit_window/2),
                                                upper_float(two_bit_window/2)],
           "four_bit_growth_threshold_enclosure": [lower_float(required_c),
                                                   upper_float(required_c)]})


if __name__ == "__main__":
    main()
