"""Bounded exact checks for WIDTH_RESET_BOUNDS.md; no files or simulations.

Checks pair-slot rounding on all small gaps, exact tiny-route expectations,
and the conservative finite K20 width-four rejection. Uses only stdlib.
Run with python -B. These are obstruction checks, not distance certificates.
"""

from fractions import Fraction as F
from itertools import combinations


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def harmonic(n):
    return sum((F(1, i) for i in range(1, n + 1)), F(0))


def deleted_mean_lower(slots, occupied, selected, width):
    half = width // 2
    return (F(slots - occupied, slots) * F(selected, occupied + 1)
            * (1 + harmonic(occupied) - harmonic(selected))
            - F(2 * (half - 1) * selected, slots))


def kept_slots(slots, half, left, right):
    require(slots % half == 0, "unaligned region")
    count = max(0, right // half - left // half - 1)
    brute = sum(1 for start in range(0, slots, half)
                if start > left and start + half <= right)
    require(count == brute, ("rounding formula", slots, half, left, right))
    gap = right - left - 1
    require(half * count >= gap - 2 * (half - 1), "rounding lower bound")
    return half * count


def tiny_checks():
    rounding = 0
    for half in (1, 2, 4):
        for slots in range(half, 33, half):
            for left in range(slots):
                for right in range(left + 1, slots + 1):
                    kept_slots(slots, half, left, right)
                    rounding += 1
    profiles = 0
    for width in (2, 4, 8):
        half = width // 2
        for slots in (4, 8):
            for occupied in range(1, slots):
                for selected in sorted({1, occupied, max(1, occupied // width)}):
                    total, routes = 0, 0
                    for positions in combinations(range(slots), occupied):
                        intervals = list(zip(positions, positions[1:] + (slots,)))
                        intervals.sort(key=lambda pair: (-(pair[1] - pair[0] - 1), pair[0]))
                        kept = sum(kept_slots(slots, half, *pair)
                                   for pair in intervals[:selected])
                        total += kept
                        routes += 1
                    require(F(total, routes * slots) >=
                            deleted_mean_lower(slots, occupied, selected, width),
                            ("tiny exact expectation", width, slots, occupied, selected))
                    profiles += 1
    return rounding, profiles


def main():
    rounding, profiles = tiny_checks()
    m = 4096
    require(harmonic(64) - harmonic(16) > F(34, 25), "harmonic 64/16")
    require(harmonic(65) - harmonic(16) > F(137, 100), "harmonic 65/16")
    require(harmonic(64) - harmonic(15) > F(71, 50), "harmonic 64/15")
    require(deleted_mean_lower(m, 64, 16, 4) > F(557, 1000), "even mean")
    require(deleted_mean_lower(m, 65, 16, 4) > F(557, 1000), "odd mean")
    require(deleted_mean_lower(m, 64, 15, 4) > F(53, 100), "exceptional mean")
    require((255 * F(557, 1000) + F(53, 100)) / 256 > F(139, 250),
            "conservative average")
    require(4 * (255 * 16 + 15) < 256 * 64, "dimension budget")
    exponent = 2 * 256 * (F(139, 250) - F(55, 100))**2
    require(exponent == F(1152, 62500), "Hoeffding exponent")
    # exp(x)>1+x, hence 1-exp(-x)>x/(1+x).
    failure_lower = exponent / (1 + exponent)
    require(failure_lower > F(18, 1000), "width-four failure lower bound")
    for width, selected in ((4, 16), (8, 8)):
        average = (127 * deleted_mean_lower(m, 64, selected, width)
                   + 128 * deleted_mean_lower(m, 65, selected, width)
                   + deleted_mean_lower(m, 64, selected - 1, width)) / 256
        print(f"w={width}; exact-expression mean-deletion lower={float(average):.12f}")
    print(f"rounding cases={rounding}; exact tiny-route expectation profiles={profiles}")
    print(f"w4 at 45%: conservative failure lower > {float(failure_lower):.12f}")
    print("PASS: finite reset obstruction checks; no distance lower-bound certificate")


if __name__ == "__main__":
    main()
