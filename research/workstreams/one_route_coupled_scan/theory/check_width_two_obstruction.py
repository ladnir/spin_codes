"""Exact arithmetic for the analytic width-two reset obstruction.

This checks the scalar inequalities in WIDTH_TWO_OBSTRUCTION.md. It does
not simulate routes or numerically certify a distance lower bound. No files
are written. Run with python -B.
"""

from fractions import Fraction
from math import factorial


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def harmonic(n):
    return sum((Fraction(1, j) for j in range(1, n + 1)), Fraction(0))


def main():
    regions, batches, selected, width = 256, 4096, 7, 2
    dimension = 32 * 128
    require(regions * selected * width < dimension,
            "reset constraints exhaust the message dimension")
    for occupied in (16, 17):
        mean_lower = (Fraction(batches - occupied, batches)
                      * Fraction(selected, occupied + 1)
                      * (1 + harmonic(occupied) - harmonic(selected)))
        require(mean_lower > Fraction(71, 100), "insufficient deletion mean")
        print(f"occupied={occupied}; exact_mean_lower={mean_lower}; "
              f"decimal={float(mean_lower):.12f}")

    deviation = Fraction(71, 100) - Fraction(65, 100)
    exponent = 2 * regions * deviation * deviation
    require(exponent == Fraction(1152, 625), "incorrect Hoeffding exponent")
    # exp(x) exceeds every finite nonnegative Taylor sum for x>0.
    exponential_lower = sum((exponent**j / factorial(j) for j in range(14)),
                            Fraction(0))
    tail_upper = 1 / exponential_lower
    require(tail_upper < Fraction(159, 1000), "tail bound did not close")
    require(1 - tail_upper > Fraction(841, 1000), "failure lower bound too small")
    print(f"constraints={regions * selected * width}; dimension={dimension}")
    print(f"rational_tail_upper={float(tail_upper):.12f}; "
          f"bad_setup_probability_lower={float(1 - tail_upper):.12f}")
    print("PASS: analytic width-two failure lower bound > 0.841")


if __name__ == "__main__":
    main()
