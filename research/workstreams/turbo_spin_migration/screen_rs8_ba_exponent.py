"""Non-certified limiting BAA spectrum screen for RS8 versus Golay bases.

This uses floating-point global optimization for discovery. It is not an
outward interval certificate and must not be cited as a proved bound.
"""

from math import exp, log

import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar

from rs8_base_spectrum import row_spectrum


def entropy(x: float) -> float:
    if x <= 0 or x >= 1:
        return 0.0
    return -x * log(x) - (1 - x) * log(1 - x)


def transition_exponent(a: float, b: float) -> float:
    if not 0 <= a <= 1 or not a / 2 <= b <= 1 - a / 2:
        return float("-inf")
    if b <= 0 or b >= 1:
        return 0.0
    return (b * entropy(a / (2 * b))
            + (1 - b) * entropy(a / (2 * (1 - b))) - entropy(a))


def base_exponent(counts: list[int], a: float) -> float:
    if a <= 0 or a >= 1:
        return 0.0
    weights = np.flatnonzero(counts).astype(float)
    log_counts = np.log(np.asarray(counts)[weights.astype(int)])
    length = len(counts) - 1

    def objective(t: float) -> float:
        terms = log_counts + weights * t
        maximum = float(max(terms))
        log_sum = maximum + log(sum(exp(v - maximum) for v in terms))
        return log_sum / length - a * t

    return float(minimize_scalar(objective, bounds=(-25, 25), method="bounded",
                                 options={"xatol": 1e-8}).fun)


def baa_exponent(counts: list[int], output_weight: float) -> tuple[float, float, float]:
    def negative(value: np.ndarray) -> float:
        a, t = map(float, value)
        c = a / 2 + (1 - a) * t
        return -(base_exponent(counts, a)
                 + transition_exponent(a, c)
                 + transition_exponent(c, output_weight))

    answer = differential_evolution(negative, [(0, 1), (0, 1)],
                                    maxiter=80, popsize=8, seed=19,
                                    polish=True, tol=1e-8, workers=1)
    a, t = map(float, answer.x)
    return -float(answer.fun), a, a / 2 + (1 - a) * t


if __name__ == "__main__":
    rs8 = row_spectrum()
    golay = [0] * 25
    for weight, count in ((0, 1), (8, 759), (12, 2576), (16, 759), (24, 1)):
        golay[weight] = count
    for output_weight in (0.104, 0.11, 0.2):
        for name, counts in (("RS8 row", rs8), ("Golay", golay)):
            value, a, c = baa_exponent(counts, output_weight)
            print(f"weight={output_weight:.3f} base={name:7s} "
                  f"exponent={value:.9f} a={a:.5f} c={c:.5f}", flush=True)
