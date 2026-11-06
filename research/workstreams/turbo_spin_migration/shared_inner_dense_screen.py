"""Binary64 diagnostic for the scalable one-bit route with Turbo's s22 inner.

This checks points of the positive-occupation exponent, not a certificate.
The 2-state transfer is an entrywise upper bound for every entering state.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit


HERE = Path(__file__).resolve().parent
OUTER = HERE.parent / "inner_design/imt_asymptotic/d11/OUTER_REFINED.json"
WIDTH = 64
STATE = 22
M = (1 << STATE) - 1
SPECTRUM = {
    0: 1,
    16: 2604,
    24: 291648,
    28: 888832,
    32: 1828134,
    36: 888832,
    40: 291648,
    48: 2604,
    64: 1,
}


def dual_spectrum() -> np.ndarray:
    """MacWilliams transform of the exact s22 expansion spectrum."""
    assert sum(SPECTRUM.values()) == 1 << STATE
    counts = []
    for weight in range(WIDTH + 1):
        total = 0
        for image_weight, multiplicity in SPECTRUM.items():
            coefficient = sum(
                (-1) ** j * math.comb(image_weight, j)
                * math.comb(WIDTH - image_weight, weight - j)
                for j in range(max(0, weight - (WIDTH - image_weight)),
                               min(weight, image_weight) + 1)
            )
            total += multiplicity * coefficient
        quotient, remainder = divmod(total, 1 << STATE)
        assert remainder == 0 and quotient >= 0
        counts.append(quotient)
    assert counts[0] == 1 and not any(counts[1:8])
    assert sum(counts) == 1 << (WIDTH - STATE)
    return np.asarray(counts, dtype=np.float64)


DUAL = dual_spectrum()
WEIGHTS = np.asarray([w for w in SPECTRUM if w], dtype=np.float64)


def entropy_kl(x: float, y: float) -> float:
    if x == 0:
        return -math.log1p(-y)
    if x == 1:
        return -math.log(y)
    return x * math.log(x / y) + (1 - x) * math.log((1 - x) / (1 - y))


@lru_cache(maxsize=128)
def outer_majorant(x: float) -> float:
    payload = json.loads(OUTER.read_text(encoding="utf-8"))
    assert payload["status"] == "proved"
    left = payload["left_segments"]
    half = min(x, 1 - x)
    for segment in left:
        lo, hi = (float(Fraction(v)) for v in segment["weight"])
        if lo - 1e-14 <= half <= hi + 1e-14:
            return (float(Fraction(segment["slope"])) * half
                    + float(Fraction(segment["intercept"])))
    assert float(Fraction(left[-1]["weight"][1])) - 1e-14 <= half <= .5
    return float(Fraction(payload["central_constant_upper"]))


@lru_cache(maxsize=100_000)
def transfer_rate(beta: float, lam: float) -> float:
    z = math.exp(-lam)
    g0 = 1 - beta + beta * z
    g1 = beta + (1 - beta) * z
    live_values = g0 ** (WIDTH - WEIGHTS) * g1 ** WEIGHTS
    fresh_live = float(np.max(live_values))
    uniform_live = sum(SPECTRUM[int(w)] * value
                       for w, value in zip(WEIGHTS, live_values)) / M
    j = np.arange(WIDTH + 1)
    zero_zero = float(np.dot(DUAL, (beta * z) ** j
                            * (1 - beta) ** (WIDTH - j)))
    zero_total = g0 ** WIDTH
    zero_live = max(0.0, zero_total - zero_zero)
    # A live state refreshes uniformly over nonzero states. Track an
    # arbitrary freshly activated state for one step, then dominate its
    # successor distribution by uniform live mass. The uniform-live
    # row uses the exact average of the map spectrum. Returning to zero
    # requires C X != 0; continuation is deliberately overcounted.
    kernel_probability = float(np.dot(DUAL, beta ** j
                                      * (1 - beta) ** (WIDTH - j)))
    active_probability = max(0.0, 1 - kernel_probability)
    matrix = np.asarray([
        [zero_zero, zero_live, 0.0],
        [min(fresh_live, active_probability) / M, 0.0, fresh_live],
        [min(uniform_live, active_probability) / M, 0.0, uniform_live],
    ])
    radius = float(max(abs(np.linalg.eigvals(matrix))))
    return math.log(radius) / WIDTH


def exponent(alpha: float, x: float, p: float, y: float,
             lam: float, delta: float) -> float:
    return (alpha * outer_majorant(x)
            + entropy_kl(alpha, p) + alpha * entropy_kl(x, y)
            + transfer_rate(p * y, lam) + delta * lam)


def optimize(alpha: float, x: float, delta: float) -> dict:
    q = alpha * x
    p0 = min(1 - 1e-12, .5 + .5 * alpha * (1 - x) / (1 - q))
    starts = [(p0, .5 / p0, math.log((1 - delta) / delta)),
              (max(alpha, 1e-5), x, max(alpha, 1e-5)),
              (min(.999, 2 * alpha), .4, max(alpha, 1e-5)),
              (.8, .5, 1.8)]
    best = None
    for p, y, lam in starts:
        p, y = np.clip([p, y], 1e-10, 1 - 1e-10)
        point = np.asarray([math.log(p / (1 - p)), math.log(y / (1 - y)),
                            math.log(lam)])

        def decode(v):
            p1, y1 = np.clip(expit(v[:2]), 1e-12, 1 - 1e-12)
            return float(p1), float(y1), math.exp(float(np.clip(v[2], -25, 2.4)))

        fit = minimize(lambda v: exponent(alpha, x, *decode(v), delta),
                       point, method="Nelder-Mead",
                       options={"maxiter": 280, "xatol": 1e-8, "fatol": 1e-12})
        if best is None or fit.fun < best["exponent"]:
            bp, by, blam = decode(fit.x)
            best = {"alpha": alpha, "x": x, "exponent": float(fit.fun),
                    "p": bp, "y": by, "lam": blam,
                    "optimizer_success": bool(fit.success)}
    assert best is not None
    return best


def fixed_one_two(delta: float = .11, growth: float = 39 / 4) -> list[dict]:
    """Floating continuum gates copied algebraically from the IMT proof."""
    r = 1 / M
    p0 = (1 << (STATE - 1)) / M
    rows = []
    for q in (1, 2):
        gamma = q * math.log(2)
        theta = gamma / p0
        a = 1 / (1 << q)
        f = (1 - a) / gamma
        gap = 2 * (gamma - 1 + a) / gamma**2
        edge = 2 * (1 - (1 + gamma) * a) / gamma**2
        k0 = np.asarray([[1., 0.], [0., a]])
        k1 = np.asarray([[0., f], [r * f, a]])
        k2 = np.asarray([[r * gap, edge], [r * edge, r * edge + a]])
        matrix = k0 + q * k1 + (k2 if q == 2 else 0)
        radius = float(max(abs(np.linalg.eigvals(matrix))))
        exponent = (math.log(radius) + delta * theta - q * math.log(2) / 2
                    + q * .01281 + q * math.log(2) / growth)
        rows.append({"q": q, "exponent_per_outer_bit": exponent,
                     "exponent_per_active_row_bit": exponent / q})
    return rows


def sparse_gate() -> dict:
    """Exact rational Collatz inequalities; log exponent is diagnostic."""
    f = Fraction
    alpha_max = f(1, 10_000)
    beta_scale, tilt_scale = f(4, 5), f(8, 5)
    live_scale, fresh_slope, contraction = f(1, 1024), f(24), f(48)
    # v=(1,live_scale*(1+fresh_slope*alpha),live_scale).
    # Bound T00 using kernel distance 8 and Bonferroni; all other
    # entries use union bounds and the expansion's minimum/mean weight.
    residual_zero = (
        -64 * beta_scale + math.comb(64, 2) * beta_scale**2 * alpha_max
        + f(math.comb(64, 8)) * beta_scale**8 * alpha_max**7
        + 64 * beta_scale * live_scale * (1 + fresh_slope * alpha_max)
        + contraction
    )
    residual_fresh = (
        64 * beta_scale / (M * live_scale) - 16 * tilt_scale
        - fresh_slope + contraction
        + (16 * beta_scale * tilt_scale
           + math.comb(16, 2) * tilt_scale**2
           + contraction * fresh_slope) * alpha_max
    )
    residual_uniform = (
        64 * beta_scale / (M * live_scale) - 32 * tilt_scale
        + contraction
        + (64 * beta_scale * tilt_scale
           + math.comb(64, 2) * tilt_scale**2) * alpha_max
    )
    assert max(residual_zero, residual_fresh, residual_uniform) < 0

    payload = json.loads(OUTER.read_text(encoding="utf-8"))
    endpoints = {f(13, 125), f(112, 125), f(1, 2)}
    for segment in payload["left_segments"]:
        for raw in segment["weight"]:
            x = f(raw)
            endpoints.add(x)
            endpoints.add(1 - x)

    def row_cost(x: float) -> float:
        entropy = -x * math.log(x) - (1 - x) * math.log(1 - x)
        return outer_majorant(x) + math.log(2) - entropy

    worst_x = max(endpoints, key=lambda x: row_cost(float(x)))
    a = float(alpha_max)
    support = math.log(2) / (39 / 4)
    exponent_per_alpha = (row_cost(float(worst_x))
                          - math.log(8 / 5)
                          + (3 / 5) / (1 - (8 / 5) * a)
                          - float(contraction) / WIDTH
                          + .11 * (8 / 5) / (1 - (8 / 5) * a)
                          + support)
    return {"alpha_range": ["0", str(alpha_max)],
            "collatz_residual_upper_over_alpha": [str(v) for v in
                                                   (residual_zero, residual_fresh,
                                                    residual_uniform)],
            "worst_outer_weight_endpoint": str(worst_x),
            "sparse_exponent_per_alpha_binary64_upper": exponent_per_alpha,
            "note": "Collatz inequalities exact rational; entropy/log bound binary64."}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha", type=float, nargs="+",
                        default=[1e-4, 1e-3, .01, .1, .5, 1.0])
    parser.add_argument("--x", type=float, nargs="+",
                        default=[.104, .16, .25, .4, .5, .6, .75, .84, .896])
    parser.add_argument("--delta", type=float, default=.11)
    parser.add_argument("--summary", action="store_true")
    parser.add_argument("--fixed-one-two", action="store_true")
    parser.add_argument("--sparse-gate", action="store_true")
    args = parser.parse_args()
    if args.sparse_gate:
        print(json.dumps(sparse_gate()))
        return
    if args.fixed_one_two:
        print(json.dumps(fixed_one_two(args.delta)))
        return
    rows = []
    for alpha in args.alpha:
        for x in args.x:
            row = optimize(alpha, x, args.delta)
            rows.append(row)
            if not args.summary:
                print(json.dumps(row), flush=True)
    if args.summary:
        print(json.dumps({"points": len(rows),
                          "nonnegative": sum(r["exponent"] >= 0 for r in rows),
                          "worst": max(rows, key=lambda r: r["exponent"])}))


if __name__ == "__main__":
    main()
