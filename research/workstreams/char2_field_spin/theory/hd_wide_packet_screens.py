"""Bounded full-Cauchy screens for wider unpunctured packets.

The rank cover is the existing high_distance_fields formula, with its
validation range extended to width/state <=64 and width+state <=128.
Only a finite tau grid is searched; exact rank roots are cached. Positive
sampled margins are not a whole-code or directed-rounding certificate.
"""

import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
import math

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from cauchy_firstmoment_screen import geometry, log_choose, log_terminal_power
from gf8_rank_history_screen import modal_log_penalty, source_row_log
from high_distance_fields import local_matrices as original_local, roots
from high_distance_window_obstruction import screen as window_screen
from higher_distance_screen import output_threshold
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients, log_terminal


# b, k, width, state, sparse, exact target, tau, fixed marks, displayed margin.
WITNESSES = (
    (128, 64, 32, 16, True, 16, 64, 6.297228538498054, 894.1465302557804),
    (128, 64, 32, 16, False, 62915, 4,
     (-2.8193095149790306, -.0780151804712123, 2.0603043646564583), -73514.16890069892),
    (128, 64, 32, 16, False, 838861, 16,
     (1.097680704597706, -1.6085040857910273, -1.267112492938427), -625141.364455111),
    (128, 64, 32, 32, True, 16, 256, 6.5696898200319085, 1360.3281759759593),
    (128, 64, 32, 32, False, 62915, 4,
     (-2.8169494899332848, -.07817846173502761, 2.0609052787469473), -73536.67059680686),
    (128, 64, 32, 32, False, 838861, 16,
     (1.0986173814229536, -1.609442753872432, -1.2671288116839812), -625140.4605299972),
    (64, 32, 32, 64, True, 16, 256, 5.876635235305124, 541.2459404582668),
    (64, 32, 32, 64, False, 62915, 4,
     (-2.8140117494382704, -.07841883724433317, 2.034433104771449), -76941.97954229922),
    (64, 32, 32, 64, False, 838861, 16,
     (1.0657356528059827, -1.575745938312668, -1.2694839791088663), -642707.2406489819),
    (128, 64, 64, 32, True, 16, 256, 6.569689828241869, 1360.4486069129364),
    (128, 64, 64, 32, False, 62915, 4,
     (-2.5426857420287434, -.10202960289116121, 2.0509221708551695), -78207.55473987752),
    (128, 64, 64, 32, False, 838861, 16,
     (1.098885172623033, -1.6097102369061826, -1.2671548078211015), -625139.5164783401),
    (128, 64, 64, 32, True, 16, 65535, 7.25817371978083, 2925.795620779394),
    (128, 64, 64, 32, False, 62915, 4096,
     (-2.5124336429405276, -.0779723306243227, -.7259922359728143), -297822.4175860373),
    (128, 64, 64, 32, False, 838861, 4096,
     (4.103278649999029, -4.134563447121628, -3.772772855900015), -3070769.3108791984),
)


@lru_cache(maxsize=128)
def local_constants(tau, width, state, bits):
    q = 1 << bits
    gamma = np.asarray(roots(F(tau), width, state, q), dtype=float)
    powers = np.asarray([float(F(1, tau) ** r)
                         for r in range(state + width + 1)])
    flats = np.cumsum([math.comb(state, k) * powers[k] for k in range(state)])
    return gamma, powers, flats


def local_matrices(y, tau, width=32, state=16, bits=16):
    if not (1 <= width <= 64 and 1 <= state <= 64 and width + state <= 128
            and 0 < y <= 1 and 1 <= tau <= (1 << bits) - 1):
        raise ValueError("Invalid full-Cauchy geometry or marks")
    gamma, powers, flats = local_constants(tau, width, state, bits)
    result = np.zeros((width + 1, 2, 2))
    result[0, 0, 0] = 1
    terms = [math.comb(width, z) * y ** (width - z) * (1 - y) ** z
             for z in range(width + 1)]
    for w in range(width + 1):
        for z in range(width + 1):
            pivots, remaining = min(w, z), max(0, z - w)
            if remaining < state:
                result[w, 1, 1] += terms[z] * gamma[remaining] * powers[pivots]
            if w and z < w:
                result[w, 0, 1] += terms[z] * powers[z]
            if w > state + z:
                result[w, 0, 0] += terms[z] * powers[state + z]
        if w:
            minimum = max(0, width + 1 - w)
            result[w, 0, 1] = min(result[w, 0, 1], y ** minimum)
            pure = y ** minimum * powers[state] * flats[min(w - 1, state - 1)]
            joint = sum(terms[z] * powers[state + z] * flats[min(w - z - 1, state - 1)]
                        for z in range(min(w, width + 1)))
            result[w, 1, 0] = min(pure, joint)
    return result


def bound(point, target, tau, sparse=True, length=128, dimension=64,
          width=32, state=16, bits=16, k_log=20, distance=F(2, 5)):
    g = geometry(k_log, length, dimension, math.lcm(64, width))
    threshold = output_threshold(g["n"], distance)
    logtau = math.log(tau)
    if sparse:
        active = int(target)
        logy = -math.exp(float(point)) * active / g["n"]
        local = local_matrices(math.exp(logy), tau, width, state, bits)[:min(width, active) + 1]
        with np.errstate(divide="ignore"):
            factor = np.log(local)
        for w in range(len(factor)):
            factor[w] += math.log(math.comb(width, w))
        result = np.array([[[0., -math.inf], [-math.inf, 0.]]])
        steps = g["rows"] // width
        while steps:
            if steps & 1:
                result = log_polynomial_product(result, factor, active)
            steps >>= 1
            if steps:
                factor = log_polynomial_product(factor, factor, active)
        regional = regional_from_coefficients(result, active, g["rows"], logtau)
        return (log_choose(g["rows"], active)
                - ((length - dimension) * active + 1) * logtau
                - threshold * logy + log_terminal(regional, length)) / math.log(2)
    logx, logy, logtilt = map(float, point)
    local = local_matrices(math.exp(logy), tau, width, state, bits)
    packet = np.einsum("w,wij->ij", [math.comb(width, w) * math.exp(w * logx)
                                     for w in range(width + 1)], local)
    return (length * modal_log_penalty(g["rows"])
            + g["rows"] * source_row_log(target / g["n"], logtau, logtilt,
                                         length, dimension, True)
            + log_terminal_power(packet, g["n"] // width)
            - target * (logx + logtilt) - threshold * logy - logtau) / math.log(2)


def search(value, sparse=True, taus=(4, 16, 64, 256), iterations=25, **options):
    g = geometry(options.get("k_log", 20), options["length"], options["dimension"],
                 math.lcm(64, options["width"]))
    target = int(value) if sparse else round(float(value) * g["n"])
    candidates = []
    previous = None
    for tau in taus:
        objective = lambda point: bound(point, target, tau, sparse=sparse, **options)
        if sparse:
            result = minimize_scalar(objective, bounds=(-2, 10), method="bounded",
                                     options={"maxiter": iterations, "xatol": 1e-7})
        else:
            starts = [(-2., -.2, 0.), (0., -1., -1.)]
            if previous is not None:
                starts.insert(0, previous)
            limits = [(-20, min(20, 640 / options["width"])), (-20, 0), (-20, 20)]
            choices = [minimize(objective, point, method="L-BFGS-B", bounds=limits,
                                options={"maxiter": iterations, "ftol": 1e-10, "maxls": 12})
                       for point in starts]
            result = min(choices, key=lambda item: item.fun if math.isfinite(item.fun) else math.inf)
            previous = result.x.tolist()
        candidates.append((result, tau))
    result, tau = min(candidates, key=lambda item: item[0].fun if math.isfinite(item[0].fun) else math.inf)
    point = float(result.x) if sparse else result.x.tolist()
    return dict(**options, mode="active" if sparse else "weight", target=target,
                tau=tau, margin_bits=-float(result.fun), fixed_tau_parameters=point,
                tau_grid=list(taus), optimizer_success=bool(result.success),
                rigorous=False, whole_code=False,
                diagonal_distribution="independent uniform nonzero full-field scalars")


def checks():
    cases = 0
    for width, state in ((4, 4), (8, 8), (16, 8), (32, 16), (32, 32)):
        for tau in (4, 16, 256):
            for y in (.2, .9, .9999):
                new = local_matrices(y, tau, width, state)
                old = np.asarray(original_local(y, float(tau), width, width, state, 16))
                if not np.allclose(new, old, rtol=2e-13, atol=1e-280):
                    raise AssertionError((width, state, tau, y, np.max(abs(new - old))))
                cases += 1
    print(json.dumps(dict(local_operator_regressions=cases, passed=True)))


def replay():
    for b, k, w, s, sparse, target, tau, marks, expected in WITNESSES:
        margin = -bound(marks, target, tau, sparse, b, k, w, s)
        if not math.isclose(margin, expected, rel_tol=2e-10, abs_tol=2e-5):
            raise AssertionError((b, k, w, s, target, margin, expected))
        print(json.dumps(dict(outer=[b, k], width=w, state=s, mode="active" if sparse else "weight",
                              target=target, tau=tau, margin_bits=margin,
                              rigorous=False, whole_code=False)), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--window", action="store_true")
    parser.add_argument("--outer", type=int, nargs=2, default=(128, 64))
    parser.add_argument("--width", type=int, default=32)
    parser.add_argument("--state", type=int, default=16)
    parser.add_argument("--bits", type=int, default=16)
    parser.add_argument("--taus", type=int, nargs="+", default=(4, 16, 64, 256))
    parser.add_argument("--active", type=int, nargs="*", default=[])
    parser.add_argument("--density", type=float, nargs="*", default=[])
    parser.add_argument("--iterations", type=int, default=25)
    args = parser.parse_args()
    if args.check:
        checks()
    if args.replay:
        replay()
    if args.window:
        row = window_screen(*args.outer, args.state, args.width)
        # The complete exact Fraction is available from the imported gate;
        # avoid printing thousands of denominator digits in routine screens.
        row.pop("failure_probability_lower_bound")
        print(json.dumps(row), flush=True)
    options = dict(length=args.outer[0], dimension=args.outer[1], width=args.width,
                   state=args.state, bits=args.bits)
    for sparse, values in ((True, args.active), (False, args.density)):
        for value in values:
            print(json.dumps(search(value, sparse, args.taus, args.iterations, **options)), flush=True)


if __name__ == "__main__":
    main()
