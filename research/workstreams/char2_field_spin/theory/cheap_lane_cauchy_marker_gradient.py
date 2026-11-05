"""Directed first derivatives of the fixed r512 Cauchy anchor, without tuning.

Keep the frozen z=7/10, alpha=119/1280 and markers u=1/225,v=2/15.
For the actual repeated matching profile (f,s)=(16,480), the evaluated
chain numerator is H=e0 P(u,v)^32767 T(u,v). Its normalized denominator
is [M_4096(f,s)u^f v^s]^256. Source and source-count factors are constant
in u,v. Therefore the derivative of log2(bound) with respect to ln(u)
is [d_ln(u)H/H-256f]/ln(2), and similarly for v and256s.

The ratios are the total full-pair/singleton counts under the positive
tilted complete-chain measure. Their regional means should equal f,s
at marker stationarity. Positive dual products compute them directly;
there is no finite difference, Perron approximation, or optimization.
Caps are raised to alpha before route averaging, exactly as in the
frozen anchor. The terminal contributes route counts but no update.

Scope is only the selected actual profile, not a full size certificate.
Source-only dependencies are hash pinned.30-second cooperative deadline,
128MiB Windows process-memory guard,20GiB disk guard, stdout only.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
ANCHOR_PIN = "319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"


def load_anchor():
    path = Path(__file__).with_name("cheap_lane_cone_cauchy_anchors.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != ANCHOR_PIN:
        raise RuntimeError("frozen Cauchy-anchor hash mismatch")
    spec = importlib.util.spec_from_file_location("marker_gradient_anchor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dual_product(x, y):
    return (x[0]*y[0], x[1]*y[0]+x[0]*y[1], x[2]*y[0]+x[0]*y[2])


def dual_matmul(left, right, zero):
    result = []
    for row in left:
        out = []
        for j in range(len(right[0])):
            value = [zero, zero, zero]
            for k, item in enumerate(row):
                term = dual_product(item, right[k][j])
                for d in range(3):
                    value[d] += term[d]
            out.append(tuple(value))
        result.append(out)
    return result


def dual_matpow(value, n, zero, one):
    result = [[(one if i == j else zero, zero, zero) for j in range(2)]
              for i in range(2)]
    while n:
        if n & 1:
            result = dual_matmul(result, value, zero)
        n //= 2
        if n:
            value = dual_matmul(value, value, zero)
    return result


def shape_moments(u, v, guard):
    """Exact value, full-count moment and singleton-count moment per shape."""
    half = defaultdict(lambda: [F(0), F(0), F(0)])
    for f in range(17):
        for s in range(17-f):
            mass = F(comb(16, f)*comb(16-f, s), 2**s)*u**f*v**s
            for left in range(s+1):
                weight = comb(s, left)*mass
                item = half[f+left, f+s-left]
                item[0] += weight; item[1] += f*weight; item[2] += s*weight
    unordered = [(a, b, tuple(x*(1 if a == b else 2) for x in value))
                 for (a, b), value in half.items() if a <= b]
    result = defaultdict(lambda: [F(0), F(0), F(0)])
    for i, (a, b, value) in enumerate(unordered):
        guard()
        for j in range(i, len(unordered)):
            c, d, other = unordered[j]
            factor = 1 if i == j else 2
            term = dual_product(value, other)
            target = result[tuple(sorted((a, b, c, d)))]
            for axis in range(3):
                target[axis] += factor*term[axis]
    sums = [sum((value[d] for value in result.values()), F(0)) for d in range(3)]
    assert sums == [(1+u+v)**32, 32*u*(1+u+v)**31, 32*v*(1+u+v)**31]
    return dict(result)


def exact_tests():
    """Noncommuting two-state polynomials versus exact dual evaluation."""
    p = [[{(0, 0): F(1, 3), (1, 0): F(2, 5)}, {(0, 1): F(3, 7)}],
         [{(1, 0): F(5, 11)}, {(0, 0): F(7, 13), (0, 1): F(11, 17)}]]
    t = [[{(1, 0): F(13, 19)}], [{(0, 1): F(17, 23), (0, 0): F(19, 29)}]]
    u, v = F(2, 7), F(3, 11)

    def poly_matmul(left, right):
        out = []
        for row in left:
            rr = []
            for j in range(len(right[0])):
                value = defaultdict(F)
                for k, item in enumerate(row):
                    for (f, s), x in item.items():
                        for (ff, ss), y in right[k][j].items():
                            value[f+ff, s+ss] += x*y
                rr.append(dict(value))
            out.append(rr)
        return out

    def evaluate(poly):
        value = [F(0), F(0), F(0)]
        for (f, s), coefficient in poly.items():
            mass = coefficient*u**f*v**s
            value[0] += mass; value[1] += f*mass; value[2] += s*mass
        return tuple(value)

    pd = [[evaluate(x) for x in row] for row in p]
    td = [[evaluate(x) for x in row] for row in t]
    expected = poly_matmul(poly_matmul(poly_matmul(p, p), p), t)
    actual = dual_matmul(dual_matpow(pd, 3, F(0), F(1)), td, F(0))
    for i in range(2):
        assert actual[i][0] == evaluate(expected[i][0])
    return 6


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+30
    anchor = load_anchor()
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("30-second gradient budget exceeded; no automatic retry")
        if peak_memory() > anchor.MEMORY_LIMIT:
            raise MemoryError("128MiB process-memory budget exceeded")

    checks = exact_tests()
    adapter = anchor.load_adapter()
    core = adapter.load_verifier()
    from flint import arb, ctx
    ctx.prec = args.precision

    class GradientVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = GradientVerifier(2, deadline)
    r, f, s, u, v, z = anchor.ANCHORS[0]
    core.Z, core.ALPHA = z, anchor.ALPHA
    moments = shape_moments(u, v, guard)
    frozen_weights = anchor.shape_weights(u, v, guard)
    assert {key: value[0] for key, value in moments.items()} == frozen_weights
    del frozen_weights
    shapes = list(combinations_with_replacement(range(17), 4))
    assert set(moments) == set(shapes)
    caps = verifier.local_caps(shapes)
    exponent = verifier.exact(anchor.ALPHA)
    entries = [[arb(0), arb(0), arb(0)] for _ in range(6)]
    for occ, weights in moments.items():
        guard()
        masses = [verifier.exact(x) for x in weights]
        for j, cap in enumerate(caps[occ]):
            powered = verifier.positive(cap)**exponent
            for axis in range(3):
                entries[j][axis] += masses[axis]*powered
    p = [[tuple(entries[0]), tuple(entries[1])], [tuple(entries[2]), tuple(entries[3])]]
    t = [[tuple(entries[4])], [tuple(entries[5])]]
    value = dual_matmul(dual_matpow(p, 32767, arb(0), arb(1)), t, arb(0))[0][0]
    total_means = [value[axis]/value[0] for axis in (1, 2)]
    targets = [256*f, 256*s]
    logtwo = arb(2).log()
    log_gradients = [(mean-target)/logtwo for mean, target in zip(total_means, targets)]
    denominator = (verifier.exact(u)**f*verifier.exact(v)**s
                   *comb(4096, f)*comb(4096-f, s))
    h, c, eta = anchor.source_factor(r)
    source = (verifier.exact(eta)*verifier.exact(core.BETA)**r
              *verifier.exact(z)**(-core.LIMIT))**exponent
    number = anchor.run_total(4096, r, s)
    bound = number*source*value[0]/denominator**256
    # Reproduce the frozen regional assembly, rather than just an equivalent
    # scalar exponent, as an independent normalization/terminal check.
    plain_p = [[x[0] for x in row] for row in p]
    plain_t = [[x[0] for x in row] for row in t]
    rr = [[x/denominator for x in row] for row in verifier.matpow(plain_p, 128)]
    tt = [[x/denominator for x in row]
          for row in verifier.matmul(verifier.matpow(plain_p, 127), plain_t)]
    regional_bound = number*source*verifier.matmul(verifier.matpow(rr, 255), tt)[0][0]
    assert bound.overlaps(regional_bound)
    guard()
    print(json.dumps({"scope": "marker-gradient diagnostic for one actual r512 profile",
        "directed": True, "whole_code_certificate": False, "optimizer_run": False,
        "parameters": {"source_weight": r, "ordered_profile": [[f, s], [f, s]],
                       "z": str(z), "alpha": str(anchor.ALPHA),
                       "full_marker": str(u), "singleton_marker": str(v)},
        "log2_bound_interval": str(bound.log()/logtwo),
        "derivative_coordinate": "ln(marker); derivative of log2(bound)",
        "full": {"tilted_total_mean": str(total_means[0]), "total_target": targets[0],
                 "tilted_mean_per_region": str(total_means[0]/256), "regional_target": f,
                 "tilted_mean_per_packet": str(total_means[0]/32768),
                 "gradient_interval": str(log_gradients[0])},
        "singleton": {"tilted_total_mean": str(total_means[1]), "total_target": targets[1],
                      "tilted_mean_per_region": str(total_means[1]/256), "regional_target": s,
                      "tilted_mean_per_packet": str(total_means[1]/32768),
                      "gradient_interval": str(log_gradients[1])},
        "interpretation": "negative gradient favors increasing that marker locally; no finite-step gain claimed",
        "tests": {"exact_dual_polynomial_equalities": checks,
                  "shape_weights_equal_frozen_anchor": True,
                  "exact_unweighted_packet_moments": True,
                  "regional_and_global_chains_overlap": True, "updates": 32767},
        "precision_bits": args.precision, "seconds": time.monotonic()-started,
        "peak_process_bytes": peak_memory(), "generated_bytes": 0,
        "anchor_pin": ANCHOR_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
