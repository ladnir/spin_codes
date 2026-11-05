"""FLOAT tile diagnostic: pool masks and route placements before alpha.

A tile contains p consecutive packets and 32p paired route slots. Let c_d
be its raw degree-d coefficient, including the final mask count D_j=C(16,j).
For each positive activity class separately, use the Jensen upper bound

    C(32p,d) D_j [c_d / (C(32p,d) D_j)]^alpha.

The terminal tile has no outgoing mask, so its normalization omits D_j.
Raw coefficients sum all internal masks before the fractional power.
The binomial factor counts every route placement, not merely placements
with nonzero scalar cap. Quiet and nonquiet classes remain separate.
Their union is positive; no subtraction or FFT is used.

The frozen grouped-mask helper supplies alpha1 packet caps, including
deterministic gamma minima, joint-live filtering, and shared-point
singleton extinction. This wrapper does not alter those caps. It first
forms raw packet polynomials, doubles to p packets, applies Jensen once,
then doubles to a 128-packet region. The terminal product is P^(p-1)T.
Degree truncation at128 is safe because every polynomial degree is
nonnegative. The full artificial chain retains32767 regular updates.

Each invocation tests one tile size2,4,8,or16 at r256,z.85,alpha.523.
Repeating the fully paired regional marginal in both shifted matchings is
not a realizable proper source support. Whole-chain scores are artificial
comparisons, NOT support-class bounds or distance certificates.
Use python -B --self-test or --tile 2. Guards:20GiB free,48seconds,
under256MiB estimated peak memory. Stdout only; no generated files.
"""

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
GROUPED_SHA256 = "4526522c47c9920482f8990bb104243eada12318774ee2d0dbb68fda3704ab36"


def load_grouped():
    path = Path(__file__).with_name("cheap_lane_grouped_masks.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != GROUPED_SHA256:
        raise RuntimeError("frozen grouped-mask helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_grouped_masks", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def raw_polynomials(m, t, degrees, weights, np, middle):
    degree = int(degrees.max())
    rows = m.shape[1]
    q = np.zeros((degree+1, rows)); n = np.zeros((degree+1, rows, rows))
    tq = np.zeros((degree+1, 1)); tn = np.zeros((degree+1, rows, 1))
    np.add.at(q, degrees, weights[:, None]*m[:, 0])
    np.add.at(n[:, 1:], degrees, weights[:, None, None]*m[:, 1:])
    np.add.at(tq, degrees, weights[:, None]*t[:, 0])
    np.add.at(tn[:, 1:], degrees, weights[:, None, None]*t[:, 1:])
    return (middle.to_log(q), middle.to_log(n)), (middle.to_log(tq), middle.to_log(tn))


def jensen(poly, slots, multiplicities, alpha, np):
    """Positive log coefficients; multiplicities=[1] for a terminal."""
    q, n = poly
    counts = np.array([math.log(math.comb(slots, d)) for d in range(len(q))])
    ends = np.log(np.asarray(multiplicities, dtype=float))
    norm = counts[:, None]+ends[None, :]
    return alpha*q+(1-alpha)*norm, alpha*n+(1-alpha)*norm[:, None, :]


def double_pair(a, b, times, limit, paired):
    for _ in range(times):
        paired.guard()
        b, a = (paired.positive_convolution(a, b, limit),
                paired.positive_convolution(a, a, limit))
    return a, b


def toy_test(grouped, paired):
    """Exact binomial route counts, endpoint counts, and alpha1 invariance."""
    np = paired.np
    counts = (1, 2, 1)
    p, t, matrices, terminals = [], [], [], []
    for degree in range(3):
        cap = [[Fraction(1+i+degree, 7+2*i+3*j+degree)
                for j in range(3)] for i in range(3)]
        end = [[Fraction(2+i+degree, 11+2*i+degree)] for i in range(3)]
        m = [[cap[i][j]*counts[j] for j in range(3)] for i in range(3)]
        matrices.append(m); terminals.append(end)
        coefficient = paired.exact_zero(3, 3)
        final = paired.exact_zero(3, 1)
        weight = math.comb(2, degree)
        for i in range(3):
            coefficient[int(i > 0)][i] = [weight*x for x in m[i]]
            final[int(i > 0)][i] = [weight*x for x in end[i]]
        p.append(coefficient); t.append(final)
    exact_p, exact_t = p, t
    logs_p, logs_t = paired.exact_to_logs(p), paired.exact_to_logs(t)
    checkpoints = {}
    comparisons = 0
    for tile in (1, 2, 4):
        if tile > 1:
            exact_t, exact_p = (paired.exact_poly_product(exact_p, exact_t, 8),
                                paired.exact_poly_product(exact_p, exact_p, 8))
            logs_p, logs_t = double_pair(logs_p, logs_t, 1, 8, paired)
        jp = jensen(logs_p, 2*tile, counts, 1.0, np)
        jt = jensen(logs_t, 2*tile, [1], 1.0, np)
        comparisons += paired.compare_logs(exact_p, jp)+paired.compare_logs(exact_t, jt)
        checkpoints[tile] = double_pair(jp, jt, 2-(tile.bit_length()-1), 8, paired)
        # Exact route normalization is C(2*tile,d), including zero terms.
        for d in range(2*tile+1):
            route_count = sum(math.prod(math.comb(2, x) for x in split)
                              for split in itertools.product(range(3), repeat=tile)
                              if sum(split) == d)
            assert route_count == math.comb(2*tile, d)
    for result in checkpoints.values():
        for expected, actual in zip(checkpoints[1], result):
            for e, a in zip(expected, actual):
                assert np.array_equal(np.isfinite(e), np.isfinite(a))
                finite = np.isfinite(e)
                assert np.allclose(e[finite], a[finite], atol=3e-12, rtol=3e-12)
    # For two packets, compare explicit occupancy powers with Jensen after
    # pooling placements of the same total degree. The latter must be larger.
    alpha = .523
    exact_route_q = np.zeros((5, 3)); exact_route_n = np.zeros((5, 3, 3))
    exact_route_tq = np.zeros((5, 1)); exact_route_tn = np.zeros((5, 3, 1))
    for a, b in itertools.product(range(3), repeat=2):
        ma = np.array(matrices[a], dtype=float)
        mb = np.array([matrices[b]], dtype=float)
        tb = np.array([terminals[b]], dtype=float)
        q, n, tq, tn = grouped.pair_batch(ma, mb, tb, alpha, np.array(counts), np)
        weight = math.comb(2, a)*math.comb(2, b)
        exact_route_q[a+b] += weight*q[0]
        exact_route_n[a+b] += weight*n[0]
        exact_route_tq[a+b] += weight*tq[0]
        exact_route_tn[a+b] += weight*tn[0]
    raw_pair = paired.exact_poly_product(p, p, 4)
    raw_tail = paired.exact_poly_product(p, t, 4)
    upper = (*jensen(paired.exact_to_logs(raw_pair), 4, counts, alpha, np),
             *jensen(paired.exact_to_logs(raw_tail), 4, [1], alpha, np))
    inequalities = 0
    for exact, bound in zip((exact_route_q, exact_route_n, exact_route_tq, exact_route_tn), upper):
        assert np.all(exact <= np.exp(bound)*(1+5e-13)+1e-200)
        inequalities += exact.size
    return {"alpha1_log_comparisons": comparisons, "tile_partition_invariance": [1, 2, 4],
            "exact_binomial_route_normalization": "passed",
            "explicit_pair_vs_degree_Jensen_inequalities": inequalities,
            "terminal_and_positive_activity_split": "passed"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--self-test", action="store_true")
    modes.add_argument("--tile", type=int, choices=(2, 4, 8, 16))
    options = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on actual source volume")
    start = time.monotonic()
    grouped = load_grouped(); singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(start+48)
    tests = {**toy_test(grouped, paired), **grouped.toy_test(paired)}
    paired.guard()
    if options.self_test:
        print(json.dumps({"tests": tests, "seconds": time.monotonic()-start, "generated_bytes": 0}, indent=2))
        return
    args = SimpleNamespace(r=256, z=.85, alpha=.523, activity_mark=1.0)
    filtered, middle = paired.imports()
    m, t, degrees, weights = grouped.raw_tables(args, singleton, paired, filtered)
    a, b = raw_polynomials(m, t, degrees, weights, paired.np, middle)
    tile_log2 = options.tile.bit_length()-1
    a, b = double_pair(a, b, tile_log2, 128, paired)
    a = jensen(a, 32*options.tile, [math.comb(16, j) for j in range(17)], args.alpha, paired.np)
    b = jensen(b, 32*options.tile, [1], args.alpha, paired.np)
    coefficient_bytes = sum(x.nbytes for x in (*a, *b))
    # The frozen engine's largest degree-batch temporary is 129*17^3 doubles.
    estimated_peak_bytes = 64*(1 << 20)
    assert coefficient_bytes < estimated_peak_bytes < 256*(1 << 20)
    result = grouped.evaluate(args, a, b, 7-tile_log2, paired, filtered, middle)
    print(json.dumps({"scope": "FLOAT artificial repeated-marginal degree-Jensen tile diagnostic",
        "whole_code_certificate": False, "actual_support_class_bound": False,
        "parameters": vars(args), "tile_packets": options.tile,
        "normalization": "C(32*tile,degree)*C(16,destination); terminal omits destination count",
        "raw_gamma_minimum": [0, .5, 1], "singleton_point_average_alpha": 1,
        "result": result, "tests": tests, "tile_coefficient_bytes": coefficient_bytes,
        "estimated_peak_bytes": estimated_peak_bytes, "seconds": time.monotonic()-start,
        "generated_bytes": 0, "grouped_source_sha256": GROUPED_SHA256}, indent=2))


if __name__ == "__main__":
    main()
