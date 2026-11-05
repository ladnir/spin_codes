"""Exact-degree pair-only route diagnostic; FLOAT, not a certificate.

This helper removes coefficient-marker slack for one regional marginal.
Let r be even and f=r/2. The region has f selected full pairs among4096
pair slots and no singleton pairs. A packet's two half-packets contain
a,b full pairs, so its four block occupancies are (a,a,b,b). Define

 P_d = sum_{a+b=d} C(16,a)C(16,b) M_(a,a,b,b),
 T_d = sum_{a+b=d} C(16,a)C(16,b) E_(a,a,b,b).

M is the frozen filtered scalar matrix majorant. E is its output-only
terminal vector. Both retain the positive incoming-state activity flag.
The regional averages are [t^f]P(t)^128/C(4096,f) and
[t^f]P(t)^127 T(t)/C(4096,f). Swapping a full pair changes no occupancy,
so there is no extra orientation factor. The coefficient counts pairs,
not records. It exactly averages a local majorant, not encoder failures.

Seven doublings (A,B)<-(A*A,A*B), starting(P,T), compute both expressions.
Positive coefficient convolution is truncated above f; discarded degrees
cannot contribute to any retained coefficient. All coefficients are stored
in logarithmic form. There is no FFT, cancellation, global scale, or
subtraction of quiet from total. The quiet flag has only its incoming0
row, which is represented separately for faster positive multiplication.

The scalar local kernel, fractional moment, actual support filters,
source prefactor, and emission-only final packet are unchanged. Frozen
source imports are checked by SHA256. Small Fraction tests compare
direct slot enumeration, sequential multiplication, doubling, activity
flags, and terminal products for1,2 and4 packets with1 or2 slots per half.

--gamma-mix allows the deterministic local choice gamma in{0,1/2,1}.
For each fixed occupancy and incoming/outgoing state size, the smallest
scalar table entry specifies a deterministic gamma. Fix that choice
before exposing D1 or any requested-zero subset.
All choices split the same global z=a*b, and use the same source bound,
alpha, destination multiplicity, and activity marker. The terminal may
make its own fixed choice. This is not a minimum of incompatible rank-tag
CDF measures after a rank counter; it instantiates the permitted local
floor/zero-expansion choice before that packet's proof relaxation.

--joint-live-filter retains actual nonzero next-state feasibility jointly
with the scalar rank cost. For a fixed requested-zero system, let G be
the E-linear image of V=Cx+u on its kernel. A live outgoing state requires
rank(R_J|G)<dim(G), where |J|=16-j. The projective D1 estimate gives
E[q^-rank(R_J|G) 1(rank<dim G)] <= F_j,
F_j=min(1,(q/(q-1))*255^-(16-j)). Consequently a regular live-destination
entry is also at most F_j times its output-only scalar kernel. This
filter is imposed before introducing the auxiliary rank count; it is
not an assumption that every auxiliary V is nonzero. The implementation
uses the terminal output-only table, with F_j^alpha and the destination
multiplicity C(16,j). Activity factors remain outside alpha; j=0 is
unchanged. Each gamma receives this cap before the optional gamma minimum.

Regional marginal scores are valid diagnostics for the specified f.
Any whole-chain score repeats this same marginal through both alternating
matchings only to compare numerical relaxations. For0<r<8192, no actual
source support is fully paired in both matchings. Consequently that
whole-chain score is artificial and is NOT a class bound, worst-case
bound, or distance certificate, even if it is negative.

No output data or caches are written. Every invocation checks20GiB free,
uses under256MiB, and has a48-second deadline. Arithmetic is FLOAT, not
directed; exact degree extraction does not make floating sums certified.
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

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

FILTERED_SHA256 = "4198a26f54c3074fc953bb3f752c236de1fb0ab0e2d9c9c4b74bb10bd9405d78"
DEADLINE = 0.0
LN2 = math.log(2)
INF = float("-inf")


def guard():
    if time.monotonic() > DEADLINE:
        raise RuntimeError("bounded48-second invocation exhausted")


def imports():
    path = Path(__file__).with_name("cheap_lane_filtered_scalar.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != FILTERED_SHA256:
        raise RuntimeError("frozen filtered helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_filtered_scalar", path)
    filtered = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(filtered)
    filtered.DEADLINE = DEADLINE
    middle = filtered.reference()
    middle.DEADLINE = DEADLINE
    return filtered, middle


def positive_convolution(left, right, limit):
    """Log coefficients (quiet row, nonquiet matrix), composed by flag OR.

    Quiet matrices have only row0. Thus the nonquiet product is the
    positive sum N_left*N_right + N_left*Q_right + Q_left*N_right.
    The latter two products use an outer product and a single row.
    """
    aq, an = left
    bq, bn = right
    degree = min(limit, len(aq)+len(bq)-2)
    rows, columns = an.shape[1], bn.shape[2]
    assert an.shape[2] == bn.shape[1] == aq.shape[1]
    cq = np.full((degree+1, columns), INF)
    cn = np.full((degree+1, rows, columns), INF)
    for d in range(degree+1):
        guard()
        lo, hi = max(0, d-len(bq)+1), min(d, len(aq)-1)
        ids = np.arange(lo, hi+1)
        js = d-ids
        cq[d] = np.logaddexp.reduce(aq[ids, 0, None]+bq[js], axis=0)
        dense = an[ids, :, :, None]+bn[js, None, :, :]
        dense = np.logaddexp.reduce(np.logaddexp.reduce(dense, axis=2), axis=0)
        cross = np.logaddexp.reduce(an[ids, :, 0, None]+bq[js, None, :], axis=0)
        cn[d] = np.logaddexp(dense, cross)
        initial = aq[ids, :, None]+bn[js]
        initial = np.logaddexp.reduce(np.logaddexp.reduce(initial, axis=1), axis=0)
        cn[d, 0] = np.logaddexp(cn[d, 0], initial)
    assert not np.isnan(cq).any() and not np.isnan(cn).any()
    return cq, cn


def exact_flag_product(left, right):
    """Small exact2-flag matrix product, independent of the log routine."""
    flags, rows, middle = len(left), len(left[0]), len(left[0][0])
    columns = len(right[0][0])
    out = [[[Fraction(0) for _ in range(columns)] for _ in range(rows)]
           for _ in range(flags)]
    for a in range(flags):
        for b in range(flags):
            for i in range(rows):
                for k in range(middle):
                    for j in range(columns):
                        out[a | b][i][j] += left[a][i][k]*right[b][k][j]
    return out


def exact_add(left, right):
    return [[[a+b for a, b in zip(lrow, rrow)]
             for lrow, rrow in zip(lflag, rflag)]
            for lflag, rflag in zip(left, right)]


def exact_zero(rows, columns):
    return [[[Fraction(0) for _ in range(columns)] for _ in range(rows)]
            for _ in range(2)]


def exact_poly_product(left, right, limit):
    out = [exact_zero(len(left[0][0]), len(right[0][0][0]))
           for _ in range(min(limit, len(left)+len(right)-2)+1)]
    for a, aa in enumerate(left):
        for b, bb in enumerate(right):
            if a+b < len(out):
                out[a+b] = exact_add(out[a+b], exact_flag_product(aa, bb))
    return out


def exact_to_logs(poly):
    rows, columns = len(poly[0][0]), len(poly[0][0][0])
    q = np.full((len(poly), columns), INF)
    n = np.full((len(poly), rows, columns), INF)
    for d, coefficient in enumerate(poly):
        assert all(x == 0 for row in coefficient[0][1:] for x in row)
        for j in range(columns):
            if coefficient[0][0][j]:
                q[d, j] = math.log(coefficient[0][0][j])
        for i in range(rows):
            for j in range(columns):
                if coefficient[1][i][j]:
                    n[d, i, j] = math.log(coefficient[1][i][j])
    return q, n


def compare_logs(poly, logs):
    exact = exact_to_logs(poly)
    comparisons = 0
    for expected, actual in zip(exact, logs):
        assert np.array_equal(np.isfinite(expected), np.isfinite(actual))
        mask = np.isfinite(expected)
        assert np.allclose(expected[mask], actual[mask], atol=3e-12, rtol=3e-12)
        comparisons += int(mask.sum())
    return comparisons


def toy_tables(a, b):
    """Two states and rational unequal kernels for small occupancy tests."""
    d = a+b
    matrix = exact_zero(2, 2)
    terminal = exact_zero(2, 1)
    matrix[0][0] = [Fraction(1, 2+d), Fraction(1+a, 7+2*d)]
    matrix[1][1] = [Fraction(1+b, 11+d), Fraction(2, 13+d)]
    terminal[0][0][0] = Fraction(1, 3+d)
    terminal[1][1][0] = Fraction(1+a+b, 17+d)
    return matrix, terminal


def gf4_multiply(a, b):
    """Polynomial basis modulo x^2+x+1, only for the exact toy test."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def gf4_rank(rows, columns):
    matrix = [row[:] for row in rows]
    rank = 0
    inverse = (0, 1, 3, 2)
    for column in range(columns):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = inverse[matrix[rank][column]]
        matrix[rank] = [gf4_multiply(x, scale) for x in matrix[rank]]
        for i in range(len(matrix)):
            if i != rank and matrix[i][column]:
                scale = matrix[i][column]
                matrix[i] = [x ^ gf4_multiply(scale, y)
                             for x, y in zip(matrix[i], matrix[rank])]
        rank += 1
    return rank


def joint_live_toy_test():
    """All GF4^2 subspaces and fresh nonzero diagonal choices, exactly."""
    subspaces = [[]]+[[(1, t)] for t in range(4)]+[[(0, 1)], [(1, 0), (0, 1)]]
    h = ((3, 2), (2, 3))
    averages = rank_cases = 0
    for basis in subspaces:
        g = len(basis)
        for mask in range(4):
            rows = [j for j in range(2) if mask >> j & 1]
            total = Fraction(0)
            for d0, d1 in itertools.product((1, 2, 3), repeat=2):
                matrix = [[gf4_multiply(h[j][0], gf4_multiply(d0, v[0]))
                           ^ gf4_multiply(h[j][1], gf4_multiply(d1, v[1]))
                           for v in basis] for j in rows]
                rank = gf4_rank(matrix, g)
                total += Fraction(4)**(-rank) if rank < g else 0
                rank_cases += 1
            mean = total/9
            cap = Fraction(4, 3)*(1-Fraction(4)**(-g))*Fraction(3)**(-len(rows))
            assert mean <= cap
            if not g:
                assert mean == cap == 0
            averages += 1
    return {"field": "GF4 modulo x^2+x+1", "subspaces": len(subspaces),
            "exact_joint_moment_averages": averages, "rank_configurations": rank_cases,
            "status": "passed"}


def self_test():
    checks = configurations = 0
    for slots_per_half in (1, 2):
        p = [exact_zero(2, 2) for _ in range(2*slots_per_half+1)]
        t = [exact_zero(2, 1) for _ in range(2*slots_per_half+1)]
        for a, b in itertools.product(range(slots_per_half+1), repeat=2):
            m, final = toy_tables(a, b)
            multiplicity = math.comb(slots_per_half, a)*math.comb(slots_per_half, b)
            scaled_m = [[[multiplicity*x for x in row] for row in flag] for flag in m]
            scaled_t = [[[multiplicity*x for x in row] for row in flag] for flag in final]
            p[a+b] = exact_add(p[a+b], scaled_m)
            t[a+b] = exact_add(t[a+b], scaled_t)
        for packets in (1, 2, 4):
            guard()
            slot_count = 2*slots_per_half*packets
            limit = min(3, slot_count)
            direct_p = [exact_zero(2, 2) for _ in range(limit+1)]
            direct_t = [exact_zero(2, 1) for _ in range(limit+1)]
            counts = [0]*(limit+1)
            for slots in itertools.product(range(2), repeat=slot_count):
                degree = sum(slots)
                if degree > limit:
                    continue
                matrices = []
                for i in range(packets):
                    start = 2*slots_per_half*i
                    a = sum(slots[start:start+slots_per_half])
                    b = sum(slots[start+slots_per_half:start+2*slots_per_half])
                    matrices.append(toy_tables(a, b))
                regular, final = matrices[0]
                for m, emitted in matrices[1:]:
                    final = exact_flag_product(regular, emitted)
                    regular = exact_flag_product(regular, m)
                direct_p[degree] = exact_add(direct_p[degree], regular)
                direct_t[degree] = exact_add(direct_t[degree], final)
                configurations += 1
                counts[degree] += 1
            a, b = p[:limit+1], t[:limit+1]
            la, lb = exact_to_logs(a), exact_to_logs(b)
            for _ in range(packets.bit_length()-1):
                b, a = exact_poly_product(a, b, limit), exact_poly_product(a, a, limit)
                lb, la = positive_convolution(la, lb, limit), positive_convolution(la, la, limit)
            assert a == direct_p and b == direct_t
            checks += compare_logs(a, la)+compare_logs(b, lb)
            for degree in range(limit+1):
                count = math.comb(slot_count, degree)
                assert count == counts[degree]
                for poly in (a, b):
                    assert all(x/count >= 0 for flag in poly[degree] for row in flag for x in row)
    return {"exact_fraction_configurations": configurations,
            "log_coefficient_comparisons": checks, "packets_checked": [1, 2, 4],
            "slots_per_half_checked": [1, 2],
            "terminal_and_activity_flags": "passed", "slot_normalization": "passed",
            "joint_live_GF4": joint_live_toy_test()}


def packet_polynomials(args, filtered, middle):
    shapes = np.array(sorted(set(tuple(sorted((a, a, b, b)))
                                 for a, b in itertools.product(range(17), repeat=2))), dtype=np.int16)
    index = {tuple(s): i for i, s in enumerate(shapes)}
    local, terminal = filtered.local_tables(shapes, args.z, args.alpha, args.gamma, args.activity_mark)
    original, original_terminal = local.copy(), terminal.copy()
    def joint_cap(matrix, final):
        if args.joint_live_filter:
            for j in range(1, 17):
                factor = min(1.0, 256/255*255.0**(-(16-j)))
                candidate = final[:, :, :, 0]*factor**args.alpha*math.comb(16, j)
                matrix[:, :, :, j] = np.minimum(matrix[:, :, :, j], candidate)
        return matrix
    local = joint_cap(local, terminal)
    gamma_effect = None
    if args.gamma_mix:
        for gamma in (0, .5, 1):
            if gamma != args.gamma:
                next_local, next_terminal = filtered.local_tables(
                    shapes, args.z, args.alpha, gamma, args.activity_mark)
                next_local = joint_cap(next_local, next_terminal)
                local = np.minimum(local, next_local)
                terminal = np.minimum(terminal, next_terminal)
    if args.gamma_mix or args.joint_live_filter:
        changed = (local > 0) & (local < original)
        gamma_effect = {"improved_local_entries": int(changed.sum()),
                        "improved_terminal_entries": int((terminal < original_terminal).sum()),
                        "maximum_local_log2_gain": float(np.log2(original[changed]/local[changed]).max())
                        if changed.any() else 0.0}
    f = args.r//2
    degree = min(32, f)
    p = np.zeros((degree+1, 4, 17, 17))
    t = np.zeros((degree+1, 4, 17, 1))
    for a, b in itertools.product(range(17), repeat=2):
        if a+b <= degree:
            coefficient = math.comb(16, a)*math.comb(16, b)
            row = index[tuple(sorted((a, a, b, b)))]
            p[a+b] += coefficient*local[row]
            t[a+b] += coefficient*terminal[row]
    assert np.count_nonzero(p[:, 0, 1:, :]) == 0
    assert np.count_nonzero(t[:, 0, 1:, :]) == 0
    return ((middle.to_log(p[:, 0, 0, :]), middle.to_log(p[:, 2])),
            (middle.to_log(t[:, 0, 0, :]), middle.to_log(t[:, 2])),
            int(local.nbytes+terminal.nbytes+p.nbytes+t.nbytes), gamma_effect)


def expand_coefficient(poly, degree, denominator):
    quiet, nonquiet = poly
    result = np.full((4, nonquiet.shape[1], nonquiet.shape[2]), INF)
    result[0, 0] = quiet[degree]-denominator
    result[2] = nonquiet[degree]-denominator
    return result


def main():
    global DEADLINE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--r", type=int, default=256)
    parser.add_argument("--z", type=float, default=.85)
    parser.add_argument("--alpha", type=float, default=.523)
    parser.add_argument("--gamma", type=float, choices=(0, .5, 1), default=1)
    parser.add_argument("--gamma-mix", action="store_true")
    parser.add_argument("--joint-live-filter", action="store_true")
    parser.add_argument("--activity-mark", type=float, default=1)
    args = parser.parse_args()
    if not (2 <= args.r <= 256 and args.r % 2 == 0 and .001 <= args.z <= 1
            and .1 <= args.alpha <= 1 and 1e-6 <= args.activity_mark <= 1e6):
        raise ValueError("parameters outside bounded pair-only interface")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on actual source volume")
    start = time.monotonic()
    DEADLINE = start+48
    tests = self_test()
    if args.self_test:
        print(json.dumps({"scope": "small exact coefficient tests only", "tests": tests,
                          "seconds": time.monotonic()-start, "generated_bytes": 0}, indent=2))
        return
    filtered, middle = imports()
    a, b, initial_bytes, gamma_effect = packet_polynomials(args, filtered, middle)
    steps = []
    for step in range(7):
        guard()
        stage = time.monotonic()
        b, a = positive_convolution(a, b, args.r//2), positive_convolution(a, a, args.r//2)
        steps.append({"packets": 2**(step+1), "degree": len(a[0])-1,
                      "seconds": time.monotonic()-stage})
    f = args.r//2
    denominator = math.log(math.comb(4096, f))
    region = expand_coefficient(a, f, denominator)
    tail = expand_coefficient(b, f, denominator)
    result = filtered.bound(args, region, tail, middle)
    print(json.dumps({"scope": "FLOAT exact pair-only REGIONAL coefficient diagnostic",
                      "whole_code_certificate": False, "actual_support_class_bound": False,
                      "repeated_marginal_whole_chain_is_artificial": True,
                      "r": args.r, "f": f, "z": args.z, "alpha": args.alpha,
                      "gamma": args.gamma, "activity_mark": args.activity_mark,
                      "gamma_mix": args.gamma_mix, "gamma_mix_effect": gamma_effect,
                      "joint_live_filter": args.joint_live_filter,
                      "regional_log2": {"quiet00": float(region[0, 0, 0]/LN2),
                                         "nonquiet00": float(region[2, 0, 0]/LN2),
                                         "nonquiet16to16": float(region[2, 16, 16]/LN2),
                                         "terminal_quiet0": float(tail[0, 0, 0]/LN2)},
                      "artificial_repeated_marginal_score": result,
                      "tests": tests, "steps": steps,
                      "initial_array_bytes": initial_bytes,
                      "coefficient_array_bytes": sum(x.nbytes for x in (*a, *b)),
                      "seconds": time.monotonic()-start, "generated_bytes": 0,
                      "filtered_source_sha256": FILTERED_SHA256}, indent=2))


if __name__ == "__main__":
    main()
