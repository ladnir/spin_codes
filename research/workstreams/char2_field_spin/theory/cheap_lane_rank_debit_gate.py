"""Bounded rank-debit proof gate for five through eight actual groups.

The default mode is only a FLOAT proposal: matrix entries and distant
paths can underflow. --verify-five instead recomputes one fixed rational
r5 witness with directed Arb arithmetic, without an optimizer. Neither
mode establishes a whole-code certificate. All output goes to stdout.

Every zero-to-zero packet marks its actual feedback-column rank
c=sum_b min(4,ell_b). No claim that the allowed inputs are all zero is
made when c<ell. Pair slots and fair orientations are counted exactly by
integer coefficients; only the numerical matrix arithmetic is floating.

The source-rank debit in one outer column is d. The largest32 of64 d's
give a valid source-rank lower bound. For threshold v=0..4r-1, retain
the nonzero-witness cutoff e<=M=128r-32v-1 with a positive mark x<=1:
sum_(e<=M) w_e tau^-e <= x^-M sum_e w_e (x/tau)^e.
All v classes are summed. Their marks may differ but cannot depend on
the sampled setup. Frozen helpers are neither imported nor modified.
"""

import argparse
from collections import defaultdict
from functools import lru_cache
import json
import math
import os
from pathlib import Path
import shutil
import time
from fractions import Fraction

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

NSTATE, NGROUP, NOUT = 17, 8192, 1 << 21
SUPPORT_SUM = np.arange(NSTATE)[:, None]+np.arange(NSTATE)[None, :]
DEADLINE = None


def guard():
    if DEADLINE is not None and time.monotonic() > DEADLINE:
        raise RuntimeError("bounded 50-second invocation budget exhausted")


@lru_cache(maxsize=81)
def occupancy_counts(singles, full):
    """Integer counts of pair slots AND all 2^singles fair orientations.

    Halfpacket pair slots join blocks0/1 or2/3 at the same RS point.
    Positions inside each block are distinct, so rank depends only on
    these four occupancies, not the particular selected points.
    """
    counts = defaultdict(int)
    for s0 in range(singles+1):
        s1 = singles-s0
        for f0 in range(full+1):
            f1 = full-f0
            if s0+f0 > 16 or s1+f1 > 16:
                continue
            slots = (math.comb(16, s0)*math.comb(16-s0, f0)
                     * math.comb(16, s1)*math.comb(16-s1, f1))
            for u0 in range(s0+1):
                for u1 in range(s1+1):
                    occ = tuple(sorted((f0+u0, f0+s0-u0,
                                        f1+u1, f1+s1-u1), reverse=True))
                    counts[occ] += slots*math.comb(s0, u0)*math.comb(s1, u1)
    expected = math.comb(32, singles)*math.comb(32-singles, full)*2**singles
    assert sum(counts.values()) == expected
    return tuple(sorted(counts.items()))


def packet_matrix(occ, tau, z):
    """Nonquiet transitions only; 00 is inserted with rank marker later."""
    ell, occupied = sum(occ), sum(x != 0 for x in occ)
    rank = sum(min(4, x) for x in occ)
    matrix = np.zeros((NSTATE, NSTATE))
    for h in range(NSTATE):
        blocks = (h+3)//4
        emitted = max(0, 13*blocks-sum(occ[:blocks]))
        for hh in range(NSTATE):
            if h == hh == 0:
                continue
            multiplicity = math.comb(16, hh)
            if h == 0:
                if ell and hh >= max(1, 17-4*occupied):
                    activation = min(1.0, ((256**rank-1)//255)*255.0**(hh-16))
                    matrix[h, hh] = multiplicity*z*activation
            elif hh == 0:
                if ell and h <= 4*occupied:
                    matrix[h, hh] = tau**(-h)*z**emitted
            elif h+hh >= 17-4*occupied:
                matrix[h, hh] = (multiplicity*tau**(-min(16-hh, h))
                                 * z**emitted)
    return rank, matrix


def add_matrix(poly, key, matrix):
    if key in poly:
        poly[key] += matrix
    else:
        poly[key] = matrix.copy()


def product(left, right, singles, full):
    guard()
    result = {}
    for (a, b, d), mat1 in left.items():
        for (aa, bb, dd), mat2 in right.items():
            if a+aa <= singles and b+bb <= full:
                add_matrix(result, (a+aa, b+bb, d+dd), mat1@mat2)
    if sum(mat.nbytes for mat in result.values()) > 64*(1 << 20):
        raise RuntimeError("bounded matrix storage exceeds64 MiB")
    return result


def regional(singles, full, tau, z, local_rank=False):
    base = {}
    for s in range(singles+1):
        for f in range(full+1):
            for occ, count in occupancy_counts(s, f):
                rank, nonquiet = packet_matrix(occ, tau, z)
                if local_rank:
                    debits = np.maximum(0, rank-SUPPORT_SUM)
                    for d in range(rank+1):
                        marked = np.where(debits == d, count*nonquiet, 0)
                        if np.any(marked):
                            add_matrix(base, (s, f, d), marked)
                else:
                    add_matrix(base, (s, f, 0), count*nonquiet)
                quiet = np.zeros((NSTATE, NSTATE))
                quiet[0, 0] = count
                add_matrix(base, (s, f, rank), quiet)
    value = base
    for _ in range(7):
        value = product(value, value, singles, full)
    normalizer = (math.comb(4096, singles)*math.comb(4096-singles, full)
                  * 2**singles)
    result = [np.zeros((NSTATE, NSTATE)) for _ in range(singles+2*full+1)]
    for (s, f, d), matrix in value.items():
        if s == singles and f == full:
            result[d] += matrix/normalizer
    return result


def column_product(left, right):
    result = [np.zeros((NSTATE, NSTATE)) for _ in range(len(left)+len(right)-1)]
    for a, aa in enumerate(left):
        for b, bb in enumerate(right):
            result[a+b] += aa@bb
    return result


def threshold_bound(coeff, groups, tau, z, v, log2x, distance):
    guard()
    limit = 128*groups-32*v
    assert limit > 0
    matrix = np.zeros((NSTATE, NSTATE))
    for d, term in enumerate(coeff):
        matrix += math.exp2((log2x-math.log2(tau))*max(d-v, 0))*term
    state = np.zeros(NSTATE)
    state[0] = 1
    logscale = 0.0
    for _ in range(64):
        state = state@matrix
        scale = max(state)
        if not (scale > 0 and math.isfinite(scale)):
            return math.inf
        state /= scale
        logscale += math.log2(scale)
    return (math.log2(math.comb(NGROUP, groups))+limit*math.log2(tau)
            -math.log2(tau-1)-math.floor(distance*NOUT)*math.log2(z)
            -(limit-1)*log2x+logscale+math.log2(sum(state)))


def last_threshold_bound(coeff, groups, tau, z, distance):
    """Exact32-counter cutoff for v=4r-1; FLOAT diagnostic only."""
    guard()
    zero = sum(coeff[:-1])
    one = coeff[-1]/tau
    state = np.zeros((32, NSTATE))
    state[0, 0] = 1
    logscale = 0.0
    for _ in range(64):
        following = state@zero
        following[1:, :] += (state@one)[:-1, :]
        state = following
        scale = np.max(state)
        if not (scale > 0 and math.isfinite(scale)):
            return math.inf
        state /= scale
        logscale += math.log2(scale)
    state *= np.array([1-tau**(-(32-e)) for e in range(32)])[:, None]
    return (math.log2(math.comb(NGROUP, groups))+32*math.log2(tau)
            -math.log2(tau-1)-math.floor(distance*NOUT)*math.log2(z)
            +logscale+math.log2(np.sum(state)))


def logsum(values):
    values = list(values)
    peak = max(values)
    return peak+math.log2(sum(math.exp2(x-peak) for x in values))


def proposal(groups, tau, zs, log2xs, distance, rank_scale=None, local_rank=False,
             last_z=None, last_counter=False):
    profiles = [(a, b) for a in range(groups//2+1)
                for b in range(groups//2+1) if a+b < groups]
    best = {profile: [None]*(4*groups) for profile in profiles}
    geometry = []
    for f in range(groups//2+1):
        s = groups-2*f
        total = sum(w for _, w in occupancy_counts(s, f))
        deficient = sum(w for occ, w in occupancy_counts(s, f)
                        if sum(min(4, x) for x in occ) < groups)
        geometry.append({"packet_singles": s, "packet_full": f,
                         "rank_deficient_count": deficient,
                         "total_oriented_slot_count": total})
    plans = ([(1-rank_scale*(128*groups-32*v), [v]) for v in range(4*groups)]
             if rank_scale is not None else [(z, range(4*groups)) for z in zs])
    if last_z:
        plans = [(z, [v for v in vs if v != 4*groups-1]) for z, vs in plans]
        plans += [(z, [4*groups-1]) for z in last_z]
    for z, selected_thresholds in plans:
        if not selected_thresholds:
            continue
        guard()
        regions = {f: regional(groups-2*f, f, tau, z, local_rank)
                   for f in range(groups//2+1)}
        for profile in profiles:
            even, odd = [regions[f] for f in profile]
            half = column_product(even, odd)
            coeff = column_product(half, half)
            for v in selected_thresholds:
                if last_counter and v == 4*groups-1:
                    bound = last_threshold_bound(coeff, groups, tau, z, distance)
                    old = best[profile][v]
                    if old is None or bound < old["log2_bound"]:
                        best[profile][v] = {"v": v, "z": z, "exact_cutoff": True,
                                            "log2_x": None, "log2_bound": bound}
                    continue
                for lx in log2xs:
                    bound = threshold_bound(coeff, groups, tau, z, v, lx, distance)
                    old = best[profile][v]
                    if old is None or bound < old["log2_bound"]:
                        best[profile][v] = {"v": v, "z": z,
                                            "log2_x": lx, "log2_bound": bound}
    reports = [{"even_full": a, "odd_full": b,
                "log2_sum_bound": logsum(row["log2_bound"] for row in best[(a, b)]),
                "thresholds": best[(a, b)]} for a, b in profiles]
    return {"status": "floating proposal only; underflow can underestimate",
            "groups": groups, "tau": tau, "distance": distance,
            "rank_scale": rank_scale, "all_transition_rank_debit": local_rank,
            "last_exact_counter": last_counter,
            "max_log2_union_bound": max(p["log2_sum_bound"] for p in reports),
            "profiles": reports, "single_packet_geometry": geometry,
            "scope": "exactly r active groups, ideal setup; NOT a certificate",
            "whole_code_certificate": False}


def verify_five(precision, target_bits):
    """Optimizer-free directed replay of one fixed r5 witness."""
    from flint import arb, arb_mat, ctx
    ctx.prec = precision
    groups, tau, n = 5, arb(65280), NSTATE
    distance, limit_weight = Fraction(7, 20), (7*NOUT)//20
    x = 1/arb(2).sqrt()
    profiles = [(a, b) for a in range(3) for b in range(3)]
    terms = {p: [] for p in profiles}

    def add(poly, key, matrix):
        if key in poly:
            poly[key] += matrix
        else:
            poly[key] = arb_mat(matrix)

    def multiply(left, right, singles, full):
        guard()
        out = {}
        for (a, b, d), aa in left.items():
            for (c, e, f), bb in right.items():
                if a+c <= singles and b+e <= full:
                    add(out, (a+c, b+e, d+f), aa*bb)
        if len(out)*n*n*128 > 64*(1 << 20):
            raise RuntimeError("directed matrix storage exceeds bounded64 MiB estimate")
        return out

    def regional_arb(singles, full, z):
        base = {}
        for s in range(singles+1):
            for f in range(full+1):
                for occ, count in occupancy_counts(s, f):
                    ell = sum(occ)
                    occupied = sum(y != 0 for y in occ)
                    rank = sum(min(4, y) for y in occ)
                    mats = {}
                    for h in range(n):
                        blocks = (h+3)//4
                        emitted = max(0, 13*blocks-sum(occ[:blocks]))
                        for hh in range(n):
                            value = arb(0)
                            if h == hh == 0:
                                value = arb(count)
                            elif h == 0:
                                if ell and hh >= max(1, 17-4*occupied):
                                    cap = arb((256**rank-1)//255)*arb(255)**(hh-16)
                                    activation = cap if cap < 1 else arb(1)
                                    value = count*math.comb(16, hh)*z*activation
                            elif hh == 0:
                                if ell and h <= 4*occupied:
                                    value = count*tau**(-h)*z**emitted
                            elif h+hh >= 17-4*occupied:
                                value = (count*math.comb(16, hh)
                                         *tau**(-min(16-hh, h))*z**emitted)
                            if value != 0:
                                debit = max(0, rank-h-hh)
                                if debit not in mats:
                                    mats[debit] = arb_mat(n, n)
                                mats[debit][h, hh] += value
                    for debit, mat in mats.items():
                        add(base, (s, f, debit), mat)
        value = base
        for _ in range(7):
            value = multiply(value, value, singles, full)
        normalizer = (math.comb(4096, singles)*math.comb(4096-singles, full)
                      *2**singles)
        result = [arb_mat(n, n) for _ in range(singles+2*full+1)]
        for (s, f, debit), mat in value.items():
            if s == singles and f == full:
                result[debit] += mat/normalizer
        return result

    def columns(left, right):
        result = [arb_mat(n, n) for _ in range(len(left)+len(right)-1)]
        for a, aa in enumerate(left):
            for b, bb in enumerate(right):
                result[a+b] += aa*bb
        return result

    for v in range(20):
        guard()
        zr = (Fraction(9993, 10000) if v == 19 else
              1-Fraction(27, 2000000)*(640-32*v))
        z = arb(zr.numerator)/zr.denominator
        regions = {f: regional_arb(5-2*f, f, z) for f in range(3)}
        for profile in profiles:
            half = columns(regions[profile[0]], regions[profile[1]])
            coeff = columns(half, half)
            limit = 640-32*v
            prefix = (math.comb(8192, 5)*tau**limit/(tau-1)*z**(-limit_weight))
            if v == 19:
                zero = sum(coeff[:-1], arb_mat(n, n))
                one = coeff[-1]/tau
                state = arb_mat(32, n)
                state[0, 0] = 1
                for _ in range(64):
                    following = state*zero
                    contribution = state*one
                    for e in range(1, 32):
                        for j in range(n):
                            following[e, j] += contribution[e-1, j]
                    state = following
                moment = sum((1-tau**(-(32-e)))*state[e, j]
                             for e in range(32) for j in range(n))
            else:
                matrix = sum(((x/tau)**max(d-v, 0)*mat
                              for d, mat in enumerate(coeff)), arb_mat(n, n))
                state = arb_mat(1, n)
                state[0, 0] = 1
                for _ in range(64):
                    state = state*matrix
                moment = x**(-(limit-1))*sum(state[0, j] for j in range(n))
            terms[profile].append(prefix*moment)
    reports, passed = [], True
    for (a, b), values in terms.items():
        bound = sum(values)
        ok = bool(bound.is_finite() and bound < arb(2)**(-target_bits))
        passed = passed and ok
        reports.append({"even_full": a, "odd_full": b, "passed": ok,
                        "log2_union_bound_interval": str(bound.log()/arb(2).log()),
                        "threshold_log2_intervals": [str(t.log()/arb(2).log()) for t in values]})
    return {"status": "directed fixed-witness replay", "groups": 5,
            "tau": "65280", "distance": str(distance), "precision": precision,
            "cutoff_mark": "1/sqrt(2) for v0..18; exact32-counter for v19",
            "z_rule": "1-(27/2000000)*(640-32v) for v0..18;9993/10000 for v19",
            "target_bits": target_bits, "passed": passed, "profiles": reports,
            "scope": "all nonzero GF128 messages in exactly5 of8192 source groups; ideal setup",
            "whole_code_certificate": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=5)
    parser.add_argument("--tau", type=float, default=65280)
    parser.add_argument("--z", type=float, nargs="+", default=[0.995])
    parser.add_argument("--log2-x", type=float, nargs="+", default=[0, -4, -16])
    parser.add_argument("--distance", type=float, default=0.35)
    parser.add_argument("--rank-scale", type=float,
                        help="use z_v=1-scale*(128*r-32*v), instead of --z")
    parser.add_argument("--local-rank", action="store_true",
                        help="mark rank max(0,c-h-h') on every transition")
    parser.add_argument("--details", action="store_true")
    parser.add_argument("--last-z", type=float, nargs="+",
                        help="override candidate marks only for the final rank threshold")
    parser.add_argument("--last-counter", action="store_true")
    parser.add_argument("--verify-five", action="store_true",
                        help="replay fixed r5 witness; proposal tuning flags are ignored")
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--target-bits", type=int, default=40)
    args = parser.parse_args()
    if not (5 <= args.groups <= 8 and 1 < args.tau <= 65280
            and 0 < args.distance < 0.5
            and all(0 < z < 1 for z in args.z)
            and all(-64 <= x <= 0 for x in args.log2_x)
            and (args.rank_scale is None or 0 < args.rank_scale < 1/(128*args.groups))
            and (args.last_z is None or (1 <= len(args.last_z) <= 4
                                        and all(0 < z < 1 for z in args.last_z)))
            and 1 <= len(args.z) <= 4 and 1 <= len(args.log2_x) <= 8):
        parser.error("require r5..8, tau1..65280, valid distance, at most4 z and8 x marks")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        parser.error("at least20 GiB free is required before numerical runs")
    global DEADLINE
    started = time.monotonic()
    DEADLINE = started+50
    if args.verify_five:
        if not (args.groups == 5 and 64 <= args.precision <= 512
                and 1 <= args.target_bits <= 10000):
            parser.error("verify-five requires group5, precision64..512, target1..10000")
        result = verify_five(args.precision, args.target_bits)
    else:
        result = proposal(args.groups, args.tau, args.z, args.log2_x, args.distance,
                          args.rank_scale, args.local_rank, args.last_z, args.last_counter)
    result["elapsed_seconds"] = time.monotonic()-started
    if not args.details:
        for profile in result["profiles"]:
            if args.verify_five:
                del profile["threshold_log2_intervals"]
            else:
                profile["worst_threshold"] = max(profile["thresholds"], key=lambda row: row["log2_bound"])
                del profile["thresholds"]
    print(json.dumps(result, sort_keys=True), flush=True)
    if args.verify_five and not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
