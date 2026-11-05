"""Bounded gates for joint activation/rank charging.

The default FLOAT screen fixes an actual source-group set and its pair
profile. Every local deterministic input rank is charged at sqrt(tau),
because the largest32 column ranks dominate half their total. At an
activation, rank and the activation event are averaged JOINTLY:
E[s^-d 1[d<c]] <= P_256(c) s^-(c-1) 255^-(16-h').
Here c=rank(C_I), d=rank(R_J restricted to im(C_I)), s=sqrt(tau),
and P_256(c)=(256^c-1)/255. This is not a product of marginal bounds.

In that default mode, only histories with a nonzero INCOMING state
at some emitted packet
are counted. A positive two-flag transfer excludes all-quiet histories;
an unused final outgoing state is not counted as physical activity.
Quiet histories need a separate proof. The rank potential is bounded
by tau^(128*r-rank)/(tau-1), without a full-rank cutoff in this screen.

--rank-distribution instead retains a separate positive coefficient for
each possible activation rank. It restores the frozen top32 threshold
operator and its nonzero-source cutoff, without multiplying marginal
activation and rank bounds. That mode includes quiet histories.

--verify-six recomputes a fixed r6 witness with directed Arb arithmetic,
without an optimizer. It includes quiet histories, activation/live
rank classes, all24 top32 threshold terms, and the nonzero-source cutoff.
Proposal tuning flags do not change this fixed witness.

The other modes use floating arithmetic: underflow can underestimate,
so their outputs are NOT certificates. Runs are time/memory bounded,
stdout only, and require20 GiB free. The frozen geometry helper is hash
checked. No mode establishes a whole-code certificate.
"""

import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from functools import lru_cache
from fractions import Fraction

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

NSTATE, NGROUP, NOUT = 17, 8192, 1 << 21
DEADLINE = None
FROZEN_HASH = "04b62fddce34dcdf60f14d9969573fbf0818a4755eeffb593845e37a955e9053"


def guard():
    if DEADLINE is not None and time.monotonic() > DEADLINE:
        raise RuntimeError("bounded50-second invocation exhausted")


def load_frozen():
    path = Path(__file__).with_name("cheap_lane_rank_debit_gate.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != FROZEN_HASH:
        raise RuntimeError("frozen geometry helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_rank_debit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.DEADLINE = DEADLINE
    return module


def local_transfer(occ, tau, z):
    """U marks no emitted nonzero state; F marks one seen already."""
    c = sum(min(4, x) for x in occ)
    occupied = sum(x != 0 for x in occ)
    s = math.sqrt(tau)
    unflagged = np.zeros((NSTATE, NSTATE))
    flagged = np.zeros((NSTATE, NSTATE))
    unflagged[0, 0] = s**(-c)
    if c:
        gamma = ((256**c-1)//255)/s**(c-1)
        for hh in range(max(1, 17-4*occupied), NSTATE):
            joint = min(1.0, gamma*255.0**(hh-16))
            unflagged[0, hh] = math.comb(16, hh)*z*joint
    for h in range(1, NSTATE):
        blocks = (h+3)//4
        emitted = max(0, 13*blocks-sum(occ[:blocks]))
        for hh in range(NSTATE):
            allowed = ((hh == 0 and c and h <= 4*occupied)
                       or (hh > 0 and h+hh >= 17-4*occupied))
            if allowed:
                debit = max(0, c-h-hh)
                flagged[h, hh] = (math.comb(16, hh)*z**emitted
                                   * tau**(-min(16-hh, h))*s**(-debit))
    return unflagged, flagged


def compose(left, right):
    u, f = left
    v, g = right
    return u@v, f@(v+g)+u@g


def add(poly, key, pair):
    if key in poly:
        poly[key][0][:] += pair[0]
        poly[key][1][:] += pair[1]
    else:
        poly[key] = pair[0].copy(), pair[1].copy()


def product(left, right, singles, full):
    guard()
    result = {}
    for (a, b), pair in left.items():
        for (c, d), other in right.items():
            if a+c <= singles and b+d <= full:
                add(result, (a+c, b+d), compose(pair, other))
    if len(result)*2*NSTATE*NSTATE*8 > 64*(1 << 20):
        raise RuntimeError("matrix storage exceeds64 MiB")
    return result


def regional(singles, full, tau, z, geometry):
    base = {}
    for a in range(singles+1):
        for b in range(full+1):
            for occ, count in geometry(a, b):
                u, f = local_transfer(occ, tau, z)
                add(base, (a, b), (count*u, count*f))
    value = base
    for _ in range(7):
        value = product(value, value, singles, full)
    norm = math.comb(4096, singles)*math.comb(4096-singles, full)*2**singles
    u, f = value[(singles, full)]
    return u/norm, f/norm


def profile_bound(regions, profile, groups, tau, z, distance):
    half = compose(regions[profile[0]], regions[profile[1]])
    u, f = compose(half, half)
    unflagged, flagged = np.zeros(NSTATE), np.zeros(NSTATE)
    unflagged[0] = 1
    logscale = 0.0
    for _ in range(64):
        guard()
        next_u = unflagged@u
        next_f = flagged@(u+f)+unflagged@f
        scale = max(np.max(next_u), np.max(next_f))
        if not (scale > 0 and math.isfinite(scale)):
            raise RuntimeError("floating underflow or overflow in profile recurrence")
        unflagged, flagged = next_u/scale, next_f/scale
        logscale += math.log2(scale)
    nonquiet = sum(flagged)
    if not nonquiet > 0:
        raise RuntimeError("floating nonquiet moment underflow")
    prefix = (math.log2(math.comb(NGROUP, groups))+128*groups*math.log2(tau)
              -math.log2(tau-1)-math.floor(distance*NOUT)*math.log2(z))
    return {"even_full": profile[0], "odd_full": profile[1], "z": z,
            "log2_nonquiet_union_bound": prefix+logscale+math.log2(nonquiet),
            "unflagged_terminal_scaled_mass": float(sum(unflagged)),
            "flagged_terminal_scaled_mass": float(nonquiet)}


@lru_cache(maxsize=2048)
def activation_rank_floor(occ, checks):
    """Rank on image directions supported inside <=checks coordinates.

    A rank-c_b subspace in four coordinates has intersection dimension
    at least max(0,c_b-4+q_b) with any q_b-coordinate subspace. The MDS
    check matrix is injective on a union of at most checks coordinates.
    """
    value = [0]+[-100]*(checks)
    for ell in occ:
        c = min(4, ell)
        next_value = [-100]*(checks+1)
        for used, rank in enumerate(value):
            if rank < 0:
                continue
            for q in range(min(4, checks-used)+1):
                next_value[used+q] = max(next_value[used+q],
                                         rank+max(0, c-4+q))
        value = next_value
    return max(value)


def rank_regional(singles, full, tau, z, frozen, live_rank=False):
    """Positive polynomial tagged by actual/guaranteed input rank.

    At h=0, tag every feasible actual rank d with its own event cap.
    At h>0, use the existing pointwise lower rank c-h-h'.
    """
    base = {}
    for a in range(singles+1):
        for b in range(full+1):
            for occ, count in frozen.occupancy_counts(a, b):
                c, matrix = frozen.packet_matrix(occ, tau, z)
                matrix[0, :] = 0
                if live_rank:
                    tagged = {}
                    for h in range(1, NSTATE):
                        for j in range(NSTATE):
                            if not matrix[h, j]:
                                continue
                            m, k = 16-j, max(0, c-h)
                            if j == 0:
                                masses = [(k, 1.0)]
                            elif h >= m or not k:
                                masses = [(0, 1.0)]
                            else:
                                n, checks = 16-h, m-h
                                masses = []
                                for d in range(max(0, checks+k-n), min(checks, k)+1):
                                    cap = (1.0 if d == k else
                                           min(1.0, ((256**k-1)//255)*255.0**(-checks)
                                               /((256**(k-d)-1)//255)))
                                    masses.append((d, cap))
                            for d, cap in masses:
                                if d not in tagged:
                                    tagged[d] = np.zeros((NSTATE, NSTATE))
                                tagged[d][h, j] += count*matrix[h, j]*cap
                    for d, marked in tagged.items():
                        frozen.add_matrix(base, (a, b, d), marked)
                else:
                    debits = np.maximum(0, c-frozen.SUPPORT_SUM)
                    for d in range(c+1):
                        marked = np.where(debits == d, count*matrix, 0)
                        if np.any(marked):
                            frozen.add_matrix(base, (a, b, d), marked)
                quiet = np.zeros((NSTATE, NSTATE))
                quiet[0, 0] = count
                frozen.add_matrix(base, (a, b, c), quiet)
                if not c:
                    continue
                occupied = sum(x != 0 for x in occ)
                for j in range(max(1, 17-4*occupied), NSTATE):
                    m = 16-j
                    lower = max(0, c-j, activation_rank_floor(occ, m))
                    for d in range(lower, min(m, c-1)+1):
                        cap = min(1.0, ((256**c-1)//255)*255.0**(-m)
                                  /((256**(c-d)-1)//255))
                        marked = np.zeros((NSTATE, NSTATE))
                        marked[0, j] = count*math.comb(16, j)*z*cap
                        frozen.add_matrix(base, (a, b, d), marked)
    value = base
    for _ in range(7):
        value = frozen.product(value, value, singles, full)
    norm = math.comb(4096, singles)*math.comb(4096-singles, full)*2**singles
    result = [np.zeros((NSTATE, NSTATE)) for _ in range(singles+2*full+1)]
    for (a, b, d), matrix in value.items():
        if a == singles and b == full:
            result[d] += matrix/norm
    return result


def rank_proposal(args, profiles, frozen):
    used = sorted({f for profile in profiles for f in profile})
    best = {p: [None]*(4*args.groups) for p in profiles}
    plans = [(1-args.rank_scale*(128*args.groups-32*v), v)
             for v in range(4*args.groups-1)]
    plans += [(z, 4*args.groups-1) for z in args.last_z]
    for z, v in plans:
        guard()
        regions = {f: rank_regional(args.groups-2*f, f, args.tau, z, frozen,
                                    args.live_rank)
                   for f in used}
        for p in profiles:
            half = frozen.column_product(regions[p[0]], regions[p[1]])
            coeff = frozen.column_product(half, half)
            if v == 4*args.groups-1:
                value = frozen.last_threshold_bound(coeff, args.groups, args.tau,
                                                    z, args.distance)
                candidates = [(value, None)]
            else:
                candidates = [(frozen.threshold_bound(coeff, args.groups,
                               args.tau, z, v, lx, args.distance), lx)
                              for lx in args.log2_x]
            value, lx = min(candidates)
            if best[p][v] is None or value < best[p][v]["log2_bound"]:
                best[p][v] = {"v": v, "z": z, "log2_x": lx,
                              "log2_bound": value}
    reports = []
    for p in profiles:
        thresholds = best[p]
        reports.append({"even_full": p[0], "odd_full": p[1],
                        "log2_union_bound": frozen.logsum(t["log2_bound"]
                                                          for t in thresholds),
                        "worst_threshold": max(thresholds,
                                               key=lambda t: t["log2_bound"]),
                        "thresholds": thresholds})
    return reports


def verify_six(frozen, precision, target_bits):
    """Directed, optimizer-free replay of the fixed r6 rank witness.

    Uses the activation knapsack floor and the live quotient-rank caps.
    All thresholds include the full source-support union. The maximum
    over matching profiles needs no extra profile multiplicity.
    """
    from flint import arb, arb_mat, ctx
    ctx.prec = precision
    tau, n, weight_limit = arb(65280), NSTATE, (7*NOUT)//20
    profiles = [(a, b) for a in range(4) for b in range(4) if a+b < 6]
    values = {p: [] for p in profiles}

    def put(poly, key, mat):
        if key in poly:
            poly[key] += mat
        else:
            poly[key] = arb_mat(mat)

    def multiply(left, right, singles, full):
        guard()
        out = {}
        for (a, b, d), aa in left.items():
            for (c, e, f), bb in right.items():
                if a+c <= singles and b+e <= full:
                    put(out, (a+c, b+e, d+f), aa*bb)
        if len(out)*n*n*128 > 64*(1 << 20):
            raise RuntimeError("directed matrix estimate exceeds64 MiB")
        return out

    def regional_arb(singles, full, z):
        base = {}
        for a in range(singles+1):
            for b in range(full+1):
                for occ, count in frozen.occupancy_counts(a, b):
                    c = sum(min(4, ell) for ell in occ)
                    occupied = sum(ell != 0 for ell in occ)
                    mats = {c: arb_mat(n, n)}
                    mats[c][0, 0] = count
                    for h in range(n):
                        blocks = (h+3)//4
                        emitted = max(0, 13*blocks-sum(occ[:blocks]))
                        for j in range(n):
                            m = 16-j
                            if h == 0:
                                if not (c and j >= max(1, 17-4*occupied)):
                                    continue
                                lo = max(0, c-j, activation_rank_floor(occ, m))
                                ranks = []
                                for d in range(lo, min(m, c-1)+1):
                                    cap = (arb((256**c-1)//255)*arb(255)**(-m)
                                           /((256**(c-d)-1)//255))
                                    ranks.append((d, cap if cap < 1 else arb(1)))
                                coefficient = count*math.comb(16, j)*z
                            else:
                                allowed = ((j == 0 and c and h <= 4*occupied)
                                           or (j > 0 and h+j >= 17-4*occupied))
                                if not allowed:
                                    continue
                                k = max(0, c-h)
                                if j == 0:
                                    ranks = [(k, arb(1))]
                                elif h >= m or not k:
                                    ranks = [(0, arb(1))]
                                else:
                                    checks, slots = m-h, 16-h
                                    ranks = []
                                    for d in range(max(0, checks+k-slots), min(checks, k)+1):
                                        cap = (arb(1) if d == k else
                                               arb((256**k-1)//255)*arb(255)**(-checks)
                                               /((256**(k-d)-1)//255))
                                        ranks.append((d, cap if cap < 1 else arb(1)))
                                coefficient = (count*math.comb(16, j)*z**emitted
                                               *tau**(-min(m, h)))
                            for d, cap in ranks:
                                if d not in mats:
                                    mats[d] = arb_mat(n, n)
                                mats[d][h, j] += coefficient*cap
                    for d, mat in mats.items():
                        put(base, (a, b, d), mat)
        value = base
        for _ in range(7):
            value = multiply(value, value, singles, full)
        norm = math.comb(4096, singles)*math.comb(4096-singles, full)*2**singles
        result = [arb_mat(n, n) for _ in range(singles+2*full+1)]
        for (a, b, d), mat in value.items():
            if a == singles and b == full:
                result[d] += mat/norm
        return result

    def columns(left, right):
        out = [arb_mat(n, n) for _ in range(len(left)+len(right)-1)]
        for a, aa in enumerate(left):
            for b, bb in enumerate(right):
                out[a+b] += aa*bb
        return out

    for v in range(24):
        guard()
        zr = (Fraction(3981, 4000) if v == 12 else
              Fraction(3997, 4000) if v == 23 else
              1-Fraction(27, 2000000)*(768-32*v))
        z = arb(zr.numerator)/zr.denominator
        x = arb(1) if v == 12 else 1/arb(2).sqrt()
        regions = {f: regional_arb(6-2*f, f, z) for f in range(4)}
        for p in profiles:
            half = columns(regions[p[0]], regions[p[1]])
            coeff = columns(half, half)
            limit = 768-32*v
            prefix = math.comb(8192, 6)*tau**limit/(tau-1)*z**(-weight_limit)
            if v == 23:
                zero = sum(coeff[:-1], arb_mat(n, n))
                one = coeff[-1]/tau
                state = arb_mat(32, n)
                state[0, 0] = 1
                for _ in range(64):
                    following, incoming = state*zero, state*one
                    for e in range(1, 32):
                        for j in range(n):
                            following[e, j] += incoming[e-1, j]
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
            values[p].append(prefix*moment)
    reports, passed = [], True
    for (a, b), terms in values.items():
        bound = sum(terms)
        ok = bool(bound.is_finite() and bound < arb(2)**(-target_bits))
        passed = passed and ok
        reports.append({"even_full": a, "odd_full": b, "passed": ok,
                        "log2_union_bound_interval": str(bound.log()/arb(2).log()),
                        "threshold_log2_intervals": [str(t.log()/arb(2).log()) for t in terms]})
    return {"status": "directed fixed-witness replay", "groups": 6,
            "tau": "65280", "distance": "7/20", "precision": precision,
            "target_bits": target_bits, "passed": passed,
            "activation_rank_floor": "four-block coordinate-intersection knapsack",
            "live_rank_caps": "quotient image dimension k0=max(0,c-h)",
            "z_rule": "v12:3981/4000; v23:3997/4000; others:1-(27/2000000)(768-32v)",
            "cutoff_rule": "v12:x1; v23:exact32-counter; others:x1/sqrt2",
            "profiles": reports, "all_quiet_histories": "included",
            "scope": "all nonzero GF128 messages supported on exactly6 of8192 groups; ideal setup",
            "whole_code_certificate": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=6)
    parser.add_argument("--tau", type=float, default=65280)
    parser.add_argument("--distance", type=float, default=.35)
    parser.add_argument("--z", type=float, nargs="+", default=[.99, .995, .998])
    parser.add_argument("--profile", type=int, nargs=2, default=[2, 3])
    parser.add_argument("--all-profiles", action="store_true")
    parser.add_argument("--rank-distribution", action="store_true")
    parser.add_argument("--live-rank", action="store_true",
                        help="tag quotient-rank caps on live transitions too")
    parser.add_argument("--rank-scale", type=float, default=.0000135)
    parser.add_argument("--log2-x", type=float, nargs="+", default=[-.5])
    parser.add_argument("--last-z", type=float, nargs="+", default=[.9993, .9995, .9997])
    parser.add_argument("--details", action="store_true")
    parser.add_argument("--verify-six", action="store_true",
                        help="replay fixed r6 witness; proposal tuning flags are ignored")
    parser.add_argument("--precision", type=int, default=192)
    parser.add_argument("--target-bits", type=int, default=40)
    args = parser.parse_args()
    if not (5 <= args.groups <= 8 and 1 < args.tau <= 65280
            and 0 < args.distance < .5 and 1 <= len(args.z) <= 4
            and all(0 < z < 1 for z in args.z)
            and 0 < args.rank_scale < 1/(128*args.groups)
            and 1 <= len(args.last_z) <= 4 and all(0 < z < 1 for z in args.last_z)
            and 1 <= len(args.log2_x) <= 4 and all(-64 <= x <= 0 for x in args.log2_x)
            and all(0 <= f <= args.groups//2 for f in args.profile)
            and sum(args.profile) < args.groups):
        parser.error("require r5..8, tau1..65280, valid profile/distance, at most4 z")
    free = shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free < 20*(1 << 30):
        parser.error("at least20 GiB free required")
    global DEADLINE
    started = time.monotonic()
    DEADLINE = started+50
    frozen = load_frozen()
    if args.verify_six:
        if not (args.groups == 6 and 64 <= args.precision <= 512
                and 1 <= args.target_bits <= 10000):
            parser.error("verify-six requires r6, precision64..512, positive bounded target")
        result = verify_six(frozen, args.precision, args.target_bits)
        if not args.details:
            for row in result["profiles"]:
                del row["threshold_log2_intervals"]
        result["elapsed_seconds"] = time.monotonic()-started
        result["free_bytes_at_start"] = free
        print(json.dumps(result, sort_keys=True), flush=True)
        if not result["passed"]:
            raise SystemExit(1)
        return
    profiles = ([(a, b) for a in range(args.groups//2+1)
                 for b in range(args.groups//2+1) if a+b < args.groups]
                if args.all_profiles else [tuple(args.profile)])
    used = sorted({f for profile in profiles for f in profile})
    if args.rank_distribution:
        reports = rank_proposal(args, profiles, frozen)
        if not args.details:
            for report in reports:
                del report["thresholds"]
    else:
        reports = []
        for z in args.z:
            regions = {f: regional(args.groups-2*f, f, args.tau, z,
                                   frozen.occupancy_counts) for f in used}
            reports.extend(profile_bound(regions, p, args.groups, args.tau, z,
                                         args.distance) for p in profiles)
    print(json.dumps({"status": "FLOAT proposal only; underflow can underestimate",
                      "groups": args.groups, "tau": args.tau,
                      "distance": args.distance, "reports": reports,
                      "rank_distribution": args.rank_distribution,
                      "live_rank_distribution": args.live_rank,
                      "all_quiet_histories": ("included" if args.rank_distribution
                                               else "excluded; separate proof required"),
                      "scope": "exactly r groups; ideal setup; mode determines quiet inclusion",
                      "whole_code_certificate": False,
                      "elapsed_seconds": time.monotonic()-started,
                      "free_bytes_at_start": free}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
