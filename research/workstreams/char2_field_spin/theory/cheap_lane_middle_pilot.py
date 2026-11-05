"""Bounded FLOAT intermediate-support screen, NOT a certificate.

This screen keeps actual paired routing, a common output mark, and positive
full/deficient-source-rank and quiet/nonquiet flags. It upgrades only quiet
0->0 and live16->16 transitions by their output-zero rank expansion; other
entries use the previously audited floor-only rank caps. The final packet
has output equations only: it pays no fictitious feedback/refresh factor.

For an outer column let e=4r-d be its source-rank deficit. The largest32
column ranks give remaining dimension at most half the total deficit. A
nonzero source also requires at least33 deficient columns. Thus local
weights tau^((ell-d)/2-k), followed by a column marker y>=1 and y^-33,
bound the source potential by 1/(tau-1). This is a deliberately lossy
cutoff, not the exact best32 rank counter. tau is a proof parameter in
(1,65280], not a change to the construction's refresh distribution.

Every paired-route profile f uses its own positive coefficient markers
(with a zero boundary marker for an exactly absent singleton/full class).
The denominator is the exact trinomial coefficient evaluated at those
markers (computed in FLOAT here). Take separate entrywise regional
envelopes for each of the four positive flags. No rare occupancy tail is
discarded. The all-quiet class is reported separately, never subtracted
from a total. A source union factor binom(8192,r) is paid once.

The full-rank regional envelope is conservatively 1 at entry00 and zero
elsewhere when r>0: under these tags any activation or nonempty live
packet is deficient; an empty packet cannot extinguish. Therefore a
full-rank regional path must stay quiet throughout. This is a support
statement after coefficient extraction, not a claim about the evaluated
polynomial, which contains the empty-input monomial.

All results are diagnostics. Arithmetic is positive but not directed;
the finite-state/marker relaxations may be very loose. No files or caches
are written. Source arrays occupy under100 MiB; runtime is capped50s.
"""

import argparse
from collections import defaultdict
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

DEADLINE = 0.0
LN2 = math.log(2)
INF = float("-inf")


def guard():
    if time.monotonic() > DEADLINE:
        raise RuntimeError("bounded50-second invocation exhausted")


def pq(k):
    return (256**k-1)//255


def activation_floor(c, checks):
    table = [0]+[-99]*16
    for cb in c:
        nxt = [-99]*17
        for used in range(17):
            for q in range(min(4, used)+1):
                nxt[used] = max(nxt[used], table[used-q]+max(0, cb-4+q))
        table = nxt
    return max(table[:checks+1])


def feedback_moment(c, h, j, tau):
    k, m = max(0, c-h), 16-j
    if j == 0:
        return tau**(-k/2)
    if h >= m or k == 0:
        return 1.0
    checks, slots = m-h, 16-h
    total = 0.0
    for d in range(max(0, checks+k-slots), min(checks, k)+1):
        cap = 1.0 if d == k else min(1.0, pq(k)*255.0**(-checks)/pq(k-d))
        total += cap*tau**(-d/2)
    return total


def flag_scalar_product(left, right):
    a, b = left
    c, d = right
    return a*c, b*(c+d)+a*d


def block_caps(z, tau):
    quiet, live, terminal = [], [], []
    for ell in range(17):
        q = [0.0, 0.0]
        t = [0.0, 0.0]
        w = [0.0, 0.0]
        for p in range(ell+1):
            weight = math.comb(ell, p)*z**(ell-p)*(1-z)**p
            d = p+min(4, ell-p)
            q[d < ell] += weight*tau**((ell-d)/2)
            t[p < ell] += weight*tau**((ell-p)/2)
            for v in range(min(3, 16-ell)+1):
                weight_out = (math.comb(ell, p)*math.comb(16-ell, v)
                              *z**(16-p-v)*(1-z)**(p+v))
                k = min(4, p+v)
                d = p+v-k
                w[d < ell] += weight_out*tau**((ell-d)/2-k)
        quiet.append(q)
        terminal.append(t)
        live.append(w)
    return quiet, live, terminal


def tables(z, tau):
    shapes = list(itertools.combinations_with_replacement(range(17), 4))
    shape_index = {s: i for i, s in enumerate(shapes)}
    n = len(shapes)
    occ = np.array(shapes, dtype=np.int16)
    ell = occ.sum(axis=1)
    c = np.minimum(occ, 4).sum(axis=1)
    occupied = (occ > 0).sum(axis=1)
    growth = np.power(tau, ell/2)
    local = np.zeros((n, 4, 17, 17))
    terminal = np.zeros((n, 4, 17, 1))
    moment = np.zeros((17, 17, 17))
    for ci in range(17):
        for h in range(1, 17):
            for j in range(17):
                moment[ci, h, j] = feedback_moment(ci, h, j, tau)
    for h in range(1, 17):
        blocks = (h+3)//4
        floor = np.maximum(0, 13*blocks-occ[:, 4-blocks:].sum(axis=1))
        emit = growth*np.power(z, floor)
        category = np.where(ell == 0, 2, 3)
        terminal[np.arange(n), category, h, 0] = emit
        for j in range(17):
            if j == 0:
                feasible = (c > 0) & (h <= 4*occupied)
            else:
                feasible = h+j >= 17-4*occupied
            value = emit*moment[c, h, j]*tau**(-min(16-j, h))*math.comb(16, j)
            local[np.arange(n), category, h, j] = value*feasible
    acache = {}
    for i, shape in enumerate(shapes):
        guard()
        key = tuple(min(4, v) for v in shape)
        if key not in acache:
            ci = sum(key)
            row = np.zeros(17)
            if ci:
                blocks = sum(v > 0 for v in key)
                for j in range(max(1, 17-4*blocks), 17):
                    m = 16-j
                    lo = max(0, ci-j, activation_floor(key, m))
                    total = 0.0
                    for d in range(lo, min(m, ci-1)+1):
                        cap = min(1.0, pq(ci)*255.0**(-m)/pq(ci-d))
                        total += cap*tau**(-d/2)
                    row[j] = math.comb(16, j)*z*total
            acache[key] = row
        local[i, 1, 0] = growth[i]*acache[key]
    quiet, live, final_quiet = block_caps(z, tau)
    for i, shape in enumerate(shapes):
        q = w = t = (1.0, 0.0)
        for l in shape:
            q = flag_scalar_product(q, quiet[l])
            w = flag_scalar_product(w, live[l])
            t = flag_scalar_product(t, final_quiet[l])
        local[i, 0, 0, 0], local[i, 1, 0, 0] = q
        local[i, 2, 16, 16], local[i, 3, 16, 16] = w
        terminal[i, 0, 0, 0], terminal[i, 1, 0, 0] = t
        terminal[i, 2, 16, 0], terminal[i, 3, 16, 0] = w
    ordered_shapes = np.fromiter((shape_index[tuple(sorted(x))]
                                  for x in itertools.product(range(17), repeat=4)),
                                 dtype=np.int32, count=17**4)
    geometry = []
    fac = math.factorial
    for a in range(17):
        for b in range(17):
            terms = []
            for f in range(max(0, a+b-16), min(a, b)+1):
                single, empty = a+b-2*f, 16-a-b+f
                coefficient = fac(16)//(fac(f)*fac(a-f)*fac(b-f)*fac(empty))
                terms.append((single, f, empty, coefficient/2**single))
            geometry.append(terms)
    return local, terminal, ordered_shapes, geometry


def to_log(a):
    result = np.full(a.shape, INF)
    np.log(a, out=result, where=a > 0)
    return result


def logmm(a, b):
    return np.logaddexp.reduce(a[:, :, None]+b[None, :, :], axis=1)


def multiply(a, b):
    out = np.full((4, a.shape[1], b.shape[2]), INF)
    for i in range(4):
        if not np.isfinite(a[i]).any():
            continue
        for j in range(4):
            if np.isfinite(b[j]).any():
                out[i | j] = np.logaddexp(out[i | j], logmm(a[i], b[j]))
    return out


def powers127_128(a):
    identity = np.full(a.shape, INF)
    np.fill_diagonal(identity[0], 0)
    p127, power = identity, a
    for bit in range(7):
        guard()
        p127 = multiply(p127, power)
        power = multiply(power, power)
    return p127, power


def profile_weights(r, f, scale, geometry, ordered):
    single, empty = r-2*f, 4096-r+f
    u, v = single/empty*scale, f/empty*scale**2
    p0 = 1/(1+u+v)
    ps, pf = u*p0, v*p0
    half = np.array([sum(coeff*ps**s*pf**f0*p0**e
                         for s, f0, e, coeff in terms) for terms in geometry])
    ordered_weights = np.multiply.outer(half, half).ravel()
    weights = np.bincount(ordered, weights=ordered_weights, minlength=4845)
    log_probability = (math.lgamma(4097)-math.lgamma(single+1)
                       -math.lgamma(f+1)-math.lgamma(empty+1)+empty*math.log(p0))
    if single:
        log_probability += single*math.log(ps)
    if f:
        log_probability += f*math.log(pf)
    assert abs(weights.sum()-1) < 2e-11
    return weights, log_probability


def regional_envelopes(args, local, terminal, ordered, geometry):
    regional = np.full((4, 17, 17), INF)
    final = np.full((4, 17, 1), INF)
    origin = np.full((4, 17, 17), -1, dtype=np.int16)
    profiles = range(args.r//2+1) if args.profiles == "all" else sorted(
        set(int(x) for x in args.profiles.split(",")))
    processed = []
    for f in profiles:
        guard()
        if f < 0 or f > args.r//2:
            raise ValueError("profile f outside support range")
        weights, denominator = profile_weights(args.r, f, args.scale, geometry, ordered)
        packet = to_log(np.tensordot(weights, local, axes=(0, 0)))
        terminal_packet = to_log(np.tensordot(weights, terminal, axes=(0, 0)))
        p127, p128 = powers127_128(packet)
        tail = multiply(p127, terminal_packet)
        p128 -= denominator
        tail -= denominator
        # Coefficient-support pruning: a nonempty region with full tags stays0.
        p128[0:3:2] = INF
        p128[0, 0, 0] = 0.0
        tail[0:3:2] = INF
        tail[0, 0, 0] = 0.0
        origin[p128 > regional] = f
        regional = np.maximum(regional, p128)
        final = np.maximum(final, tail)
        processed.append(f)
    leading = sorted((float(regional[3, h, h]/LN2), h, int(origin[3, h, h]))
                     for h in range(17))[-5:][::-1]
    provenance = {"largest_nonquiet_deficient_diagonals": [
        {"state_weight": h, "log2_mass": value, "profile_f": f}
        for value, h, f in leading],
        "quiet_deficient00_log2_mass": float(regional[1, 0, 0]/LN2),
        "quiet_deficient00_profile_f": int(origin[1, 0, 0])}
    return regional, final, processed, provenance


def evaluate(args, region, terminal, logy):
    column = multiply(multiply(region, region), multiply(region, region))
    end = multiply(multiply(multiply(region, region), region), terminal)
    a0 = np.logaddexp(column[0], column[1]+logy)
    a1 = np.logaddexp(column[2], column[3]+logy)
    total = np.logaddexp(a0, a1)
    v0 = np.full((1, 17), INF)
    v0[0, 0] = 0.0
    v1 = np.full((1, 17), INF)
    for _ in range(63):
        v1 = np.logaddexp(logmm(v1, total), logmm(v0, a1))
        v0 = logmm(v0, a0)
    last0 = np.logaddexp(end[0], end[1]+logy)
    last1 = np.logaddexp(end[2], end[3]+logy)
    quiet = float(logmm(v0, last0)[0, 0])
    nonquiet = float(np.logaddexp(logmm(v1, np.logaddexp(last0, last1)),
                                logmm(v0, last1))[0, 0])
    limit = (35*(1 << 21))//100
    prefix = (math.log(math.comb(8192, args.r))-math.log(args.tau-1)
              -33*logy-limit*math.log(args.z))
    return {"log2_y": logy/LN2, "nonquiet_log2_bound": (prefix+nonquiet)/LN2,
            "quiet_log2_bound": (prefix+quiet)/LN2}


def main():
    global DEADLINE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r", type=int, choices=(256, 512), default=256)
    parser.add_argument("--z", type=float, default=.8)
    parser.add_argument("--tau", type=float, default=256)
    parser.add_argument("--scale", type=float, default=1)
    parser.add_argument("--profiles", default="all")
    parser.add_argument("--log2-y", default="0,8,16,32")
    args = parser.parse_args()
    if not (.1 <= args.z < 1 and 1 < args.tau <= 65280 and 1e-6 <= args.scale <= 1e6):
        raise ValueError("bounded parameter range violated")
    ys = [float(x) for x in args.log2_y.split(",")]
    if len(ys) > 8 or any(not 0 <= y <= 256 for y in ys):
        raise ValueError("at most8 nonnegative marker exponents<=256")
    free = shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free < 20*(1 << 30):
        raise RuntimeError("requires20 GiB free on actual source volume")
    start = time.monotonic()
    DEADLINE = start+48
    local, terminal, ordered, geometry = tables(args.z, args.tau)
    prepared = time.monotonic()-start
    region, final, profiles, provenance = regional_envelopes(args, local, terminal, ordered, geometry)
    results = [evaluate(args, region, final, y*LN2) for y in ys]
    print(json.dumps({"scope": "FLOAT actual-route rank-cap diagnostic; not a certificate",
                      "whole_code_certificate": False, "r": args.r, "z": args.z,
                      "tau": args.tau, "marker_scale": args.scale,
                      "profiles_count": len(profiles), "all_profiles": len(profiles) == args.r//2+1,
                      "profile_range": [min(profiles), max(profiles)],
                      "positive_activity_split": True, "deficient_column_cutoff": 33,
                      "local_upgrades": ["quiet00", "live16to16", "output_only_terminal"],
                      "array_bytes": local.nbytes+terminal.nbytes+ordered.nbytes,
                      "table_seconds": prepared, "seconds": time.monotonic()-start,
                      "generated_bytes": 0, "regional_provenance": provenance,
                      "results": results}, indent=2))


if __name__ == "__main__":
    main()
