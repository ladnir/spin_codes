"""Bounded FLOAT actual-support filtered scalar pilot, not a certificate.

The actual active source set T is retained for r<=q=256. On the initial
E-defined kernel, the hyperplane-union lemma gives the prefactor
q^2/((q-1)^2*(q-r+1)) times the exact-T scalar count. Each active group's
outer density costs beta=q^160/(q-1)^32. The source-support union remains
outside the fractional power: C(8192,r)*(prefactor*beta^r*z^-L)^alpha.

The per-mask scalar kernel is the positive output-zero sum of
q^-rank[M,N_I]. Actual-path impossibilities are removed BEFORE relaxing
to auxiliary E vectors. Empty input cannot activate or extinguish; live
paths require h+j>=17-4*occupied when j>0. In an active RS state block,
requested zeros outside actual I must number at most3. Auxiliary states
may still have zero coordinates. No nonzero-state assumption is made
about the auxiliary scalar vectors.

For live states the local bound is a minimum of three bounds on this SAME
filtered scalar kernel: the unfiltered g-padding bound; the filtered
output-only moment; and the feedback-only bound times retained zero-set
mass. The latter is valid termwise before averaging D1. Maximum-product
DP over h_b in0..4 gives cardinality-uniform bounds for h=sum h_b.

For an active block with ell inputs and n=16-ell outside positions, put
R=sum_{v<=3} C(n,v)z^(n-v)(1-z)^v. Its output-only moment is bounded by
g^ell times the minimum of sum_v weight*q^-max(0,h_b-4+v) and
q^-h_b R+(1-q^-h_b)z^max(0,n-3). The second expression separates the
zero auxiliary block from nonzero RS words; the latter have at most3
roots. For h_b=0 the factors are g^ell and1.
For h_b=1 the implemented basis column has at most two roots. Requiring
the outside-I zero subset to lie in that actual root set gives the
additional bound g^ell*z^max(0,16-ell-2), without a zero-state residue.
The other caps still bound this further-restricted positive object.

Activation uses the JOINT cap after requested input zeros P are removed:
E[q^-d 1(d<c)] <= min(1,q/(q-1)*(1-q^-c)*255^-(16-j)). The p eliminated
inputs contribute q^-p. Remaining block occupancies use weights
C(ell,k)z^k*((1-z)/q)^(ell-k). A positive DP tracks q^-c and1-q^-c,
together with the number of remaining occupied blocks. Thus no activation
and rank marginal probabilities are multiplied. Quiet00 is exact.

Every destination multiplicity C(16,j) is OUTSIDE K^alpha. Routing slot
multiplicities also stay outside alpha. One common z and alpha are fixed
before setup/zero-set averaging. The last packet is output-only, with no
next-state count or feedback factor. A positive activity flag reports
quiet and nonquiet paths separately, without subtracting close values.

The optional fixed gamma uses z=a*b, a=z^gamma and b=z^(1-gamma).
Actual weight w>=floor gives z^w<=a^floor*b^w BEFORE scalar relaxation.
Only b is expanded into requested zeros. Active-block factors include
a^max(0,13-ell), and activation pays a because its input/output is nonzero.
The global prefactor remains z^-L. gamma1 is the floor-only endpoint;
it is not obtained by multiplying two incompatible relaxed measures.

All paired occupancy tails are retained. Each f-profile has independent
coefficient markers, and its normalized denominator is the trinomial
profile probability. Optional bounded marker search minimizes a Perron
surrogate only; reported values use the actual128/127 packet powers.
--refine-quiet-markers separately minimizes the exact regional quiet00
coefficient using its own markers, before the profile envelope. This is
the same positive coefficient, so entry-specific marker choices are
valid. The full regional/terminal transfer is also evaluated at those
markers; every coefficient is replaced by the smaller of its two valid
marker bounds. The output-only final quiet region is at most
q^(16*alpha) times the regular quiet region, because the omitted feedback
has at most16 rows. This gives an additional bound for terminal00.

--activity-mark v weights each packet whose actual INCOMING state is
nonzero, including the final packet. Its factor is OUTSIDE alpha. The
reported bounds are raw marked sums, with no automatic event claim. For
0<v<1, a restriction A<=B on the number A of such packets costs v^-B;
for v>1, a restriction A>=B+1 gains v^-(B+1). These statements are valid
only with the corresponding restriction. z=1 is allowed for an existence
count with no output-weight restriction. Marked first-active partitions
retain their interpretation as positive parts of the marked sum.
Selected-profile results are NOT full support bounds. Entrywise regional
envelopes over all f are available, but remain FLOAT diagnostics.

No data/caches are written; source arrays<100MiB, free-space floor20GiB,
hard invocation budget48s. Frozen helpers are read only after hash check.
"""

import argparse
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

Q = 256
LN2 = math.log(2)
INF = float("-inf")
DEADLINE = 0.0
REFERENCE_SHA256 = "1b73c7736e1b8be40c860928d2d71454555ed2b283c3357fcb4f57f78be425b6"


def guard():
    if time.monotonic() > DEADLINE:
        raise RuntimeError("bounded48-second invocation exhausted")


def reference():
    path = Path(__file__).with_name("cheap_lane_middle_pilot.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != REFERENCE_SHA256:
        raise RuntimeError("frozen geometry helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_middle_pilot", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def geometry():
    shapes = list(itertools.combinations_with_replacement(range(17), 4))
    index = {shape: i for i, shape in enumerate(shapes)}
    ordered = np.fromiter((index[tuple(sorted(x))]
                           for x in itertools.product(range(17), repeat=4)),
                          dtype=np.int32, count=17**4)
    terms = []
    fact = math.factorial
    for a in range(17):
        for b in range(17):
            cell = []
            for f in range(max(0, a+b-16), min(a, b)+1):
                s, e = a+b-2*f, 16-a-b+f
                coefficient = fact(16)//(fact(f)*fact(a-f)*fact(b-f)*fact(e))
                cell.append((s, f, e, coefficient/2**s))
            terms.append(cell)
    return np.array(shapes, dtype=np.int16), ordered, terms


def max_product(occ, block):
    """Max over all 4-block state sizes; no unproved fixed profile."""
    n = len(occ)
    dp = np.zeros((n, 17))
    dp[:, 0] = 1.0
    for b in range(4):
        nxt = np.zeros_like(dp)
        for h in range(4*b+1):
            for hb in range(5):
                nxt[:, h+hb] = np.maximum(nxt[:, h+hb], dp[:, h]*block[occ[:, b], hb])
        dp = nxt
    return dp


def local_tables(occ, z, alpha, gamma, activity_mark):
    n = len(occ)
    ell = occ.sum(axis=1)
    c = np.minimum(occ, 4).sum(axis=1)
    occupied = (occ > 0).sum(axis=1)
    floor_mark = z**gamma
    z = z**(1-gamma)
    g = z+(1-z)/Q
    block_out = np.zeros((17, 5))
    block_keep = np.ones((17, 5))
    residual = np.zeros((17, 2, 2))  # occupied flag, Q^-c / 1-Q^-c
    for l in range(17):
        outside = 16-l
        rv = [math.comb(outside, v)*z**(outside-v)*(1-z)**v
              for v in range(min(3, outside)+1)]
        retained = sum(rv)
        block_out[l, 0] = g**l
        for hb in range(1, 5):
            rank_cap = sum(w*Q**(-max(0, hb-4+v)) for v, w in enumerate(rv))
            split_cap = Q**(-hb)*retained+(1-Q**(-hb))*z**max(0, outside-3)
            block_out[l, hb] = g**l*min(rank_cap, split_cap)
            if hb == 1:
                block_out[l, hb] = min(block_out[l, hb], g**l*z**max(0, outside-2))
            block_keep[l, hb] = retained
            block_out[l, hb] *= floor_mark**max(0, 13-l)
            block_keep[l, hb] *= floor_mark**max(0, 13-l)
        for k in range(l+1):
            w = math.comb(l, k)*z**k*((1-z)/Q)**(l-k)
            qc = Q**(-min(4, k))
            residual[l, int(k > 0), 0] += w*qc
            residual[l, int(k > 0), 1] += w*(1-qc)
    out = max_product(occ, block_out)
    retained = max_product(occ, block_keep)
    # Positive residual-rank DP, retaining the active-block count.
    qdp = np.zeros((n, 5))
    ddp = np.zeros_like(qdp)
    qdp[:, 0] = 1.0
    for b in range(4):
        nq, nd = np.zeros_like(qdp), np.zeros_like(ddp)
        for active in range(b+1):
            for add in (0, 1):
                qb = residual[occ[:, b], add, 0]
                db = residual[occ[:, b], add, 1]
                nq[:, active+add] += qdp[:, active]*qb
                nd[:, active+add] += ddp[:, active]*(qb+db)+qdp[:, active]*db
        qdp, ddp = nq, nd
    per_mask = np.zeros((n, 17, 17))
    per_mask[:, 0, 0] = qdp.sum(axis=1)
    for j in range(1, 17):
        minimum_blocks = max(1, (17-j+3)//4)
        if j == 16:
            value = (qdp[:, minimum_blocks:]+ddp[:, minimum_blocks:]).sum(axis=1)
        else:
            value = Q/(Q-1)*255.0**(-(16-j))*ddp[:, minimum_blocks:].sum(axis=1)
        per_mask[:, 0, j] = floor_mark*value
    for h in range(1, 17):
        guard()
        a0 = np.maximum(0, h-np.maximum(0, occ-12).sum(axis=1))
        blocks = (h+3)//4
        minimum_floor = np.maximum(0, 13-occ[:, 4-blocks:]).sum(axis=1)
        floor_factor = np.power(floor_mark, minimum_floor)
        k = np.maximum(0, c-h)
        for j in range(17):
            m, rf = 16-j, min(16-j, h)
            if j == 0:
                fm = np.power(float(Q), -k)
                feasible = (ell > 0) & (h <= 4*occupied)
            elif h >= m:
                fm = np.ones(n)
                feasible = h+j >= 17-4*occupied
            else:
                qk = np.power(float(Q), -k)
                fm = np.minimum(1.0, qk+(1-qk)*255.0**(-(m-h)))
                feasible = h+j >= 17-4*occupied
            feedback = Q**(-rf)*fm
            d0 = np.maximum(0, c-j-rf)
            e0 = np.maximum(0, a0-rf)
            padding = (floor_factor*Q**(-rf)*np.power(g, ell)
                       *np.power(Q*g, -d0)*np.power(g, e0))
            value = np.minimum(padding, np.minimum(out[:, h], feedback*retained[:, h]))
            per_mask[:, h, j] = value*feasible
    assert np.isfinite(per_mask).all() and (per_mask >= 0).all()
    # Destination mask multiplicity MUST remain outside the fractional power.
    matrix = np.power(per_mask, alpha)
    matrix *= np.array([math.comb(16, j) for j in range(17)])[None, None, :]
    local = np.zeros((n, 4, 17, 17))
    local[:, 0, 0, :] = matrix[:, 0, :]
    local[:, 2, 1:, :] = activity_mark*matrix[:, 1:, :]
    terminal = np.zeros((n, 4, 17, 1))
    terminal[:, 0, 0, 0] = np.power(g, alpha*ell)
    terminal[:, 2, 1:, 0] = activity_mark*np.power(out[:, 1:], alpha)
    return local, terminal


def weights(r, f, log_factors, ordered, terms):
    s, e = r-2*f, 4096-r+f
    u = s/e*math.exp(log_factors[0]) if s else 0.0
    v = f/e*math.exp(log_factors[1]) if f else 0.0
    p0 = 1/(1+u+v)
    ps, pf = u*p0, v*p0
    half = np.array([sum(coef*ps**single*pf**full*p0**empty
                         for single, full, empty, coef in cell) for cell in terms])
    shape = np.bincount(ordered, weights=np.multiply.outer(half, half).ravel(), minlength=4845)
    assert abs(shape.sum()-1) < 3e-11
    lp = math.lgamma(4097)-math.lgamma(s+1)-math.lgamma(f+1)-math.lgamma(e+1)+e*math.log(p0)
    if s:
        lp += s*math.log(ps)
    if f:
        lp += f*math.log(pf)
    return shape, lp


def select_markers(args, f, local, ordered, terms, quiet=False):
    logs = np.array([math.log(args.u_scale), math.log(args.v_scale)])
    calls = 0
    def objective(point):
        nonlocal calls
        guard()
        w, denom = weights(args.r, f, point, ordered, terms)
        if quiet:
            rho = float(np.dot(w, local[:, 0, 0, 0]))
        else:
            matrix = np.tensordot(w, local[:, 0]+local[:, 2], axes=(0, 0))
            rho = max(abs(np.linalg.eigvals(matrix)))
        calls += 1
        return 128*math.log(float(rho))-denom
    best = objective(logs)
    if args.optimize_markers or quiet:
        step = math.log(4)
        for iteration in range(8):
            changed = False
            for coordinate, enabled in ((0, args.r-2*f > 0), (1, f > 0)):
                if not enabled:
                    continue
                for sign in (-1, 1):
                    trial = logs.copy()
                    trial[coordinate] += sign*step
                    if not -18 <= trial[coordinate] <= 18:
                        continue
                    candidate = objective(trial)
                    if candidate < best:
                        best, logs, changed = candidate, trial, True
            if not changed:
                step *= .5
    return logs, calls, best/LN2


def regional(args, local, terminal, ordered, terms, ref):
    regions = np.full((4, 17, 17), INF)
    final = np.full((4, 17, 1), INF)
    profile_results = []
    profiles = range(args.r//2+1) if args.profiles == "all" else sorted(set(
        int(x) for x in args.profiles.split(",")))
    for f in profiles:
        guard()
        if not 0 <= f <= args.r//2:
            raise ValueError("profile outside allowed range")
        logs, calls, surrogate = select_markers(args, f, local, ordered, terms)
        w, denominator = weights(args.r, f, logs, ordered, terms)
        packet = ref.to_log(np.tensordot(w, local, axes=(0, 0)))
        last = ref.to_log(np.tensordot(w, terminal, axes=(0, 0)))
        p127, p128 = ref.powers127_128(packet)
        end = ref.multiply(p127, last)
        p128 -= denominator
        end -= denominator
        quiet_refinement = None
        if args.refine_quiet_markers:
            qlogs, qcalls, qbound = select_markers(args, f, local, ordered, terms, quiet=True)
            qw, qdenominator = weights(args.r, f, qlogs, ordered, terms)
            qpacket = ref.to_log(np.tensordot(qw, local, axes=(0, 0)))
            qlast = ref.to_log(np.tensordot(qw, terminal, axes=(0, 0)))
            qp127, qp128 = ref.powers127_128(qpacket)
            qend = ref.multiply(qp127, qlast)
            p128 = np.minimum(p128, qp128-qdenominator)
            end = np.minimum(end, qend-qdenominator)
            p128[0, 0, 0] = min(p128[0, 0, 0], qbound*LN2)
            end[0, 0, 0] = min(end[0, 0, 0], qbound*LN2+16*args.alpha*math.log(Q))
            quiet_refinement = {"u_scale": math.exp(qlogs[0]), "v_scale": math.exp(qlogs[1]),
                                "marker_calls": qcalls, "regular_log2_bound": qbound}
        regions = np.maximum(regions, p128)
        final = np.maximum(final, end)
        profile_results.append({"f": f, "u_scale": math.exp(logs[0]),
                                "v_scale": math.exp(logs[1]), "marker_calls": calls,
                                "perron_surrogate_log2": surrogate,
                                "live16_regional_log2": float(p128[2, 16, 16]/LN2),
                                "quiet00_regional_log2": float(p128[0, 0, 0]/LN2),
                                "quiet_refinement": quiet_refinement})
    return regions, final, profile_results


def bound(args, region, end, ref):
    column = ref.multiply(ref.multiply(region, region), ref.multiply(region, region))
    terminal_column = ref.multiply(ref.multiply(ref.multiply(region, region), region), end)
    quiet = np.full((1, 17), INF)
    quiet[0, 0] = 0.0
    live = np.full((1, 17), INF)
    total = np.logaddexp(column[0], column[2])
    for _ in range(63):
        live = np.logaddexp(ref.logmm(live, total), ref.logmm(quiet, column[2]))
        quiet = ref.logmm(quiet, column[0])
    quiet_result = float(ref.logmm(quiet, terminal_column[0])[0, 0])
    live_result = float(np.logaddexp(ref.logmm(live, np.logaddexp(terminal_column[0], terminal_column[2])),
                                   ref.logmm(quiet, terminal_column[2]))[0, 0])
    prefactor = Q*Q/((Q-1)**2*(Q-args.r+1))
    logbeta = 160*math.log(Q)-32*math.log(Q-1)
    limit = 35*(1 << 21)//100
    prefix = math.log(math.comb(8192, args.r))+args.alpha*(
        math.log(prefactor)+args.r*logbeta-limit*math.log(args.z))
    # Positive partition by the first column containing a nonzero incoming
    # packet state. Activation at the preceding column's final update is
    # first flagged here, not in the column that performed that update.
    suffix = [None]*64
    suffix[63] = np.logaddexp(terminal_column[0], terminal_column[2])
    for t in range(62, -1, -1):
        suffix[t] = ref.logmm(total, suffix[t+1])
    qprefix = np.full((1, 17), INF)
    qprefix[0, 0] = 0.0
    first = []
    for t in range(63):
        value = float(ref.logmm(ref.logmm(qprefix, column[2]), suffix[t+1])[0, 0])
        first.append(((prefix+value)/LN2, t))
        qprefix = ref.logmm(qprefix, column[0])
    value = float(ref.logmm(qprefix, terminal_column[2])[0, 0])
    first.append(((prefix+value)/LN2, 63))
    partition = float(np.logaddexp.reduce(np.array([value*LN2 for value, _ in first]))/LN2)
    assert math.isclose(partition, (prefix+live_result)/LN2, rel_tol=2e-12, abs_tol=2e-8)
    return {"nonquiet_log2_bound": (prefix+live_result)/LN2,
            "quiet_log2_bound": (prefix+quiet_result)/LN2,
            "source_log2_prefactor": prefix/LN2,
            "largest_first_active_columns": [
                {"column_zero_based": t, "log2_contribution": value}
                for value, t in sorted(first, reverse=True)[:5]],
            "sampled_first_active_columns": [
                {"column_zero_based": t, "log2_contribution": first[t][0]}
                for t in (*range(0, 64, 8), 63)],
            "first_active_partition_log2": partition}


def main():
    global DEADLINE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r", type=int, default=256)
    parser.add_argument("--z", type=float, default=.95)
    parser.add_argument("--alpha", type=float, default=.523)
    parser.add_argument("--gamma", type=float, choices=(0, .5, 1), default=0)
    parser.add_argument("--activity-mark", type=float, default=1)
    parser.add_argument("--profiles", default="0,64,128")
    parser.add_argument("--u-scale", type=float, default=1)
    parser.add_argument("--v-scale", type=float, default=1)
    parser.add_argument("--optimize-markers", action="store_true")
    parser.add_argument("--refine-quiet-markers", action="store_true")
    args = parser.parse_args()
    if not (1 <= args.r <= 256 and .001 <= args.z <= 1 and .1 <= args.alpha <= 1):
        raise ValueError("parameters outside bounded proof interface")
    if not 1e-6 <= args.activity_mark <= 1e6:
        raise ValueError("activity marker outside bounded range")
    if not (1e-6 <= args.u_scale <= 1e6 and 1e-6 <= args.v_scale <= 1e6):
        raise ValueError("coefficient marker scales out of bounds")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on actual volume")
    start = time.monotonic()
    DEADLINE = start+48
    ref = reference()
    ref.DEADLINE = DEADLINE
    occ, ordered, terms = geometry()
    local, terminal = local_tables(occ, args.z, args.alpha, args.gamma, args.activity_mark)
    prepared = time.monotonic()-start
    region, last, profiles = regional(args, local, terminal, ordered, terms, ref)
    result = bound(args, region, last, ref)
    print(json.dumps({"scope": "FLOAT filtered actual-support scalar pilot; not a certificate",
                      "whole_code_certificate": False, "r": args.r, "z": args.z,
                      "alpha": args.alpha, "gamma": args.gamma,
                      "activity_mark": args.activity_mark,
                      "marked_bound_requires_matching_activity_restriction": args.activity_mark != 1,
                      "all_profiles": len(profiles) == args.r//2+1,
                      "profiles": profiles, "result": result,
                      "array_bytes": local.nbytes+terminal.nbytes+ordered.nbytes+occ.nbytes,
                      "table_seconds": prepared, "seconds": time.monotonic()-start,
                      "generated_bytes": 0, "geometry_source_sha256": REFERENCE_SHA256}, indent=2))


if __name__ == "__main__":
    main()
