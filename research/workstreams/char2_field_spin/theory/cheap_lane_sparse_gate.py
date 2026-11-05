"""Bounded one-actual-group proof gate; stdout only.

The scalar-shadow mode preserves the failed 70-state diagnostic. The
density-free mode uses P_tau(d)=(tau**d-1)/(tau-1), actual state supports,
and four source ranks per wholly quiet outer column. It has 17 states.
Its --verify mode evaluates fixed rational parameters with directed Arb
arithmetic. This gate covers all GF128 messages in exactly one source
group, unioned over8192 groups; it is not a whole-code certificate.

No optimizer, bulk output, benchmark, or files are generated. Each call
evaluates at most 64 pairs. Floating-point output is a proposal only:
very small intermediate paths can underflow and make it optimistic.
Only --verify supplies a directed arithmetic check.
"""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
import shutil

import numpy as np

Q = 256
SHAPES = list(itertools.combinations_with_replacement(range(5), 4))
MULT = np.array([
    math.factorial(4) / math.prod(math.factorial(p.count(j)) for j in range(5))
    * math.prod(math.comb(4, j) for j in p) for p in SHAPES
])
assert sum(MULT) == 65536 and SHAPES[0] == (0, 0, 0, 0)


def subset_tail(n, z):
    return sum(math.comb(n, j) * z ** (n-j) * (1-z) ** j
               for j in range(4))


def matrices(alpha, z):
    size = len(SHAPES)
    empty = np.zeros((size, size))
    occupied = np.zeros((size, size))
    final_empty = np.zeros(size)
    final_occupied = np.zeros(size)
    a16, a15 = subset_tail(16, z), subset_tail(15, z)
    for i, p in enumerate(SHAPES):
        h = sum(p)
        blocks = sum(x != 0 for x in p)
        zero = Q ** (-h)
        moment = math.prod(Q**(-x) + (1-Q**(-x))*z**13 for x in p)
        inp_moment = (moment / Q + (1-1/Q) *
                      (zero*z + (moment-zero)/z))
        final_empty[i] = min(moment, a16**blocks)**alpha
        final_occupied[i] = sum(
            min(inp_moment, a16**(blocks-(b != 0)) * a15**(b != 0))**alpha
            for b in p) / 4
        for j, pp in enumerate(SHAPES):
            hh = sum(pp)
            if h == 0:
                if hh == 0:
                    empty[i, j] = 1
                    occupied[i, j] = Q**(-alpha)
                elif hh >= 13:
                    # Actual activation requires nonzero x, so requested
                    # output zeros exclude its position (factor z).
                    k = z * (1/Q + (1-1/Q)*255**(hh-16))
                    occupied[i, j] = MULT[j] * k**alpha
                continue
            if hh != 0 and h+hh >= 17:
                raw = zero + 255**(hh-16)*(moment-zero)
                restricted = Q**(hh-16) * a16**blocks
                empty[i, j] = MULT[j]*min(raw, restricted)**alpha
            if hh == 0:
                # Cancellation by one input is possible only inside its
                # one state block; omit additional column-pattern tests.
                if blocks == 1:
                    occupied[i, j] = (Q**(-h)*a15)**alpha / 4
                continue
            raw = zero + 255**(hh-16)*inp_moment
            for b in p:
                if h-b+hh < 13:
                    continue
                restricted = (Q**(-min(16-hh, h)) *
                              a16**(blocks-(b != 0)) * a15**(b != 0))
                occupied[i, j] += MULT[j]*min(raw, restricted)**alpha / 4
    return empty, occupied, final_empty, final_occupied


def powers_and_sum(e, i, n):
    """Return E^n and sum_{j<n} E^j I E^(n-1-j)."""
    if n == 0:
        return np.eye(len(e)), np.zeros_like(e)
    if n == 1:
        return e.copy(), i.copy()
    p, s = powers_and_sum(e, i, n//2)
    pp, ss = p@p, s@p+p@s
    if n % 2:
        return pp@e, ss@e+pp@i
    return pp, ss


def gate(alpha, z, distance):
    e, i, ve, vi = matrices(alpha, z)
    p127, s127 = powers_and_sum(e, i, 127)
    region = (s127@e+p127@i)/128
    final = (s127@ve+p127@vi)/128
    v = np.zeros(len(e))
    v[0] = 1
    logscale = 0.0
    for _ in range(255):
        v = v@region
        scale = np.max(v)
        if not math.isfinite(scale) or scale <= 0:
            return {"alpha": alpha, "z": z, "distance": distance,
                    "error": "numeric range exceeded"}
        v /= scale
        logscale += math.log2(scale)
    logmoment = logscale+math.log2(v@final)
    threshold = math.floor(distance*(1 << 21))
    logbeta = 160*8-32*math.log2(255)
    bound = 13+alpha*(logbeta-math.log2(255)-threshold*math.log2(z))+logmoment
    return {"alpha": alpha, "z": z, "distance": distance,
            "threshold": threshold, "log2_moment": logmoment,
            "log2_gate": bound, "margin_diagnostic": -bound,
            "states": len(SHAPES),
            "status": "floating proposal only; intermediate underflow can underestimate"}


def density_free_gate(tau, z, distance, activation_d1=False, truncate_quiet=False,
                      extinction_position=False, partial_quiet=False):
    """Audited projective rank potential, with no scalar outer density.

    State h counts the actual active state coordinates. A whole quiet
    outer column imposes four deterministic source ranks. The fictitious
    terminal state records H16 D1(Q+CX); no terminal lambda is charged.
    """
    size = 17
    e, i = np.zeros((size, size)), np.zeros((size, size))
    e[0, 0] = i[0, 0] = 1
    for h in range(size):
        for hh in range(size):
            if h == hh == 0:
                continue
            mult = math.comb(16, hh)
            if h and hh and h+hh >= 17:
                e[h, hh] = mult*tau**(-(16-hh))*z**(13*((h+3)//4))
            if h == 0:
                if hh >= 13:
                    i[h, hh] = mult*z*(255**(hh-16) if activation_d1 else 1)
            elif hh == 0:
                if h <= 4:
                    ext = ({1: 1/64, 2: 1/64, 4: 14/64}.get(h, 0)
                           if extinction_position else 1)
                    i[h, hh] = ext*tau**(-h)*z**(13*((h+3)//4)-1)
            elif h+hh >= 13:
                i[h, hh] = (mult*tau**(-min(16-hh, h))*
                            z**(13*((h+3)//4)-1))
    _, s128 = powers_and_sum(e, i, 128)
    region = s128/128
    assert region[0, 0] == 1
    # Compute G^4-e00 without cancellation of a tiny quiet-column weight.
    remainder = region.copy()
    remainder[0, 0] = 0
    first_row = np.zeros_like(region)
    first_row[0, :] = remainder[0, :]
    for _ in range(3):
        remainder = remainder@region+first_row
    column = remainder.copy()
    column[0, 0] += tau**(-4)
    if partial_quiet:
        quiet = np.zeros_like(region)
        quiet[0, 0] = 1
        nonquiet = region.copy()
        nonquiet[0, 0] = 0
        coeff = [np.eye(size)]
        for _ in range(4):
            nxt = [np.zeros_like(region) for _ in range(len(coeff)+1)]
            for k, mat in enumerate(coeff):
                nxt[k] += mat@nonquiet
                nxt[k+1] += mat@quiet
            coeff = nxt
        v = np.zeros((32, 4, size))
        v[0, 0, 0] = 1
    elif truncate_quiet:
        v = np.zeros((32, size))
        v[0, 0] = 1
    else:
        v = np.zeros(size)
        v[0] = 1
    logscale = 0.0
    for _ in range(64):
        if partial_quiet:
            vv = np.zeros_like(v)
            for k in range(4):
                product = v@coeff[k]
                for m in range(4):
                    vv[:, max(m, k), :] += product[:, m, :]*tau**(-max(0, k-m))
            vv[1:, :, 0] += v[:-1, :, 0]*tau**(-4)
            v = vv
        elif truncate_quiet:
            vv = v@remainder
            vv[1:, 0] += v[:-1, 0]*tau**(-4)
            v = vv
        else:
            v = v@column
        scale = np.max(v)
        v /= scale
        logscale += math.log2(scale)
    if partial_quiet:
        v *= np.array([[1-tau**(-(128-4*c-m)) for m in range(4)]
                       for c in range(32)])[:, :, None]
    elif truncate_quiet:
        v *= np.array([1-tau**(-(128-4*c)) for c in range(32)])[:, None]
    logmoment = logscale+math.log2(np.sum(v))
    threshold = math.floor(distance*(1 << 21))
    bound = (13+128*math.log2(tau)-math.log2(tau-1)
             -threshold*math.log2(z)+logmoment)
    return {"method": "density-free", "tau": tau, "z": z,
            "distance": distance, "threshold": threshold,
            "log2_moment": logmoment, "log2_gate": bound,
            "margin_diagnostic": -bound, "states": size,
            "activation_d1": activation_d1, "truncate_quiet": truncate_quiet,
            "extinction_position": extinction_position,
            "partial_quiet": partial_quiet,
            "status": "floating proposal only; intermediate underflow can underestimate"}


def verify_density_free(tau_text, z_text, distance_text, precision, target,
                        activation_d1=False, truncate_quiet=False,
                        extinction_position=False, partial_quiet=False):
    """Directed evaluation of fixed rational parameters; no optimizer."""
    from flint import arb, arb_mat, ctx
    ctx.prec = precision
    tr, zr, dr = map(Fraction, (tau_text, z_text, distance_text))
    if not (1 < tr <= 65280 and 0 < zr < 1 and 0 < dr < Fraction(1, 2)):
        raise ValueError("invalid rational gate parameters")
    tau = arb(tr.numerator)/tr.denominator
    z = arb(zr.numerator)/zr.denominator
    n = 17
    e, i = arb_mat(n, n), arb_mat(n, n)
    e[0, 0] = i[0, 0] = 1
    for h in range(n):
        for hh in range(n):
            if h == hh == 0:
                continue
            mult = math.comb(16, hh)
            if h and hh and h+hh >= 17:
                e[h, hh] = mult*tau**(-(16-hh))*z**(13*((h+3)//4))
            if h == 0:
                if hh >= 13:
                    i[h, hh] = mult*z*(arb(255)**(hh-16) if activation_d1 else 1)
            elif hh == 0:
                if h <= 4:
                    ext = (arb({1: 1, 2: 1, 4: 14}.get(h, 0))/64
                           if extinction_position else arb(1))
                    i[h, hh] = ext*tau**(-h)*z**(13*((h+3)//4)-1)
            elif h+hh >= 13:
                i[h, hh] = (mult*tau**(-min(16-hh, h))*
                            z**(13*((h+3)//4)-1))

    def ps(k):
        if k == 1:
            return e, i
        p, s = ps(k//2)
        pp, ss = p*p, s*p+p*s
        return (pp*e, ss*e+pp*i) if k % 2 else (pp, ss)

    _, s128 = ps(128)
    region = s128/128
    if region[0, 0] != arb(1):
        raise ArithmeticError("the unique all-quiet regional path must have exact weight one")
    remainder = arb_mat(region)
    remainder[0, 0] = 0
    first = arb_mat(n, n)
    for j in range(n):
        first[0, j] = remainder[0, j]
    for _ in range(3):
        remainder = remainder*region+first
    column = arb_mat(remainder)
    column[0, 0] += tau**(-4)
    if partial_quiet:
        quiet, nonquiet, identity = arb_mat(n, n), arb_mat(region), arb_mat(n, n)
        quiet[0, 0] = 1
        nonquiet[0, 0] = 0
        for j in range(n):
            identity[j, j] = 1
        coeff = [identity]
        for _ in range(4):
            nxt = [arb_mat(n, n) for _ in range(len(coeff)+1)]
            for k, mat in enumerate(coeff):
                nxt[k] += mat*nonquiet
                nxt[k+1] += mat*quiet
            coeff = nxt
    v = arb_mat(128 if partial_quiet else (32 if truncate_quiet else 1), n)
    v[0, 0] = 1
    for _ in range(64):
        if partial_quiet:
            vv = arb_mat(128, n)
            for k in range(4):
                product = v*coeff[k]
                for c in range(32):
                    for m in range(4):
                        target_row = 4*c+max(m, k)
                        factor = tau**(-max(0, k-m))
                        for j in range(n):
                            vv[target_row, j] += product[4*c+m, j]*factor
            for c in range(1, 32):
                for m in range(4):
                    vv[4*c+m, 0] += v[4*(c-1)+m, 0]*tau**(-4)
            v = vv
        elif truncate_quiet:
            vv = v*remainder
            for c in range(1, 32):
                vv[c, 0] += v[c-1, 0]*tau**(-4)
            v = vv
        else:
            v = v*column
    if partial_quiet:
        moment = sum((1-tau**(-(128-4*c-m)))*v[4*c+m, j]
                     for c in range(32) for m in range(4) for j in range(n))
    elif truncate_quiet:
        moment = sum((1-tau**(-(128-4*c)))*v[c, j]
                     for c in range(32) for j in range(n))
    else:
        moment = sum(v[0, j] for j in range(n))
    threshold = int(dr*(1 << 21))
    bound = 8192*tau**128/(tau-1)*z**(-threshold)*moment
    passed = bool(bound.is_finite() and bound < arb(2)**(-target))
    return {"method": "density-free-directed", "tau": str(tr),
            "z": str(zr), "distance": str(dr), "threshold": threshold,
            "precision": precision, "target_bits": target, "passed": passed,
            "log2_bound_interval": str(bound.log()/arb(2).log()),
            "scope": "all nonzero GF128 messages in exactly one of 8192 source groups",
            "activation_d1": activation_d1, "truncate_quiet": truncate_quiet,
            "extinction_position": extinction_position,
            "partial_quiet": partial_quiet,
            "states": n}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, nargs="+", default=[1.0])
    ap.add_argument("--z", type=float, nargs="+", default=[0.999])
    ap.add_argument("--distance", type=float, default=0.35)
    ap.add_argument("--method", choices=["scalar-shadow", "density-free"],
                    default="scalar-shadow")
    ap.add_argument("--tau", type=float, nargs="+", default=[65280])
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--precision", type=int, default=192)
    ap.add_argument("--target-bits", type=int, default=180)
    ap.add_argument("--activation-d1", action="store_true")
    ap.add_argument("--truncate-quiet", action="store_true")
    ap.add_argument("--extinction-position", action="store_true")
    ap.add_argument("--partial-quiet", action="store_true")
    args = ap.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        ap.error("at least 20 GiB of local free space is required for numerical runs")
    if args.partial_quiet:
        args.truncate_quiet = True
    if not (64 <= args.precision <= 512 and 0 < args.target_bits <= 10000):
        ap.error("precision64..512 and target bits1..10000 required")
    if not (0 < args.distance < 0.5):
        ap.error("distance must lie strictly between zero and one half")
    if args.verify and args.method != "density-free":
        ap.error("directed verification is available only for density-free mode")
    parameters = args.alpha if args.method == "scalar-shadow" else args.tau
    if len(parameters)*len(args.z) > 64:
        ap.error("at most 64 bounded parameter pairs")
    for a in parameters:
        for z in args.z:
            if args.method == "density-free":
                if not (1 < a <= 65280 and 0 < z < 1):
                    ap.error("1 < tau <= 65280 and 0 < z < 1 required")
                if args.verify:
                    result = verify_density_free(str(a), str(z), str(args.distance),
                                                 args.precision, args.target_bits,
                                                 args.activation_d1, args.truncate_quiet,
                                                 args.extinction_position, args.partial_quiet)
                    print(json.dumps(result, sort_keys=True))
                    if not result["passed"]:
                        raise SystemExit(1)
                else:
                    print(json.dumps(density_free_gate(a, z, args.distance,
                                                      args.activation_d1, args.truncate_quiet,
                                                      args.extinction_position, args.partial_quiet), sort_keys=True))
                continue
            if not (0 < a <= 1 and 0 < z < 1):
                ap.error("0 < alpha <= 1 and 0 < z < 1 required")
            print(json.dumps(gate(a, z, args.distance), sort_keys=True))


if __name__ == "__main__":
    main()
