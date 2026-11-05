"""FLOAT singleton-extinction diagnostic for one pair-only regional marginal.

Start with the frozen block-support polynomial. Replace only regular live
transitions to zero; retain quiet00, live destinations, and the terminal.
The two blocks in a half-packet share their selected point. Their singleton
caps must be multiplied at that point BEFORE averaging. Incoming masks are
fixed during this average and maximized afterward, without a mask count.

For point a, let v(a)=[1,a,a(a+1)/6,a^2(a+1)/6] over AES GF256.
An active singleton block can extinguish only with H=supp(v(a)). Its
retained rank cap is q^-|H| z^wt((I+AA^t)e_a), q=256. There is no factor
1-1/q: the surviving scalar kernel still contains zero. A quiet singleton
block contributes q^-1. Requested zeros at nonzero positions of the emitted
column make the actual branch impossible and are removed before relaxation.

For ell>=2, an active block uses q^-max(h,min(4,ell))*z^max(0,13-ell).
A quiet block uses its exact positive MDS zero-expansion moment. For ell=0,
only the quiet state is possible. Define T_ell[s] as the maximum, over two
fixed incoming masks of total size s, of the route average of the product
of block caps raised to alpha. The half-packets are independent conditional
on their counts, so G_ab[h]=max_s T_a[s]T_b[h-s] bounds the packet moment.
Take its minimum with the old uniform powered cap and apply activity once.

Rebuild live-to-zero coefficients by a positive sum C(16,a)C(16,b), keyed
by the two half counts. No subtraction or independently averaged singleton
factors are used. Degree extraction and the boundary reuse frozen engines.
The whole-chain score repeats a pair-only marginal in both shifted
matchings, which is impossible for proper nonempty supports. It is an
artificial diagnostic, not a class bound or certificate, even if negative.

Only the fixed r256,z.85,alpha.523 point is available. Use python -B with
--self-test or --run. No generated data; guards are 20 GiB free,48 seconds,
arrays below256 MiB. All frozen imports are hash checked and unmodified.
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
BLOCK_SHA256 = "6c826ae6ff8e5cbfb256313b85f8155ab34b6f348a97604c4b0309a02e540b71"
Q = 256


def load_block():
    path = Path(__file__).with_name("cheap_lane_block_support.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != BLOCK_SHA256:
        raise RuntimeError("frozen block-support helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_block_support", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gf_mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 256:
            a ^= 0x11B
    return out


def singleton_columns():
    inverse6 = next(a for a in range(1, Q) if gf_mul(6, a) == 1)
    rows = []
    for a in range(16):
        third = gf_mul(gf_mul(a, a ^ 1), inverse6)
        rows.append((1, a, third, gf_mul(a, third)))
    def dot(x, y):
        value = 0
        for a, b in zip(x, y):
            value ^= gf_mul(a, b)
        return value
    for i in range(4):
        for j in range(4):
            assert dot(tuple(row[i] for row in rows),
                       tuple(row[j] for row in rows)) == 0
    columns = []
    for a, row in enumerate(rows):
        mask = sum(1 << i for i, x in enumerate(row) if x)
        weight = sum((dot(row, other) ^ int(a == b)) != 0
                     for b, other in enumerate(rows))
        columns.append((mask, weight))
    assert columns == [(1, 15), (3, 16), (15, 15), (15, 15)]+[(15, 16)]*12
    return columns


def singleton_formula(qpower, zpower):
    """Power parameters are q^-alpha,z^alpha; Fraction works at alpha=1."""
    v = [0]*9
    v[0] = qpower**2
    v[1] = qpower**2*zpower**15/16
    v[2] = max(qpower**2*zpower**30, qpower**3*zpower**16)/16
    v[4] = max(qpower**4*zpower**32,
               qpower**5*(2*zpower**15+12*zpower**16))/16
    v[8] = qpower**8*(2*zpower**30+12*zpower**32)/16
    return v


def quiet_mass(ell, z):
    return sum(math.comb(ell, k)*z**k*((1-z)/Q)**(ell-k)
               *Fraction(1, Q**min(4, k)) for k in range(ell+1))


def generic_block(ell, h, z):
    if h == 0:
        return quiet_mass(ell, z)
    if ell == 0:
        return 0
    if ell == 1:
        raise ValueError("a singleton requires its shared point")
    return Fraction(1, Q**max(h, min(4, ell)))*z**max(0, 13-ell)


def half_tables(z, alpha):
    table = [[0.0]*9 for _ in range(17)]
    table[0][0] = 1.0
    table[1] = singleton_formula(Q**(-alpha), z**alpha)
    for ell in range(2, 17):
        blocks = [float(generic_block(ell, h, z))**alpha for h in range(5)]
        for h1, h2 in itertools.product(range(5), repeat=2):
            table[ell][h1+h2] = max(table[ell][h1+h2], blocks[h1]*blocks[h2])
    return table


def packet_cap(halves, a, b, h):
    return max(halves[a][s]*halves[b][h-s]
               for s in range(max(0, h-8), min(8, h)+1))


def self_test():
    columns = singleton_columns()
    z = Fraction(17, 20)
    point_caps = []
    for mask, weight in columns:
        row = [Fraction(0)]*16
        row[0] = Fraction(1, Q)
        row[mask] = Fraction(1, Q**mask.bit_count())*z**weight
        point_caps.append(row)
    direct = [Fraction(0)]*9
    floating = [0.0]*9
    alpha = .523
    for h1, h2 in itertools.product(range(16), repeat=2):
        s = h1.bit_count()+h2.bit_count()
        products = [row[h1]*row[h2] for row in point_caps]
        direct[s] = max(direct[s], sum(products)/16)
        floating[s] = max(floating[s], sum(float(x)**alpha for x in products)/16)
    assert direct == singleton_formula(Fraction(1, Q), z)
    for actual, expected in zip(floating, singleton_formula(Q**(-alpha), float(z)**alpha)):
        assert math.isclose(actual, expected, rel_tol=3e-13, abs_tol=1e-100)
    quiet_checks = 0
    for ell in range(17):
        assert quiet_mass(ell, Fraction(0)) == Fraction(1, Q**ell)
        assert quiet_mass(ell, Fraction(1)) == Fraction(1, Q**min(4, ell))
        if ell <= 4:
            assert quiet_mass(ell, z) == Fraction(1, Q**ell)
        quiet_checks += 1
    halves = half_tables(float(z), alpha)
    checks = 0
    for a, b, h in itertools.product(range(17), range(17), range(17)):
        assert math.isclose(packet_cap(halves, a, b, h), packet_cap(halves, b, a, h),
                            rel_tol=1e-14, abs_tol=1e-100)
        checks += 1
    return {"AES_GF256_columns_checked": columns,
            "fixed_mask_shared_point_products": 16**3,
            "exact_Fraction_alpha1_half_table": "passed",
            "float_alpha0523_half_table": "passed",
            "quiet_moment_endpoint_checks": quiet_checks,
            "half_exchange_checks": checks}


def update_extinction(args, polynomial, paired, filtered, middle):
    np = paired.np
    shapes = np.array(sorted(set(tuple(sorted((a, a, b, b)))
                      for a, b in itertools.product(range(17), repeat=2))), dtype=np.int16)
    index = {tuple(shape): i for i, shape in enumerate(shapes)}
    # Frozen joint-live and block-support refinements leave j=0 unchanged.
    # Recompute only its old uniform cap to take a per-(a,b) minimum.
    old = None
    for gamma in (0, .5, 1):
        local, _ = filtered.local_tables(shapes, args.z, args.alpha, gamma, args.activity_mark)
        col = local[:, 2, :, 0]
        old = col.copy() if old is None else np.minimum(old, col)
    degree = len(polynomial[0])-1
    rebuilt = np.zeros((degree+1, 17))
    halves = half_tables(args.z, args.alpha)
    changed = removed = 0
    gain = 0.0
    for a, b in itertools.product(range(17), repeat=2):
        paired.guard()
        if a+b > degree:
            continue
        row = index[tuple(sorted((a, a, b, b)))]
        coefficient = math.comb(16, a)*math.comb(16, b)
        for h in range(1, 17):
            cap = args.activity_mark*packet_cap(halves, a, b, h)
            before = old[row, h]
            after = min(before, cap)
            if after < before:
                changed += 1
                if after == 0:
                    removed += 1
                else:
                    gain = max(gain, math.log2(before/after))
            rebuilt[a+b, h] += coefficient*after
    quiet, nonquiet = polynomial
    updated = nonquiet.copy()
    updated[:, 1:, 0] = middle.to_log(rebuilt[:, 1:])
    assert np.array_equal(updated[:, :, 1:], nonquiet[:, :, 1:])
    assert np.array_equal(updated[:, 0], nonquiet[:, 0])
    return (quiet, updated), {"changed_ordered_half_entries": changed,
        "removed_ordered_half_entries": removed, "largest_finite_log2_gain": gain,
        "quiet00_unchanged": True, "live_destinations_unchanged": True,
        "terminal_unchanged": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--self-test", action="store_true")
    mode.add_argument("--run", action="store_true")
    options = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    start = time.monotonic()
    block = load_block()
    paired = block.load_paired(start+48)
    tests = {**self_test(), **paired.self_test()}
    if options.self_test:
        print(json.dumps({"tests": tests, "seconds": time.monotonic()-start,
                          "generated_bytes": 0}))
        return
    args = SimpleNamespace(r=256, z=.85, alpha=.523, activity_mark=1.0)
    filtered, middle = paired.imports()
    a, b, prior_effect, array_bytes = block.local_polynomials(args, paired, filtered, middle)
    a, effect = update_extinction(args, a, paired, filtered, middle)
    for _ in range(7):
        b, a = (paired.positive_convolution(a, b, args.r//2),
                paired.positive_convolution(a, a, args.r//2))
    denominator = math.log(math.comb(4096, args.r//2))
    region = paired.expand_coefficient(a, args.r//2, denominator)
    tail = paired.expand_coefficient(b, args.r//2, denominator)
    result = filtered.bound(args, region, tail, middle)
    regional = {"quiet00": float(region[0, 0, 0]/paired.LN2),
                "nonquiet00": float(region[2, 0, 0]/paired.LN2),
                "nonquiet16to16": float(region[2, 16, 16]/paired.LN2),
                "terminal_quiet0": float(tail[0, 0, 0]/paired.LN2)}
    print(json.dumps({"scope": "FLOAT pair-only REGIONAL singleton diagnostic",
        "whole_code_certificate": False, "actual_support_class_bound": False,
        "repeated_marginal_whole_chain_is_artificial": True,
        "parameters": vars(args), "singleton_effect": effect, "prior_effect": prior_effect,
        "regional_log2": regional, "artificial_repeated_marginal_score": result,
        "artificial_restart_floor_log2": result["source_log2_prefactor"]+256*regional["nonquiet00"],
        "tests": tests, "initial_array_bytes": array_bytes,
        "seconds": time.monotonic()-start, "generated_bytes": 0,
        "block_source_sha256": BLOCK_SHA256}, indent=2))


if __name__ == "__main__":
    main()
