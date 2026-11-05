"""Directed exact all-singleton profile, not a complete source-size cover.

Fix r in {512,1024,2048,4096}. Both alternating matchings have full-pair
degree zero and singleton degree r. Set the full variable to zero exactly;
this removes only coefficients irrelevant to this face. Build the positive
singleton polynomial from the actual two16-slot halves and fair flips.
The frozen local caps are raised to alpha BEFORE route averaging.

Seven positive univariate doublings give P^128 and P^127 T, truncated at
degree r. Divide the selected coefficients by C(4096,r)*v^r, then use
e0 R^255 V. The terminal includes input degrees but no final refresh.
The indexed source count is8192/r*C(8192-r-1,r-1), crosschecked against
the frozen general profile counter. Apply the own-r source factor once.

This coefficient is independent of the positive singleton scale v. No
residue aliases, Cauchy fallback, new encoder, or full-field refresh law
are used. All bounds use Arb; only timing and memory metadata use FLOAT.
180s cooperative limit,256MiB process guard,20GiB free, stdout only.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
ANCHOR_PIN = "319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
COUNTS_PIN = "4ce61f179a8bf41e5d769e955a9efeb437e46af81f57b0d784c4cfc57958d8e8"
CORE_PIN = "7ccf0cdab823ef15fdd5edd6156d4b49d9a000f2ae28c37d5d43f4b6ea341f01"
ADAPTER_PIN = "83e00b6f6d148c8058428295b4132d401df08cfccd20cadcda1c1ba82df7b478"
SINGLE_CORRECTION = F(1140210, 880327)
MEMORY_LIMIT, SECONDS = 256*(1 << 20), 180


def load(name, pin, label):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen helper hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(label, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def local_records():
    """Integer orientation numerators over a common2^32 denominator."""
    half = [(s, left, s-left,
             comb(16, s)*comb(s, left)*(1 << (16-s)))
            for s in range(17) for left in range(s+1)]
    records = defaultdict(int)
    for s, a, b, mass in half:
        for ss, c, d, other in half:
            records[s+ss, tuple(sorted((a, b, c, d)))] += mass*other
    totals = [0]*33
    for (s, _), mass in records.items():
        totals[s] += mass
    for s, mass in enumerate(totals):
        assert mass == comb(32, s)*(1 << 32)
    return dict(records)


def multiply(v, left, right, maximum):
    """Carry-free univariate matrix product, with a positive downset."""
    assert len(left[0]) == len(right)
    out = []
    for row in left:
        result = []
        for j in range(len(right[0])):
            v.guard()
            value = sum((a*right[k][j] for k, a in enumerate(row)), v.poly([]))
            value = value.truncate(maximum+1)
            result.append(v.poly([v.positive(value[d]) for d in range(len(value))]))
        out.append(result)
    return out


def exact_tests(v, adapter):
    """Compare truncating products to an independent Fraction convolution."""
    p = [[[F(1, 3), F(2, 7)], [F(3, 11), F(4, 13), F(1, 5)]],
         [[F(2, 9), F(5, 17)], [F(7, 19), F(1, 4)]]]
    t = [[[F(1, 5), F(2, 3)]], [[F(3, 7), F(5, 11)]]]

    def exact_product(left, right):
        out = []
        for row in left:
            rr = []
            for j in range(len(right[0])):
                length = max(len(a)+len(right[k][j])-1 for k, a in enumerate(row))
                value = [F(0)]*length
                for k, a in enumerate(row):
                    for i, x in enumerate(a):
                        for d, y in enumerate(right[k][j]):
                            value[i+d] += x*y
                rr.append(value)
            out.append(rr)
        return out

    def wrapped(table):
        return [[v.poly([v.exact(x) for x in item]) for item in row] for row in table]

    er, et = p, t
    for _ in range(2):
        et = exact_product(er, et)
        er = exact_product(er, er)
    checks = 0
    for maximum in (0, 1, 3, 6, 9):
        rr, tt = wrapped(p), wrapped(t)
        for _ in range(2):
            tt = multiply(v, rr, tt, maximum)
            rr = multiply(v, rr, rr, maximum)
        for actual, expected in ((rr, er), (tt, et)):
            for i, row in enumerate(actual):
                for j, entry in enumerate(row):
                    for d in range(maximum+1):
                        point = expected[i][j][d] if d < len(expected[i][j]) else F(0)
                        adapter.exact_point(entry[d], point, v.fmpq)
                        checks += 1
    # Exact target coefficients do not change under a positive variable scale.
    unscaled = wrapped(p)
    for scale in (F(2, 3), F(7, 5)):
        scaled = [[v.poly([v.exact(x*scale**d) for d, x in enumerate(item)])
                   for item in row] for row in p]
        aa = multiply(v, scaled, scaled, 3)
        bb = multiply(v, unscaled, unscaled, 3)
        for i in range(2):
            for j in range(2):
                for d in range(4):
                    assert (aa[i][j][d]/v.exact(scale**d)).overlaps(bb[i][j][d])
                    checks += 1
    # The closed-form count includes both orientations and does not quotient rotations.
    for n in (6, 8):
        for r in range(1, n//2+1):
            count = sum(sum(word) == r and all(not(word[i] and word[(i+1) % n])
                        for i in range(n)) for word in product((0, 1), repeat=n))
            numerator = n*comb(n-r-1, r-1)
            assert numerator % r == 0 and count == numerator//r
            checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight", type=int, choices=(512, 1024, 2048, 4096), default=1024)
    parser.add_argument("--z", type=F, default=F(304397809, 625000000))
    parser.add_argument("--alpha", type=F, default=F(119, 1280))
    parser.add_argument("--v", type=F)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=256)
    args = parser.parse_args()
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    r = args.weight
    scale = args.v if args.v is not None else (
        SINGLE_CORRECTION*F(r, 4096-r) if r < 4096 else F(1))
    if scale <= 0:
        parser.error("singleton scale must be positive")
    if max(x.bit_length() for a in (args.z, args.alpha, scale)
           for x in (a.numerator, a.denominator)) > 64:
        parser.error("rational witnesses limited to64 bits")
    free = shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+SECONDS
    anchor = load("cheap_lane_cone_cauchy_anchors.py", ANCHOR_PIN, "singleton_anchor")
    assert anchor.CORE_PIN == CORE_PIN and anchor.ADAPTER_PIN == ADAPTER_PIN
    counts = load("cheap_lane_fast_profile_counts.py", COUNTS_PIN, "singleton_counts")
    counts.DEADLINE = deadline
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("180-second singleton-face budget exceeded")
        if peak_memory() > MEMORY_LIMIT:
            raise MemoryError("256MiB singleton-face budget exceeded")

    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z, core.ALPHA = args.z, args.alpha
    from flint import arb, ctx
    ctx.prec = args.precision

    class SingletonVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = SingletonVerifier(2, deadline)
    tests = exact_tests(verifier, adapter)
    records = local_records()
    shapes = sorted({occ for _, occ in records})
    caps = verifier.local_caps(shapes)
    powered = {occ: tuple(verifier.positive(x)**verifier.exact(args.alpha) for x in row)
               for occ, row in caps.items()}
    entries = [[arb(0) for _ in range(33)] for _ in range(6)]
    for (s, occ), numerator in records.items():
        guard()
        weight = verifier.exact(F(numerator, 1 << 32)*scale**s)
        for k, value in enumerate(powered[occ]):
            entries[k][s] += weight*value
    polys = [verifier.poly(row) for row in entries]
    rr, tt = [[polys[0], polys[1]], [polys[2], polys[3]]], [[polys[4]], [polys[5]]]
    del caps, powered, entries, polys
    print(json.dumps({"phase": "singleton_local_polynomials_ready", "weight": r,
        "shapes": len(shapes), "records": len(records), "seconds": time.monotonic()-started}), flush=True)
    packets, updates = 1, 0
    for stage in range(7):
        tt = multiply(verifier, rr, tt, r)
        rr = multiply(verifier, rr, rr, r)
        packets, updates = 2*packets, packets+updates
        print(json.dumps({"phase": "univariate_doubling", "stage": stage+1,
            "packets": packets, "terminal_updates": updates,
            "seconds": time.monotonic()-started}), flush=True)
    assert packets == 128 and updates == 127
    denominator = arb(comb(4096, r))*verifier.exact(scale)**r
    assert denominator > 0
    regular = [[verifier.positive(item[r])/denominator for item in row] for row in rr]
    terminal = [[verifier.positive(item[r])/denominator for item in row] for row in tt]
    chain = verifier.positive(verifier.matmul(verifier.matpow(regular, 255), terminal)[0][0])
    numerator = 8192*comb(8192-r-1, r-1)
    assert numerator % r == 0
    number = numerator//r
    assert number == counts.profile_count(4096, r, 0, 0)
    h, c_h, eta = anchor.source_factor(r)
    source = (verifier.exact(eta)*verifier.exact(core.BETA)**r
              *verifier.exact(args.z)**(-core.LIMIT))**verifier.exact(args.alpha)
    bound = number*source*chain
    assert bound.is_finite() and bound > 0
    guard()
    logtwo = arb(2).log()
    print(json.dumps({"scope": "one ordered all-singleton profile at the requested actual source weight",
        "directed": True, "exact_coefficients": True, "whole_weight_certificate": False,
        "whole_code_certificate": False, "source_weight": r,
        "ordered_profiles": [[0, r], [0, r]], "full_variable": "exact zero face",
        "parameters": {"z": str(args.z), "alpha": str(args.alpha), "singleton_scale": str(scale),
                       "distance_threshold": core.LIMIT},
        "log2_bound_interval": str(bound.log()/logtwo),
        "below_2_pow_minus_40": bool(bound < arb(2)**(-40)),
        "log2_conditional_chain_interval": str(chain.log()/logtwo),
        "log2_source_factor_interval": str(source.log()/logtwo),
        "log2_indexed_support_count_interval": str(arb(number).log()/logtwo),
        "indexed_support_count_sha256": hashlib.sha256(str(number).encode()).hexdigest(),
        "source_factor": {"h": h, "c_h": str(c_h), "eta": str(eta)},
        "tests": {"tiny_exact_coefficients_and_counts": tests, "local_degree_normalizers": 33,
                  "selected_indexed_count_matches_general_counter": True,
                  "regional_packets": packets, "regional_terminal_updates": updates,
                  "total_updates": 255*128+updates, "terminal_emissions": 1},
        "precision_bits": args.precision, "nonnegative_intersections": verifier.intersections,
        "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory(),
        "free_bytes_preflight": free, "generated_bytes": 0,
        "dependencies": {"anchor": ANCHOR_PIN, "counter": COUNTS_PIN,
                         "adapter": ADAPTER_PIN, "core": CORE_PIN},
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
