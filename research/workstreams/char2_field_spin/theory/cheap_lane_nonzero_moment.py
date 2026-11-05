"""One FLOAT diagnostic using nonzero RS block moments, not a certificate.

For an incoming block mask of size h>0, n=16-ell output positions lie
outside the allowed inputs. Its nonzero auxiliary output moment is at most

 M_h(n,z)=sum_{v=0}^{min(3,n)} C(n,v)z^(n-v)(1-z)^v
          q^-h (q^(h-max(0,h-4+v))-1), q=256.

Requested sets of at least four zeros force the block state to zero.
For smaller sets, RS rank gives the displayed upper count of nonzero
solutions. This is exact for h=4. For h=1 take the minimum with
(1-q^-1)z^max(0,n-2); for h=2,3 use (1-q^-h)z^max(0,n-3).
Set M_0=1. No additional output floor or ((q-1)/q)^b factor is applied.

For four blocks, maximize the product of these moments over their state
sizes, tracking total h and the number b of active blocks. Multiply by
g(z)^ell and q/(q-b-1) for a live outgoing mask with inactive rows, or
q/(q-b) for the output-only terminal. Inactive outgoing equations cost
255^(-(16-j)). For j=16 there are none: the output-only cap suffices.
The support comparison and the joint D1 calculation justify these bounds.
Maxima precede alpha; destination multiplicity and activity stay outside.

Compute all old gamma/joint caps with their original terminals first.
Then minimize local caps, build positive pair-only polynomials, and apply
the frozen singleton extinction refinement. Never feed a modified terminal
back into the old joint cap. The frozen degree-convolution engine is reused.

The repeated-marginal chain score is artificial: a proper nonempty support
cannot be fully paired in both shifted matchings. It is not a class bound,
even if negative. --oracle-no-extinction additionally deletes regular
live-to-zero terms and is NOT a valid probability upper bound.

Fixed r256,z.85,alpha.523 only. Use python -B --self-test or --run.
Guards:20GiB free,48seconds,<256MiB arrays. No data, caches, or old edits.
"""

import argparse
from collections import Counter
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
SINGLETON_SHA256 = "992e6b7bd02e4785b9bfc4ae61437625afabb28a430e085fe32d60a3ab4a3453"


def load_singleton():
    path = Path(__file__).with_name("cheap_lane_singleton_extinction.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != SINGLETON_SHA256:
        raise RuntimeError("frozen singleton helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_singleton", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def moment(q, h, n, z):
    if h == 0:
        return 1
    value = sum(math.comb(n, v)*z**(n-v)*(1-z)**v
                *Fraction(q**(h-max(0, h-4+v))-1, q**h)
                for v in range(min(3, n)+1))
    if h <= 3:
        floor = max(0, n-(2 if h == 1 else 3))
        value = min(value, (1-Fraction(1, q**h))*z**floor)
    return value


def gf8_mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 8:
            a ^= 11
    return out


def moment_test():
    # RS[8,4] over GF8, all16 masks and all256 puncturings. The auxiliary
    # state is uniform on the mask's coordinate space, excluding only zero.
    q, z = 8, Fraction(17, 20)
    rows = [(1, a, gf8_mul(a, a), gf8_mul(gf8_mul(a, a), a)) for a in range(q)]
    powers = [z**w for w in range(q+1)]
    comparisons = equalities = words = 0
    for mask in range(16):
        coordinates = [i for i in range(4) if mask >> i & 1]
        h = len(coordinates)
        if h == 0:
            continue
        histogram = Counter()
        for u in itertools.product(range(q), repeat=h):
            if not any(u):
                continue
            support = 0
            for a, row in enumerate(rows):
                value = 0
                for coordinate, coefficient in zip(coordinates, u):
                    value ^= gf8_mul(row[coordinate], coefficient)
                support |= int(value != 0) << a
            histogram[support] += 1
            words += 1
        assert sum(histogram.values()) == q**h-1
        for outside in range(1 << q):
            weights = Counter()
            for support, count in histogram.items():
                weights[(support & outside).bit_count()] += count
            exact = sum(count*powers[weight] for weight, count in weights.items())/q**h
            cap = moment(q, h, outside.bit_count(), z)
            assert exact <= cap
            if h == 4:
                assert exact == cap
                equalities += 1
            comparisons += 1
    # Endpoint formulas are exact for a full-dimensional block.
    for n in range(9):
        assert moment(q, 4, n, Fraction(1)) == 1-Fraction(1, q**4)
        assert moment(q, 4, n, Fraction(0)) == Fraction(q**max(0, 4-n)-1, q**4)
    return {"GF8_nonzero_words": words, "mask_puncture_comparisons": comparisons,
            "full_dimension_exact_equalities": equalities,
            "positive_Fraction_MDS_moment": "passed"}


def local_polynomials(args, paired, filtered, middle):
    np = paired.np
    shapes = np.array(sorted(set(tuple(sorted((a, a, b, b)))
                      for a, b in itertools.product(range(17), repeat=2))), dtype=np.int16)
    index = {tuple(shape): i for i, shape in enumerate(shapes)}
    local = terminal = None
    for gamma in (0, .5, 1):
        candidate, final = filtered.local_tables(shapes, args.z, args.alpha,
                                                  gamma, args.activity_mark)
        for j in range(1, 17):
            factor = min(1.0, 256/255*255.0**(-(16-j)))
            candidate[:, :, :, j] = np.minimum(candidate[:, :, :, j],
                final[:, :, :, 0]*factor**args.alpha*math.comb(16, j))
        local = candidate if local is None else np.minimum(local, candidate)
        terminal = final if terminal is None else np.minimum(terminal, final)
    # Max-product DP tracks active block count for the support penalty.
    masses = np.array([[float(moment(256, h, 16-ell, args.z)) for h in range(5)]
                       for ell in range(17)])
    dp = np.zeros((len(shapes), 17, 5))
    dp[:, 0, 0] = 1.0
    for block in range(4):
        nxt = np.zeros_like(dp)
        for h in range(4*block+1):
            for b in range(block+1):
                for hb in range(5):
                    nb = b+int(hb > 0)
                    nxt[:, h+hb, nb] = np.maximum(nxt[:, h+hb, nb],
                        dp[:, h, b]*masses[shapes[:, block], hb])
        dp = nxt
    original, original_terminal = local.copy(), terminal.copy()
    g = args.z+(1-args.z)/256
    for row, shape in enumerate(shapes):
        paired.guard()
        ell = int(shape.sum())
        floors = sorted(max(0, 13-int(x)) for x in shape)
        for h in range(17):
            old_live = old_terminal = 0.0
            new_live = new_terminal = 0.0
            for b in range((h+3)//4, min(h, 4)+1):
                common = (255/256)**b*g**ell*args.z**sum(floors[:b])
                old_live = max(old_live, 256/(255-b)*common)
                old_terminal = max(old_terminal, 256/(256-b)*common)
                new_live = max(new_live, 256/(255-b)*g**ell*dp[row, h, b])
                new_terminal = max(new_terminal, 256/(256-b)*g**ell*dp[row, h, b])
            flag, activity = (0, 1.0) if h == 0 else (2, args.activity_mark)
            for j in range(1, 17):
                new = new_terminal if j == 16 else new_live
                mass = min(old_live, new)*255.0**(-(16-j))
                cap = activity*math.comb(16, j)*mass**args.alpha
                local[row, flag, h, j] = min(local[row, flag, h, j], cap)
            cap = activity*min(old_terminal, new_terminal)**args.alpha
            terminal[row, flag, h, 0] = min(terminal[row, flag, h, 0], cap)
    assert np.array_equal(local[:, :, :, 0], original[:, :, :, 0])
    changed = (local > 0) & (local < original)
    effect = {"improved_local_entries": int(changed.sum()),
              "improved_terminal_entries": int((terminal < original_terminal).sum()),
              "maximum_local_log2_gain": float(np.log2(original[changed]/local[changed]).max())
              if changed.any() else 0.0}
    degree = min(32, args.r//2)
    p = np.zeros((degree+1, 4, 17, 17))
    t = np.zeros((degree+1, 4, 17, 1))
    for a, b in itertools.product(range(17), repeat=2):
        if a+b <= degree:
            coefficient = math.comb(16, a)*math.comb(16, b)
            row = index[tuple(sorted((a, a, b, b)))]
            p[a+b] += coefficient*local[row]
            t[a+b] += coefficient*terminal[row]
    assert not np.count_nonzero(p[:, 0, 1:])
    assert not np.count_nonzero(t[:, 0, 1:])
    return ((middle.to_log(p[:, 0, 0]), middle.to_log(p[:, 2])),
            (middle.to_log(t[:, 0, 0]), middle.to_log(t[:, 2])), effect,
            int(local.nbytes+terminal.nbytes+dp.nbytes+p.nbytes+t.nbytes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--self-test", action="store_true")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--oracle-no-extinction", action="store_true")
    options = parser.parse_args()
    if options.oracle_no_extinction and not options.run:
        raise ValueError("oracle deletion requires --run")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    start = time.monotonic()
    singleton = load_singleton()
    block = singleton.load_block()
    paired = block.load_paired(start+48)
    tests = {**moment_test(), **singleton.self_test(), **paired.self_test()}
    paired.guard()
    if options.self_test:
        print(json.dumps({"tests": tests, "seconds": time.monotonic()-start,
                          "generated_bytes": 0}))
        return
    args = SimpleNamespace(r=256, z=.85, alpha=.523, activity_mark=1.0)
    filtered, middle = paired.imports()
    a, b, effect, array_bytes = local_polynomials(args, paired, filtered, middle)
    a, extinction = singleton.update_extinction(args, a, paired, filtered, middle)
    if options.oracle_no_extinction:
        a[1][:, 1:, 0] = float("-inf")
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
    print(json.dumps({"scope": "FLOAT artificial repeated-marginal diagnostic",
        "whole_code_certificate": False, "actual_support_class_bound": False,
        "oracle_deletion_NOT_a_probability_upper_bound": options.oracle_no_extinction,
        "parameters": vars(args), "moment_effect": effect, "singleton_effect": extinction,
        "regional_log2": regional, "artificial_repeated_marginal_score": result,
        "artificial_restart_floor_log2": result["source_log2_prefactor"]+256*regional["nonquiet00"],
        "tests": tests, "initial_array_bytes": array_bytes,
        "seconds": time.monotonic()-start, "generated_bytes": 0,
        "singleton_source_sha256": SINGLETON_SHA256}, indent=2))


if __name__ == "__main__":
    main()
