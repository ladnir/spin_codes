"""Carry-safe exact-coefficient convolution for the directed cheap-lane cone.

Only polynomial multiplication changes. With full-degree cap F, encode
(f,s) as f+(2F+1)s. Products have f<=2F, so no carry aliases a forbidden
full degree into a singleton degree. After each positive product discard
f>F, s>S, and2f+s>256. All retained source profiles have2f+s=256.

Default: reproduce the complete1..16-run source union. Balanced modes
screen only f0=f1=128-c/2,s0=s1=c. Until a faster exact profile counter is
available, their source count is the total c-run count, an UPPER bound for
the balanced subset. It does NOT certify other c-run profiles. Exact
dyadic scales cancel in extraction; there is no Cauchy relaxation.

Frozen local tables, pre-route alpha, interval tests, terminal and source
factor are reused unchanged.48-second cooperative guard,20GiB free,
estimated256MiB peak, stdout only. No encoder change or generated files.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
CORE_PIN = "7ccf0cdab823ef15fdd5edd6156d4b49d9a000f2ae28c37d5d43f4b6ea341f01"
TEST_PIN = "83e00b6f6d148c8058428295b4132d401df08cfccd20cadcda1c1ba82df7b478"
COUNT_PIN = "25e0d9221a62ffdb79d1d625e6e858ffbb8c0b671078d2bde5d2c6fb27e43ed3"


def load(filename, digest, name):
    path = Path(__file__).with_name(filename)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("frozen source hash mismatch: "+filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def packed_class(core):
    class PackedVerifier(core.Verifier):
        packed_full = 128

        def multiply(self, left, right, maximum=128):
            assert len(left[0]) == len(right)
            full = min(maximum, self.packed_full)
            radix = 2*full+1

            def encode(entry):
                data = [self.arb(0)]*(radix*(self.single+1))
                for s, polynomial in enumerate(entry):
                    limit = min(full, len(polynomial)-1, (256-s)//2)
                    for f in range(limit+1):
                        data[f+radix*s] = polynomial[f]
                return self.poly(data)

            ll = [[encode(entry) for entry in row] for row in left]
            rr = [[encode(entry) for entry in row] for row in right]
            result = []
            for row in ll:
                output = []
                for j in range(len(rr[0])):
                    self.guard()
                    value = self.poly([])
                    for k, polynomial in enumerate(row):
                        value += (polynomial*rr[k][j]).truncate(radix*(self.single+1))
                    entry = []
                    for s in range(self.single+1):
                        limit = min(full, (256-s)//2)
                        entry.append(self.poly([self.positive(value[f+radix*s])
                                                for f in range(limit+1)]))
                    output.append(entry)
                result.append(output)
            return result
    return PackedVerifier


def tiny_tests(core, packed, exact_point, deadline):
    """Packed/unpacked matrices independently enclose exact Fraction products."""
    raw = [[{(0, 0): F(1, 3), (1, 0): F(2, 5), (0, 1): F(3, 7)},
            {(1, 0): F(4, 9), (0, 1): F(1, 11)}],
           [{(0, 0): F(2, 7), (1, 1): F(3, 13)},
            {(0, 0): F(5, 17), (1, 0): F(2, 19)}]]
    checks = 0
    for maximum in (0, 1, 2):
        v = packed(2, deadline); v.packed_full = maximum
        ordinary = core.Verifier(2, deadline)
        encoded = []
        for row in raw:
            encoded.append([[v.poly([v.exact(item.get((f, s), F(0))) for f in range(2)])
                             for s in range(3)] for item in row])
        # Initial truncation is required in the exact reference too.
        reference = [[{key: x for key, x in item.items() if key[0] <= maximum}
                      for item in row] for row in raw]
        trimmed = [[[v.poly([v.exact(item.get((f, s), F(0))) for f in range(maximum+1)])
                     for s in range(3)] for item in row] for row in reference]
        got = v.multiply(encoded, encoded)
        old = ordinary.multiply(trimmed, trimmed, maximum)
        for i in range(2):
            for j in range(2):
                target = defaultdict(F)
                for k in range(2):
                    for (f, s), x in reference[i][k].items():
                        for (ff, ss), y in reference[k][j].items():
                            if f+ff <= maximum and s+ss <= 2:
                                target[f+ff, s+ss] += x*y
                for s in range(3):
                    for f in range(maximum+1):
                        exact_point(got[i][j][s][f], target[f, s], v.fmpq)
                        exact_point(old[i][j][s][f], target[f, s], v.fmpq)
                        checks += 1
    v = packed(2, deadline); v.packed_full = 2
    assert v.toy_test() == 6
    return checks


def ordered_chain(v, r0, r1, terminal1):
    value = v.matmul(v.matpow(v.matmul(r0, r1), 127), r0)
    return v.positive(v.matmul(value, terminal1)[0][0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--balanced-runs", type=int, choices=(32, 64, 128, 256))
    parser.add_argument("--precision", type=int, choices=(192, 256, 384, 512), default=256)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+48
    test_adapter = load("cheap_lane_cone_interval_tests.py", TEST_PIN, "packed_cone_tests")
    assert test_adapter.PIN == CORE_PIN
    core = test_adapter.load_verifier()
    counts = load("cheap_lane_run_profile_counts.py", COUNT_PIN, "packed_cone_counts")
    counts.DEADLINE = deadline
    from flint import arb, ctx
    ctx.prec = args.precision
    packed = packed_class(core)
    exact_tests = tiny_tests(core, packed, test_adapter.exact_point, deadline)
    if args.balanced_runs is None:
        single, full = 32, 128
        indices = list(range(17))
    else:
        single = args.balanced_runs
        full = 128-single//2
        indices = [single//2]
        # Only in-memory constants on this freshly loaded module change.
        scales = {32: (F(1, 32), F(1, 128)), 64: (F(1, 32), F(1, 64)),
                  128: (F(1, 64), F(1, 32)), 256: (F(1, 32), F(1, 16))}
        core.USCALE, core.VSCALE = scales[single]
    v = packed(single, deadline); v.packed_full = full
    records, route_checks = core.route_records(single)
    # Count checks above include all relevant local degrees; dropping local
    # monomials that cannot reach any requested profile is then exact.
    records = {key: weight for key, weight in records.items()
               if key[0] <= full and 2*key[0]+key[1] <= 256}
    v.guard()
    shapes = sorted({key[2] for key in records})
    caps = v.local_caps(shapes)
    regular, terminal = v.regions(*v.polynomials(records, caps))
    rr = {a: v.extract(regular, 128-a, 2*a) for a in indices}
    tt = {a: v.extract(terminal, 128-a, 2*a) for a in indices}
    source = (v.exact(core.ETA)*v.exact(core.BETA)**256
              *v.exact(core.Z)**(-core.LIMIT))**v.exact(core.ALPHA)
    logtwo = arb(2).log()
    rows = []
    total = arb(0)
    indexed_sources = 0
    tests = {"exact_packed_unpacked_coefficients": exact_tests,
             "exact_route_normalizers": route_checks, "terminal_updates": 32767}
    if args.balanced_runs is None:
        for runs in range(1, 17):
            subtotal = arb(0); number_sum = 0
            for a in range(runs+1):
                b = runs-a
                number = counts.profile_count(4096, 256, 128-a, 128-b)
                subtotal += number*source*ordered_chain(v, rr[a], rr[b], tt[b])
                number_sum += number
                v.guard()
            assert number_sum == counts.run_total(4096, 256, runs)
            total += subtotal; indexed_sources += number_sum
            rows.append({"runs": runs, "profiles": runs+1, "indexed_sources": number_sum,
                         "log2_bound_interval": str(subtotal.log()/logtwo)})
        # Saved192-bit interval is only a replay check, not an input to the
        # proof. No large unpacked polynomial product is repeated here.
        control = arb("-133.814902217474350180257723338778012258971492919659374",
                      "3.89e-52")
        assert (total.log()/logtwo).overlaps(control)
        tests["unpacked16_control"] = "overlaps saved192-bit directed interval"
        count_kind = "exact indexed counts for all profiles with1..16 runs"
        scope = "all state paths; all size256 source supports with1..16 cyclic runs"
    else:
        runs = args.balanced_runs; a = runs//2
        conditional = source*ordered_chain(v, rr[a], rr[a], tt[a])
        number = counts.run_total(4096, 256, runs)
        total = number*conditional; indexed_sources = number
        rows.append({"runs": runs, "a": a, "b": a,
                     "full_per_matching": 128-a, "single_per_matching": 2*a,
                     "conditional_log2_bound_interval": str(conditional.log()/logtwo),
                     "source_count_upper": number,
                     "log2_bound_interval": str(total.log()/logtwo)})
        count_kind = "total c-run count; an upper bound for the balanced profile only"
        scope = "all state paths; one balanced size256 source profile only"
    v.guard()
    passed = bool(total < arb(2)**(-40))
    print(json.dumps({"scope": scope, "directed": True, "whole_code_certificate": False,
        "source_count_kind": count_kind, "source_count_or_upper": indexed_sources,
        "parameters": {"z": str(core.Z), "alpha": str(core.ALPHA),
            "full_degree": full, "singleton_degree": single, "packing_radix": 2*full+1,
            "full_scale": str(core.USCALE), "singleton_scale": str(core.VSCALE)},
        "precision_bits": args.precision, "rows": rows,
        "log2_bound_interval": str(total.log()/logtwo),
        "proved_below_2_pow_minus_40": passed,
        "proved_below_2_pow_minus_100": bool(total < arb(2)**(-100)),
        "tests": tests, "route_records": len(records), "occupancy_shapes": len(shapes),
        "nonnegative_intersections": v.intersections, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 256*(1 << 20), "generated_bytes": 0,
        "core_pin": CORE_PIN, "test_pin": TEST_PIN, "count_pin": COUNT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))
    if not passed:
        raise SystemExit("directed interval does not prove the fixed40-bit target")


if __name__ == "__main__":
    main()
