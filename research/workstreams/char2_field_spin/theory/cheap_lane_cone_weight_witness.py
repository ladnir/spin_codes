"""A bounded directed cone witness for one requested source weight.

Reuse the frozen independent scalar caps and exact source counts. Only
the rational witness z,alpha and positive coefficient truncation change.
For r<=256, the initial support factor is
 eta_r=256^2/[255^2(257-r)].
For257<=r<=384, explicit --allow-missing-one uses the proved h=1 factor
 [1+r*256^256/beta]/[255*(1-r/512)]. Pad missing input columns to the
ACTUAL support T before applying the same actual-T cone caps. Neither
local filters nor routing profiles are evaluated on a shadow support.
The common factor is (support_factor*beta^r*z^-734003)^alpha.

Carry-safe packing has radix2F+1 and retains2f+s<=r, not the frozen256
cutoff. Deterministic rectangles cover all(f,r-2f),0<=f<=floor(r/2).
Exact dyadic scales cancel before composing regions. All ordered indexed
source profiles are counted; the final packet has no update. No weight
interpolation, whole-code claim, construction change, or automatic scan.

Default z=1-r/1600 and alpha=119/1280. Positive bounds are failed fixed
witnesses, not verifier failures.180-second cooperative guard,20GiB free,
estimated512MiB peak, bounded stdout, no generated files.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
COVER_PIN = "84ac66e9e3b49b854ce10adf2b2f244d74a530524252736d35d1e72cdfe6e7fb"
CONTROL_CENTER = "-93.5574055068572906550312346546009957059259358406723846743061455591754432859448382230484510871964770748120362074"


def load_cover():
    path = Path(__file__).with_name("cheap_lane_cone_profile_cover.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != COVER_PIN:
        raise RuntimeError("frozen profile-cover hash mismatch")
    spec = importlib.util.spec_from_file_location("weight_witness_cover", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rectangles(weight):
    """At most four balanced full-count intervals; no numerical tuning."""
    size = weight//2+1
    bands = min(4, (size+32)//33)
    result = []
    for i in reversed(range(bands)):
        low = i*size//bands
        high = (i+1)*size//bands-1
        middle = (low+high+1)//2
        result.append((low, high, F(max(1, middle), 4096),
                       F(max(1, weight-2*middle), 4096)))
    return result


def make_verifier(core, weight):
    class WeightVerifier(core.Verifier):
        packed_full = weight//2

        def guard(self):
            if time.monotonic() > self.deadline:
                raise TimeoutError("180-second single-weight witness budget exceeded")

        def multiply(self, left, right, maximum=None):
            assert len(left[0]) == len(right)
            full = min(self.packed_full, weight//2 if maximum is None else maximum)
            radix = 2*full+1

            def encode(entry):
                data = [self.arb(0)]*(radix*(self.single+1))
                for s, polynomial in enumerate(entry):
                    limit = min(full, len(polynomial)-1, (weight-s)//2)
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
                    entries = []
                    for s in range(self.single+1):
                        limit = min(full, (weight-s)//2)
                        entries.append(self.poly([self.positive(value[f+radix*s])
                                                  for f in range(limit+1)]))
                    output.append(entries)
                result.append(output)
            return result
    return WeightVerifier


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight", type=int, required=True)
    parser.add_argument("--z", help="exact rational; default1-weight/1600")
    parser.add_argument("--alpha", default="119/1280", help="exact rational in(0,1]")
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    parser.add_argument("--allow-missing-one", action="store_true")
    args = parser.parse_args()
    if not 7 <= args.weight <= 384:
        parser.error("weight must be7..384")
    if args.weight > 256 and not args.allow_missing_one:
        parser.error("weights above256 require explicit --allow-missing-one")
    z = F(args.z) if args.z is not None else 1-F(args.weight, 1600)
    alpha = F(args.alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if max(x.bit_length() for x in (z.numerator, z.denominator,
                                    alpha.numerator, alpha.denominator)) > 64:
        parser.error("witness rational numerators/denominators must fit64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+180

    def progress(phase, **values):
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started, **values}), flush=True)

    cover = load_cover()
    adapter = cover.load("cheap_lane_cone_interval_tests.py", cover.TEST_PIN, "single_weight_tests")
    assert adapter.PIN == cover.CORE_PIN
    core = adapter.load_verifier()
    # Modify constants only on this fresh in-memory module; no frozen edit.
    core.Z, core.ALPHA, core.MAX_FULL = z, alpha, args.weight//2
    packed = cover.load("cheap_lane_cone_packed.py", cover.PACK_PIN, "single_weight_packed")
    assert packed.CORE_PIN == cover.CORE_PIN
    counts = cover.load("cheap_lane_fast_profile_counts.py", cover.COUNT_PIN, "single_weight_counts")
    counts.DEADLINE = deadline
    from flint import arb, ctx
    ctx.prec = args.precision
    WeightVerifier = make_verifier(core, args.weight)
    toy_checks = packed.tiny_tests(core, WeightVerifier, adapter.exact_point, deadline)
    # Exercise the actual variable cutoff, including full degrees above128.
    boundary = WeightVerifier(args.weight % 2+1, deadline)
    boundary.packed_full = args.weight//2
    single = args.weight % 2
    monomial = [boundary.poly([]) for _ in range(boundary.single+1)]
    monomial[single] = boundary.poly([0]*(args.weight//2)+[1])
    one_plus_single = [boundary.poly([1]), boundary.poly([1])]
    one_plus_single += [boundary.poly([]) for _ in range(boundary.single-1)]
    product = boundary.multiply([[monomial]], [[one_plus_single]])[0][0]
    adapter.exact_point(product[single][args.weight//2], F(1), boundary.fmpq)
    assert product[single+1][args.weight//2].is_zero()
    count_checks = counts.checks()
    setup = WeightVerifier(min(32, args.weight), deadline)
    records, route_checks = core.route_records(min(32, args.weight))
    # Local terms exceeding the total regional input budget cannot contribute.
    records = {key: value for key, value in records.items() if 2*key[0]+key[1] <= args.weight}
    shapes = sorted({key[2] for key in records})
    caps = setup.local_caps(shapes)
    exponent = setup.exact(alpha)
    powered = {occ: tuple(setup.positive(x)**exponent for x in cap)
               for occ, cap in caps.items()}
    del caps
    progress("local_caps_ready", source_weight=args.weight, records=len(records),
             occupancy_shapes=len(shapes))
    rr, tt = {}, {}
    band_reports = []
    intersections = 0
    for low, high, uscale, vscale in rectangles(args.weight):
        core.USCALE, core.VSCALE = uscale, vscale
        v = WeightVerifier(args.weight-2*low, deadline); v.packed_full = high
        regular, terminal = v.regions(*cover.local_polynomials(v, core, records, powered))
        for f in range(low, high+1):
            single = args.weight-2*f
            rr[f] = v.extract(regular, f, single)
            tt[f] = v.extract(terminal, f, single)
        intersections += v.intersections
        band_reports.append({"full_min": low, "full_max": high, "singleton_max": v.single,
                             "full_scale": str(uscale), "singleton_scale": str(vscale)})
        progress("band_ready", full_min=low, full_max=high, retained_profiles=len(rr))
        del regular, terminal
    assert set(rr) == set(tt) == set(range(args.weight//2+1))
    if args.weight <= 256:
        support = F(256**2, 255**2*(257-args.weight))
        support_details = {"comparison": "exact support", "eta_exact": str(support)}
    else:
        c1 = 1-F(args.weight, 512)
        assert c1 > 0
        support = (1+args.weight*F(256**256)/core.BETA)/(255*c1)
        support_details = {"comparison": "at most one missing group; pad to actual T first",
                           "c1_exact": str(c1), "maximum_missing_groups": 1,
                           "factor": "(1+r*256^256/beta)/(255*c1)"}
    source = (setup.exact(support)*setup.exact(core.BETA)**args.weight
              *setup.exact(z)**(-core.LIMIT))**exponent
    logtwo = arb(2).log()
    conditional = arb(0); total_sources = 0; chains = 0
    numbers = {}; run_numbers = {}; run_terms = {}
    full_max = args.weight//2
    for f0 in range(full_max+1):
        for f1 in range(full_max+1):
            if f0 <= f1:
                number = counts.profile_count(4096, args.weight, f0, f1)
                numbers[f0, f1] = number
            else:
                number = numbers[f1, f0]
            total_sources += number
            runs = args.weight-f0-f1
            run_numbers[runs] = run_numbers.get(runs, 0)+number
            if number:
                value = number*packed.ordered_chain(setup, rr[f0], rr[f1], tt[f1])
                conditional += value
                run_terms[runs] = run_terms.get(runs, arb(0))+value
                chains += 1
        setup.guard()
    assert total_sources == comb(8192, args.weight)
    assert all(run_numbers.get(c, 0) == counts.run_total(4096, args.weight, c)
               for c in range(args.weight+1))
    bound = source*conditional
    bits = bound.log()/logtwo
    control = "not the frozen r256 witness"
    if args.weight == 256 and z == F(21, 25) and alpha == F(119, 1280):
        assert bits.overlaps(arb(CONTROL_CENTER, "3.25e-110"))
        control = "overlaps frozen384-bit r256 interval"
    worst_run = max(run_terms, key=lambda c: run_terms[c].upper())
    setup.guard()
    print(json.dumps({"scope": "all state paths and all indexed supports at one source weight",
        "directed": True, "verifier_completed": True, "whole_code_certificate": False,
        "other_source_weights_covered": False,
        "parameters": {"source_weight": args.weight, "z": str(z), "alpha": str(alpha),
                       "distance_threshold": core.LIMIT, "bands": band_reports},
        "support_comparison": support_details, "precision_bits": args.precision,
        "log2_bound_interval": str(bits),
        "proved_below_2_pow_minus_40": bool(bound < arb(2)**(-40)),
        "proved_below_2_pow_minus_100": bool(bound < arb(2)**(-100)),
        "worst_run": {"runs": worst_run, "selection": "largest upper endpoint",
                      "log2_bound_interval": str((source*run_terms[worst_run]).log()/logtwo)},
        "indexed_sources": total_sources, "ordered_profile_chains": chains,
        "tests": {**count_checks, "exact_packed_unpacked_coefficients": toy_checks,
            "actual_weight_boundary_kept_and_overflow_dropped": "passed",
            "exact_route_normalizers_before_downset": route_checks,
            "source_count_sum": "C(8192,r)", "per_run_totals": "passed",
            "terminal_updates": 32767, "frozen_control": control},
        "nonnegative_intersections": intersections+setup.intersections,
        "seconds": time.monotonic()-started, "estimated_peak_bytes": 512*(1 << 20),
        "generated_bytes": 0, "cover_pin": COVER_PIN, "core_pin": cover.CORE_PIN,
        "test_pin": cover.TEST_PIN, "packed_pin": cover.PACK_PIN, "count_pin": cover.COUNT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
