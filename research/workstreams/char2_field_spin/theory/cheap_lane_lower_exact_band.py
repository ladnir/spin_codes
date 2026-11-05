"""Complete small-source bands with current pair caps and exact downsets.

The pinned local constructor supplies six positive coefficient upper bounds:
four entries of H=M*M and two entries of J=M*T. It uses the current actual-
support cone, residual-rank returns, and ordered concavity correction. This
driver does not change that operator or apply another local fractional power.

Each rectangle retains exact full and singleton degrees after every positive
integer product. Six updates T=R*T, R=R*R produce H^64 and H^63*J. Discarded
degrees cannot return to a retained target. Integer ceilings enlarge retained
coefficients; thus "exact" describes degree extraction, not unrounded values.
Positive rational markers cancel against each target's own denominator.

For every requested actual r and every f=0..floor(r/2), extract s=r-2*f and
divide by C(4096,f) C(4096-f,s) u^f v^s. The frozen join retains both terminal
orders and charges exact indexed run masses to their maximum conditional
chain. It supports odd r explicitly. Each actual r receives its own current
joint-missing eta_r and beta^r before the one source power alpha. The legacy
stronger exact-support eta for some r<=256 is not silently substituted.

Only a completed sum covers the requested band. A timeout or resource failure
produces no partial certificate. Numerical receipts are indexed separately;
source alone is not evidence. Import performs no numerical work. Cooperative
guards check300/600 seconds and512MiB peak process memory between operations;
they do not interrupt an in-progress FLINT call. A run requires20GiB free
space. Output is stdout only; no arrays or receipt files are generated.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import re
import shutil
import sys
import time

sys.dont_write_bytecode = True
LOCAL_PIN = "1ed1d59d9beed36b9a1bd1a2bc65af043efe342e95e2afd5047f90c696c04c47"
INTEGER_PIN = "9dc583ae54caed0a78f71e409cc37aa398713b419f65538370f10c3f876f40bf"
REGIONS_PIN = "ab74014cbd1e2e6334ab15acaf2fcbc9dc6d9bc97b19e5e6c495f0d2eb8a5852"
JOIN_PIN = "54f2ff06a1476ccf84fa4db1b757f1d5f8e4e814a1f0f733dc73ccf18b2ddabb"
ANCHOR_PIN = "319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
POINTS, SLOTS, MEMORY = 8192, 4096, 512*(1 << 20)
PACKED_LIMIT, MAX_TABLES = 280000, 64


def load(name, pin, label):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen source mismatch: "+name)
    spec = importlib.util.spec_from_file_location(label, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rational(text):
    """Bound parsing before constructing integers or Fractions."""
    if not re.fullmatch(r"[0-9]{1,20}(?:/[0-9]{1,20})?", text):
        raise argparse.ArgumentTypeError("use a positive integer ratio, at most20 digits per part")
    try:
        value = F(text)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError("invalid rational") from error
    if value <= 0 or max(value.numerator.bit_length(), value.denominator.bit_length()) > 64:
        raise argparse.ArgumentTypeError("positive rational with64-bit numerator/denominator required")
    return value


def rectangle_schedule(last):
    """Disjoint assigned full ranges; each table keeps a downward rectangle.

    The nominal four-way split is shortened when its packed polynomial
    would exceed280000 positions. This is below the engine's300000 limit.
    The upper domain may still exhaust the cooperative resource guards.
    """
    if type(last) is not int or not 225 <= last <= 1024:
        raise ValueError("last weight must be225..1024")
    maximum = last//2
    width = (maximum+4)//4
    result, lo = [], 0
    while lo <= maximum:
        single = last-2*lo
        allowed = (PACKED_LIMIT-single-1)//(2*single+1)
        hi = min(maximum, lo+width-1, allowed)
        if hi < lo:
            raise ValueError("one target exceeds the bounded exact rectangle layout")
        middle = (lo+hi)//2
        marker = (F(max(1, middle), SLOTS), F(max(1, last-2*middle), SLOTS))
        length = single*(2*hi+1)+hi+1
        assert length <= PACKED_LIMIT < 300000
        result.append((lo, hi, single, marker))
        lo = hi+1
        if len(result) > MAX_TABLES:
            raise ValueError("adaptive rectangle count exceeds64")
    return tuple(result)


def coefficient_hash(degrees):
    digest = hashlib.sha256()
    for key in sorted(degrees):
        digest.update((str(key)+":"+":".join(map(str, degrees[key]))+"\n").encode())
    return digest.hexdigest()


def one_rectangle(v, degrees, item, weights, integer, normalizer, bits, guard, progress=None):
    """Return actual-r -> actual-f -> (R4,V2), without global factors."""
    lo, hi, single, marker = item
    u, w = marker
    if min(u, w) <= 0 or any(type(r) is not int for r in weights):
        raise ValueError("positive markers and integer actual weights required")
    if single != max(weights)-2*lo or not 0 <= lo <= hi <= max(weights)//2:
        raise ValueError("rectangle does not match the requested maximum weight")
    emit = (lambda phase, **values: None) if progress is None else progress
    guard()
    up, wp = [u**k for k in range(65)], [w**k for k in range(65)]
    scaled = {}
    for (f, s), row in degrees.items():
        guard()
        if type(f) is not int or type(s) is not int or min(f, s) < 0 or f+s > 64 or len(row) != 6:
            raise ValueError("expected six local pair channels at total degree<=64")
        if any(not isinstance(x, (F, int)) or isinstance(x, bool) or x < 0 for x in row):
            raise ValueError("exact nonnegative rational local coefficients required")
        scaled[f, s] = tuple(up[f]*wp[s]*x for x in row)
    arithmetic = integer.IntegerArithmetic(hi, single, cyclic_full=False,
        cyclic_single=False, bits=bits, guard=guard)
    if arithmetic.length > PACKED_LIMIT:
        raise ValueError("adaptive packed-length guard failed")
    rr, tt = arithmetic.pack(scaled)
    del scaled, up, wp
    emit("exact_rectangle_packed", first_full=lo, last_full=hi, singleton_limit=single,
         u=str(u), v=str(w), packed_length=arithmetic.length)
    for step in range(6):
        guard()
        tt = arithmetic.multiply(rr, tt)
        rr = arithmetic.multiply(rr, rr)
        emit("exact_rectangle_doubling", first_full=lo, last_full=hi, doubling=step+1,
             regular_exponent=rr.exponent, terminal_exponent=tt.exponent)
    # Truncation decreases DC whereas positive quantization increases it.
    # No comparison to an untruncated DC is asserted for these exact axes.
    result, coverage = {}, []
    for r in weights:
        current = {}
        for f in range(lo, min(hi, r//2)+1):
            guard()
            s = r-2*f
            assert 0 <= s <= single and f <= hi
            denominator = normalizer.target_denominator(r, f, marker)
            assert denominator > 0
            current[f] = tuple(tuple(v.exact(x/denominator)
                for row in arithmetic.extract(matrix, f, s) for x in row)
                for matrix in (rr, tt))
            if tuple(map(len, current[f])) != (4, 2):
                raise ArithmeticError("regional extraction shape mismatch")
        result[r] = current
        coverage.append([r, len(current)])
    report = {"assigned_full_range": [lo, hi], "exact_full_limit": hi,
        "exact_singleton_limit": single, "u": str(u), "v": str(w),
        "coverage_weight_target_counts": coverage, "arithmetic": arithmetic.report(),
        "regular_exponent": rr.exponent, "terminal_exponent": tt.exponent,
        "regular_object": "H^64", "terminal_object": "H^63 J",
        "own_actual_denominators": True, "positive_aliases": False,
        "rounded_coefficients_are_exact": False, "alpha_applied": False,
        "source_factor_applied": False, "support_count_applied": False}
    del rr, tt, arithmetic
    gc.collect()
    guard()
    return result, report


def exact_controls(v, anchor, join, normalizer, guard):
    """Schedule, exact scale cancellation, odd joins, and source boundaries."""
    from flint import fmpq
    checks = 0
    for last in (225, 232, 254, 255, 256, 257, 507, 512, 544, 1024):
        guard()
        schedule = rectangle_schedule(last)
        covered = [f for lo, hi, _, _ in schedule for f in range(lo, hi+1)]
        assert covered == list(range(last//2+1))
        for r in range(max(225, last-32), last+1):
            targets = []
            for lo, hi, single, marker in schedule:
                assert single*(2*hi+1)+hi+1 <= PACKED_LIMIT
                assert min(marker) > 0
                for f in range(lo, min(hi, r//2)+1):
                    assert 0 <= r-2*f <= single
                    targets.append(f)
            assert targets == list(range(r//2+1))
            checks += 1
    # Cancellation is exact before positive quantization. A neighboring r
    # changes the singleton coefficient and cannot reuse the old numerator.
    for slots in (5, 9):
        for r in (3, 4, 5):
            for f in range(r//2+1):
                s = r-2*f
                coefficient = F(7+2*f+3*s, 31)
                multiplicity = comb(slots, f)*comb(slots-f, s)
                for marker in ((F(2, 7), F(3, 11)), (F(5, 13), F(7, 17))):
                    denominator = normalizer.target_denominator(r, f, marker, slots=slots)
                    assert coefficient*marker[0]**f*marker[1]**s/denominator == coefficient/multiplicity
                    checks += 1
    # Test the actual exported complete join, not an inferred even-weight
    # interface. Exact enumeration in the frozen join tests also checks
    # that this run envelope dominates all indexed tiny-cycle profiles.
    for r in (3, 4):
        guard()
        regions = [((F(1+f, 9), F(2, 11), F(3, 13), F(4+f, 17)),
                    (F(5+f, 19), F(7, 23+f))) for f in range(r//2+1)]
        counts = join.run_counts(8, r)
        expected = F(0)
        for c in range(1, r+1):
            candidates = []
            for a in range(r//2+1):
                b = r-c-a
                if 0 <= b <= r//2:
                    candidates.append(join.two_orientations(regions[a], regions[b])[0])
            expected += counts[c]*max(candidates)
        directed = [tuple(tuple(v.exact(x) for x in row) for row in region) for region in regions]
        actual = join.directed_join(directed, guard, weight=r, points=8)
        point = fmpq(expected.numerator, expected.denominator)
        assert actual["bound"].lower().fmpq() <= point <= actual["bound"].upper().fmpq()
        assert actual["report"]["ordered_profiles"] == (r//2+1)**2-int(r % 2 == 0)
        checks += 2
    for r in (225, 255, 256, 257, 511, 512, 1024):
        h, c_h, eta = anchor.source_factor(r)
        expected_c = 1-F(r, 256*(r//256+1))
        expected_eta = (1+sum((F(comb(r, j), 255) for j in range(1, h+1)), F(0)))/(255*expected_c)
        assert h == r//256 and c_h == expected_c > 0 and eta == expected_eta
        if r < 256:
            old_eta = F(256**2, 255**2*(257-r))
            assert eta/old_eta == F(255*(257-r), 256*(256-r)) > 1
        checks += 1
    return {"schedule_scale_odd_join_source_checks": checks,
        "odd_and_even_directed_join_tested": True, "no_legacy_eta_substitution": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=225)
    parser.add_argument("--last", type=int, default=232)
    parser.add_argument("--z", type=rational, default=F(267, 311))
    parser.add_argument("--alpha", type=rational, default=F(119, 1280))
    parser.add_argument("--bits", type=int, choices=(128, 192, 256), default=192)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=384)
    parser.add_argument("--seconds", type=int, choices=(300, 600), default=300)
    args = parser.parse_args()
    if not 225 <= args.first <= args.last <= 1024 or args.last-args.first >= 33:
        parser.error("require225<=first<=last<=1024 and at most33 actual sizes")
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    volume = Path(__file__).resolve().anchor
    free = shutil.disk_usage(volume).free
    if free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on the source volume")
    started = time.monotonic()
    deadline = started+args.seconds
    integer = load("cheap_lane_integer_regions.py", INTEGER_PIN, "lower_exact_integer")
    peak = integer.peak_memory_reader()
    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("lower exact-band deadline; no partial certificate")
        if peak() > MEMORY:
            raise MemoryError("lower exact-band512MiB peak memory budget exceeded")
    def progress(phase, **values):
        guard()
        print(json.dumps({"phase": phase, "elapsed_seconds": time.monotonic()-started,
            "peak_process_bytes": peak(), **values}), flush=True)
    guard()
    anchor = load("cheap_lane_cone_cauchy_anchors.py", ANCHOR_PIN, "lower_exact_anchor")
    local = load("cheap_lane_ordered_local_directed.py", LOCAL_PIN, "lower_exact_local")
    normalizer = load("cheap_lane_band_regions.py", REGIONS_PIN, "lower_exact_normalizer")
    join = load("cheap_lane_directed_run_mass.py", JOIN_PIN, "lower_exact_join")
    weights = tuple(range(args.first, args.last+1))
    schedule = rectangle_schedule(args.last)
    progress("lower_exact_band_preflight", actual_weights=weights, z=str(args.z), alpha=str(args.alpha),
        bits=args.bits, rectangles=len(schedule), packed_limit=PACKED_LIMIT,
        maximum_packed_length=max(s*(2*h+1)+h+1 for _, h, s, _ in schedule))
    built = local.build_local_pairs(args.z, args.alpha, anchor, guard, deadline,
                                   args.precision, progress=progress)
    v, core = built["verifier"], built["core"]
    degrees = built["rational_coefficients"]
    before = coefficient_hash(degrees)
    assert before == built["report"]["rational_coefficient_sha256"]
    tests = {"integer": integer.exact_tests(guard), "normalizer": normalizer.exact_tests(guard),
        "join": join.exact_tests(guard), "wrapper": exact_controls(v, anchor, join, normalizer, guard)}
    regions = {r: [None]*(r//2+1) for r in weights}
    table_reports = []
    for item in schedule:
        partial, report = one_rectangle(v, degrees, item, weights, integer, normalizer,
                                        args.bits, guard, progress)
        for r in weights:
            for f, region in partial[r].items():
                if regions[r][f] is not None:
                    raise ArithmeticError("duplicate assigned marginal")
                regions[r][f] = region
        table_reports.append(report)
        progress("exact_rectangle_complete", assigned_full_range=report["assigned_full_range"],
            extracted_targets=sum(len(partial[r]) for r in weights))
        del partial
    assert before == coefficient_hash(degrees)
    if any(x is None for r in weights for x in regions[r]):
        raise ArithmeticError("incomplete actual-weight marginal coverage")
    logtwo = v.arb(2).log()
    beta_log, z_log, alpha = v.exact(core.BETA).log(), v.exact(args.z).log(), v.exact(args.alpha)
    total = v.arb(0)
    reports = []
    exact_total_count = 0
    for r in weights:
        guard()
        h, c_h, eta = anchor.source_factor(r)
        source_log = alpha*(v.exact(eta).log()+r*beta_log-core.LIMIT*z_log)
        source = source_log.exp()
        joined = join.directed_join(regions[r], guard, progress, weight=r, points=POINTS)
        contribution = v.positive(source*joined["bound"])
        if not contribution > 0:
            raise ArithmeticError("positive complete per-weight contribution required")
        total += contribution
        count = comb(POINTS, r)
        exact_total_count += count
        worst = max(range(1, r+1), key=lambda c: joined["run_bounds"][c-1].upper())
        row = {"weight": r, "log2_union_bound_interval": str(contribution.log()/logtwo),
            "below_2_pow_minus_40": bool(contribution < v.arb(2)**(-40)),
            "source_factor": {"method": "joint bounded-missing", "h": h, "c_h": str(c_h),
                "eta": str(eta), "log2_interval": str(source_log/logtwo)},
            "marginals": len(regions[r]), "ordered_profiles": joined["report"]["ordered_profiles"],
            "exact_run_count_sum_equals_binomial": joined["report"]["all_run_counts_sum_to_binomial"],
            "indexed_support_count_sha256": hashlib.sha256(str(count).encode()).hexdigest(),
            "regional_input_intersections": joined["report"]["input_nonnegative_intersections"],
            "worst_run": {"runs": worst, "log2_interval": str((source*joined["run_bounds"][worst-1]).log()/logtwo)}}
        reports.append(row)
        progress("complete_actual_weight", **row)
        del joined
    guard()
    report = {"scope": "complete actual source-size band; conservative exact run-mass envelopes",
        "directed": True, "all_requested_weights_complete": True, "whole_code_certificate": False,
        "parameters": {"first": args.first, "last": args.last, "sizes": len(weights),
            "z": str(args.z), "alpha": str(args.alpha), "requested_precision": args.precision,
            "working_precision": built["report"]["working_precision"], "regional_bits": args.bits,
            "seconds_budget": args.seconds, "memory_budget_bytes": MEMORY,
            "packed_limit": PACKED_LIMIT, "preflight_free_bytes": free},
        "log2_band_union_bound_interval": str(total.log()/logtwo),
        "band_below_2_pow_minus_40": bool(total < v.arb(2)**(-40)),
        "band_below_2_pow_minus_60": bool(total < v.arb(2)**(-60)),
        "per_weight": reports, "tables": table_reports, "tests": tests,
        "local": {key: built["report"][key] for key in ("rational_coefficient_sha256",
            "raw_table_sha256", "powered_freeze", "tests", "field_tests", "residual_rank_receipt")},
        "coefficient_table_unchanged": before == coefficient_hash(degrees),
        "source_factor_scope": "own current joint-missing eta_r; no legacy<=256 eta substitution",
        "regional_order": "e0(Ra Rb)^127 Ra Vb, with reverse terminal retained",
        "regular_region": "H^64", "terminal_region": "H^63 J; no unused final refresh",
        "source_count_outside_alpha": True, "new_local_alpha_applied": False,
        "exact_indexed_count_sum_sha256": hashlib.sha256(str(exact_total_count).encode()).hexdigest(),
        "nonnegative_intersections": v.intersections, "elapsed_seconds": time.monotonic()-started,
        "peak_process_bytes": peak(), "generated_bytes": 0,
        "resource_guards": "cooperative checks between operations; no external kill or interruption of FLINT calls",
        "pins": {"local": LOCAL_PIN, "integer": INTEGER_PIN, "normalizer": REGIONS_PIN,
                 "join": JOIN_PIN, "anchor": ANCHOR_PIN},
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
