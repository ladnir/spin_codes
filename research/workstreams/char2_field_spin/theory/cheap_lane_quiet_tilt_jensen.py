"""A fixed-path weighted-Jensen diagnostic; not a whole-code certificate.

Keep the frozen boundary-straddling 0/full16 path and its actual paired
profiles (128,0)/(127,2). For each regional route write its raw product as
B*V: B contains quiet, birth and return caps; V contains sustained-live
caps and the output-only terminal cap. For 0<alpha<1, positive weighted
Jensen gives

  E[(B*V)^alpha] <= E[B^alpha]^(1-alpha) E[B^alpha*V]^alpha.

Both expectations use exact positive regional coefficient extraction.
Thus the fixed-profile coefficient normalizer occurs once overall. The
quiet/birth/return factors retain alpha BEFORE route averaging; this is
not the old whole-region Jensen step. No different binary paths are
pooled here, and neither score is an actual code-failure lower bound.

One fixed FLOAT anchor, z=.85 and alpha=.05; no search. At most48 seconds,
20GiB free, estimated peak<128MiB, stdout only, no generated files.
"""

from fractions import Fraction
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
PIN = "e5bd5beda090ca63e6648a07106181d733327a4d5614a9dfb8a340a98f07c483"
Z, ALPHA = .85, .05
LN2 = math.log(2)
FIELDS = ("regional_first_log2", "regional_middle_log2", "terminal_log2")


def exact_tilt_tests():
    """Integer-power tests avoid any floating comparison of radicals."""
    counts = (2, 3, 5)
    root_b = (Fraction(1, 2), Fraction(1, 3), Fraction(3, 4))
    root_v = (Fraction(2, 3), Fraction(1, 5), Fraction(4, 7))
    checks = 0
    for n in (2, 20):
        # B_i=root_b_i**n, V_i=root_v_i**n, alpha=1/n.
        lhs = sum(w*b*v for w, b, v in zip(counts, root_b, root_v))
        base = sum(w*b for w, b in zip(counts, root_b))
        mixed = sum(w*b*v**n for w, b, v in zip(counts, root_b, root_v))
        assert lhs**n <= base**(n-1)*mixed
        mass = sum(counts)
        assert (lhs/mass)**n <= (base/mass)**(n-1)*(mixed/mass)
        # A constant V makes the tilted Jensen bound exact.
        v = Fraction(3, 7)
        assert (base*v)**n == base**(n-1)*(base*v**n)
        checks += 3
    return checks


def free_packet(frozen, np):
    """Raw cap1, retaining every oriented paired-route multiplicity."""
    totals = {}
    for (f, s, _), weight in frozen.packet_terms().items():
        totals[f, s] = totals.get((f, s), Fraction(0))+weight
    out = np.full((3, 33), float("-inf"))
    checks = 0
    for f in range(33):
        for s in range(min(2, 32-f)+1):
            expected = math.comb(32, f)*math.comb(32-f, s)
            assert totals[f, s] == expected
            out[s, f] = math.log(expected)
            checks += 1
    assert out[0, 0] == 0
    return out, checks


def source_bits():
    limit = 35*(1 << 21)//100
    return 12+ALPHA*(math.log2((256/255)**2)
        +256*(1280-32*math.log2(255))-limit*math.log2(Z))


def combine_regions(left, right, left_weight, right_weight, np):
    regional = {key: [left_weight*a+right_weight*b
                     for a, b in zip(left[key], right[key])]
                for key in FIELDS}
    first, middle, terminal = (regional[key] for key in FIELDS)
    # First and last profiles differ; the254 middle regions contain127
    # of each matching. No outgoing mask or refresh is added at terminal.
    scores = [source_bits()+first[i]+127*sum(middle)+terminal[1-i]
              for i in (0, 1)]
    total = float(np.logaddexp(*(x*LN2 for x in scores))/LN2)
    return {"regional_log2": regional, "orientation_scores": scores,
            "summed_log2_score": total}


def main():
    if len(sys.argv) != 1:
        raise ValueError("this fixed-witness diagnostic takes no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    path = Path(__file__).with_name("cheap_lane_path_marks.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen path helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_quiet_tilt_path", path)
    frozen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(frozen)
    grouped = frozen.load_grouped()
    singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(started+48)
    np = paired.np
    tests = frozen.self_test(np, paired.guard)
    tests["exact_weighted_Jensen_tests"] = exact_tilt_tests()
    free, tests["unconstrained_packet_coefficients"] = free_packet(frozen, np)
    filtered, _ = paired.imports()
    args = SimpleNamespace(r=256, z=Z, alpha=ALPHA, activity_mark=1.)
    powered, powered_stats = frozen.local_polynomials(
        grouped, singleton, paired, filtered, args)
    raw_args = SimpleNamespace(r=256, z=Z, alpha=1., activity_mark=1.)
    raw, raw_stats = frozen.local_polynomials(
        grouped, singleton, paired, filtered, raw_args)
    before = frozen.straddling(powered, args, paired, {})
    unpowered = frozen.straddling(raw, args, paired, {})
    # B contains exactlyQ,A,E. The straddling path never uses quietT.
    base_polys, mixed_polys = dict(powered), dict(powered)
    for name in ("L", "U"):
        base_polys[name] = free
        mixed_polys[name] = raw[name]
    base = frozen.straddling(base_polys, args, paired, {})
    mixed = frozen.straddling(mixed_polys, args, paired, {})
    tilted = combine_regions(base, mixed, 1-ALPHA, ALPHA, np)
    old_jensen = combine_regions(unpowered, unpowered, ALPHA, 0., np)
    comparisons = 0
    for field in FIELDS:
        for old, new, unraised in zip(before[field], tilted["regional_log2"][field],
                                      unpowered[field]):
            assert old <= new+1e-8
            assert old <= ALPHA*unraised+1e-8
            comparisons += 2
    assert before["summed_log2_counterterm"] <= tilted["summed_log2_score"]+1e-7
    assert math.isclose(before["summed_log2_counterterm"], -4598.4240, abs_tol=.002)
    assert math.isclose(old_jensen["summed_log2_score"], 10679.5574, abs_tol=.002)
    tests["regional_upper_comparisons"] = comparisons
    paired.guard()
    print(json.dumps({
        "scope": "FLOAT weighted-Jensen upper bound for one fixed scalar-majorant path",
        "whole_code_certificate": False, "code_failure_lower_bound": False,
        "directed": False, "parameters": {"z": Z, "alpha": ALPHA},
        "path_counts": {"births": 256, "quiet": 24320, "active_incoming": 8192,
                        "regular_updates": 32767, "terminal_incoming": 16},
        "base_packet_caps": ["Q", "A", "E"],
        "residual_packet_caps": ["L", "U (output-only terminal)"],
        "matching_profiles": [[128, 0], [127, 2]],
        "source_supports_per_orientation": 4096,
        "per_orientation_source_bits": source_bits(),
        "before_route_Jensen_log2": before["summed_log2_counterterm"],
        "old_route_Jensen_log2": old_jensen["summed_log2_score"],
        "weighted_Jensen": tilted,
        "weighted_Jensen_loss_bits": tilted["summed_log2_score"]-before["summed_log2_counterterm"],
        "recovered_bits_vs_old_Jensen": old_jensen["summed_log2_score"]-tilted["summed_log2_score"],
        "base_regional_log2": {key: base[key] for key in FIELDS},
        "mixed_regional_log2": {key: mixed[key] for key in FIELDS},
        "local_tables": {"powered": powered_stats, "raw": raw_stats},
        "tests": tests, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "path_helper_pin": PIN}, indent=2))


if __name__ == "__main__":
    main()
