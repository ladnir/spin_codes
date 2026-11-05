"""Isolate within-region route Jensen on one budget-admissible mask path.

Use the frozen boundary-straddling path, actual consecutive256-group
profiles, z=.85 and alpha=.05. Its256 births,24320 quiet transitions and
8192 active incoming packets satisfy the intended complementary branch
k<=512, q<=24575 and4a<=61503. It is a mask path, NOT a proven codeword.

BEFORE averages products of raw caps**alpha over exact regional routes.
AFTER averages products of raw caps over the SAME routes, then takes
alpha once per region. Source factors and two orientations stay outside
that comparison. Both use the same gamma/joint/shared-point caps. There
is no intermediate-mask summation. AFTER is a positive contribution to
the whole-region pooled majorant; BEFORE is its Jensen lower comparison.
Neither positive value is evidence of actual code failure.

For allowed raw marks u,v,w<=1, the corrected path gains factors
u**(-256alpha), v**(-255alpha), w**(-28735alpha), all>=1. Thus a positive
AFTER score rejects the selected pooled majorant at these z,alpha for
every such mark, even though the BEFORE score may be negative.

One48-second invocation,20GiB free,<128MiB estimated peak, stdout only.
No parameter search, output files, caches or frozen source edits.
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
Z = .85
ALPHA = .05
LN2 = math.log(2)


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed-path A/B accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    path = Path(__file__).with_name("cheap_lane_path_marks.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen fixed-path helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_path_jensen", path)
    frozen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(frozen)
    grouped = frozen.load_grouped()
    singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(started+48)
    np = paired.np
    tests = frozen.self_test(np, paired.guard)
    filtered, _ = paired.imports()
    args = SimpleNamespace(r=256, z=Z, alpha=ALPHA, activity_mark=1.)
    powered, _ = frozen.local_polynomials(grouped, singleton, paired, filtered, args)
    before = frozen.straddling(powered, args, paired, {})
    raw_args = SimpleNamespace(r=256, z=Z, alpha=1., activity_mark=1.)
    raw, _ = frozen.local_polynomials(grouped, singleton, paired, filtered, raw_args)
    # Pass the actual alpha for the source term; raw regional coefficients
    # themselves do not depend on the alpha argument of straddling().
    raw_result = frozen.straddling(raw, args, paired, {})
    limit = 35*(1 << 21)//100
    source = 12+ALPHA*(math.log2((256/255)**2)
        +256*(1280-32*math.log2(255))-limit*math.log2(Z))
    orientations = []
    for i in (0, 1):
        regional_sum = (raw_result["regional_first_log2"][i]
            +127*sum(raw_result["regional_middle_log2"])
            +raw_result["terminal_log2"][1-i])
        orientations.append(source+ALPHA*regional_sum)
    after = float(np.logaddexp(*(v*LN2 for v in orientations))/LN2)
    assert before["summed_log2_counterterm"] <= after+1e-8
    comparisons = 0
    for field in ("regional_first_log2", "regional_middle_log2", "terminal_log2"):
        for a, b in zip(before[field], raw_result[field]):
            assert a <= ALPHA*b+1e-10
            comparisons += 1
    for u, v, w in ((Fraction(1),)*3,
                    (Fraction(1,2), Fraction(1,3), Fraction(1,4))):
        assert u**(-256)*v**(-255)*w**(-28735) >= 1
    paired.guard()
    print(json.dumps({"scope": "FLOAT before/after route Jensen for a fixed budget-admissible mask path",
        "whole_code_certificate": False, "code_failure_lower_bound": False,
        "directed": False, "parameters": {"z": Z, "alpha": ALPHA},
        "path_counts": {"births": 256, "quiet": 24320, "active_incoming": 8192,
            "regular_updates": 32767, "terminal_incoming": 16},
        "branch_limits": {"births_max": 512, "quiet_max": 24575, "active_blocks_sum_max": 61503},
        "before_route_Jensen_log2": before["summed_log2_counterterm"],
        "after_route_Jensen_log2": after,
        "route_Jensen_loss_bits": after-before["summed_log2_counterterm"],
        "per_orientation_source_bits": source,
        "before_orientation_scores": before["orientation_log2_counterterms"],
        "after_orientation_scores": orientations,
        "before_regional_log2": {key: before[key] for key in
            ("regional_first_log2", "regional_middle_log2", "terminal_log2")},
        "raw_regional_log2": {key: raw_result[key] for key in
            ("regional_first_log2", "regional_middle_log2", "terminal_log2")},
        "tests": {**tests, "regional_Jensen_comparisons": comparisons,
            "exact_corrected_mark_checks": 2},
        "seconds": time.monotonic()-started, "estimated_peak_bytes": 128*(1 << 20),
        "generated_bytes": 0, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "path_helper_pin": PIN}, indent=2))


if __name__ == "__main__":
    main()
