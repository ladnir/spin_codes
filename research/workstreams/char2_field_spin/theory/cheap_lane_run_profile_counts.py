"""Exact indexed support counts for two alternating matchings of a cycle.

For M pairs and r=2w active groups, put v=M-w. Matching full-pair counts
f0=w-a,f1=w-b give c=a+b one-runs. Supports remain indexed, not rotated.

Mark one 01 boundary and list the parities of the 2c alternating 01/10
boundaries. Exactly2a are even. If o one-runs and t zero-gaps have odd
lengths, their positive length compositions number
 C(w+o/2-1,c-1) C(v+t/2-1,c-1),
provided o,t are even. Let H_c(d,o,t) count these ordered parity words.
With D=diag(u,1) and E(x)=[[1,x],[x,1]], its positive transfer is
 [u^d x^o y^t] tr((D E(x) D E(y))^c).
The one-pair matrix is [[u^2+uxy,u^2y+ux],[ux+y,uxy+1]]. Thus
 N(f0,f1) = (M/c) sum_{i,j} H_c(2a,2i,2j)
                            C(w+i-1,c-1) C(v+j-1,c-1).
Each starting parity has M locations. Division by c removes the marked
run, NOT rotations. The final integer divisibility is checked.

Summing over profiles with c runs gives
 (2M/c) C(r-1,c-1) C(2M-r-1,c-1).
Constant supports are separate; invalid profiles count zero; odd weights
raise ValueError. Exact tests cover every even-weight profile on cycles
6,8,10,12. CLI: size256,1..4 runs (optionally6). Counts only, not distance
certificates.48-second/20GiB guards; stdout only, no generated files.
"""

import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
DEADLINE = float("inf")


def guard():
    if time.monotonic() > DEADLINE:
        raise RuntimeError("48-second deadline exhausted; no automatic retry")


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(maxsize=32)
def parity_histogram(runs):
    """Positive pair-step recurrence for the displayed 2x2 transfer trace."""
    if not isinstance(runs, int) or isinstance(runs, bool) or runs < 1:
        raise ValueError("a positive integer run count is required")
    answer = defaultdict(int)
    for first in (0, 1):
        # Count the initial 01 boundary once.
        states = {(first, int(first == 0), 0, 0): 1}
        for step in range(runs):
            following = defaultdict(int)
            for (current, even, odd_one, odd_zero), count in states.items():
                for end in (0, 1):  # The next 10 boundary.
                    d = even + int(end == 0)
                    o = odd_one + int(current != end)
                    if step + 1 == runs:
                        # Close without recounting the first boundary.
                        answer[d, o, odd_zero + int(end != first)] += count
                    else:
                        for start in (0, 1):
                            key = (start, d + int(start == 0), o,
                                   odd_zero + int(end != start))
                            following[key] += count
            states = following
            guard()
    assert sum(answer.values()) == 1 << (2*runs)
    return dict(answer)


def profile_count(pairs, weight, first_full, second_full):
    """Count one ordered matching profile, including indexed rotations."""
    values = (pairs, weight, first_full, second_full)
    if any(not isinstance(x, int) or isinstance(x, bool) for x in values):
        raise TypeError("all count parameters must be integers")
    if pairs < 1 or not 0 <= weight <= 2*pairs:
        raise ValueError("require M>=1 and 0<=r<=2M")
    if weight % 2:
        raise ValueError("this interface supports even source weights only")
    if weight == 0:
        return int(first_full == second_full == 0)
    if weight == 2*pairs:
        return int(first_full == second_full == pairs)
    w, v = weight//2, pairs-weight//2
    a, b = w-first_full, w-second_full
    if not (0 <= a <= min(w, v) and 0 <= b <= min(w, v)):
        return 0
    runs = a+b
    if not 1 <= runs <= min(weight, 2*pairs-weight):
        return 0
    total = 0
    for (even, odd_one, odd_zero), count in parity_histogram(runs).items():
        if even != 2*a or odd_one % 2 or odd_zero % 2:
            continue
        total += (count * choose(w+odd_one//2-1, runs-1)
                  * choose(v+odd_zero//2-1, runs-1))
    numerator = pairs*total
    assert numerator % runs == 0
    return numerator//runs


def run_total(pairs, weight, runs):
    if not 0 < weight < 2*pairs:
        return int(runs == 0 and weight in (0, 2*pairs))
    if not 1 <= runs <= min(weight, 2*pairs-weight):
        return 0
    numerator = (2*pairs * comb(weight-1, runs-1)
                 * comb(2*pairs-weight-1, runs-1))
    assert numerator % runs == 0
    return numerator//runs


def checks():
    masks = comparisons = totals = 0
    for pairs in (3, 4, 5, 6):
        exact = defaultdict(int)
        for mask in range(1 << (2*pairs)):
            masks += 1
            weight = mask.bit_count()
            if weight % 2:
                continue
            first = sum(((mask >> (2*i)) & 3) == 3 for i in range(pairs))
            second = sum(((mask >> (2*i+1)) & 1)
                         * ((mask >> ((2*i+2) % (2*pairs))) & 1)
                         for i in range(pairs))
            exact[weight, first, second] += 1
        for weight in range(0, 2*pairs+1, 2):
            by_runs = defaultdict(int)
            for first in range(pairs+1):
                for second in range(pairs+1):
                    got = profile_count(pairs, weight, first, second)
                    assert got == exact[weight, first, second]
                    assert got == profile_count(pairs, weight, second, first)
                    by_runs[weight-first-second] += got
                    comparisons += 1
            assert sum(by_runs.values()) == comb(2*pairs, weight)
            if 0 < weight < 2*pairs:
                for runs in range(1, min(weight, 2*pairs-weight)+1):
                    assert by_runs[runs] == run_total(pairs, weight, runs)
                    totals += 1
        guard()
    assert {key[1:]: value for key, value in parity_histogram(3).items()
            if key[0] == 2} == {(0, 2): 3, (2, 0): 3, (2, 2): 9}
    assert profile_count(4096, 256, 128, 127) == 4096
    assert profile_count(4096, 256, 128, 126) == 1031800832
    assert profile_count(4096, 256, 127, 127) == 6224347136
    assert profile_count(4096, 256, 128, 125) == 85934564093952
    assert profile_count(4096, 256, 127, 126) == 1305913168232448
    return {"exhaustive_indexed_masks": masks, "exact_profile_comparisons": comparisons,
            "exact_run_total_checks": totals, "matching_symmetry": "passed",
            "three_run_parity_histogram": "passed", "large_exact_edge_checks": 5}


def main():
    global DEADLINE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-runs", type=int, choices=range(1, 7), default=4)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on the source volume")
    started = time.monotonic(); DEADLINE = started+48
    tests = checks(); rows = []
    for runs in range(1, args.max_runs+1):
        profiles = []
        for a in range(runs+1):
            b = runs-a; first, second = 128-a, 128-b
            profiles.append({"full0": first, "full1": second,
                             "single0": 2*a, "single1": 2*b,
                             "indexed_sources": profile_count(4096, 256, first, second)})
        total = sum(row["indexed_sources"] for row in profiles)
        assert total == run_total(4096, 256, runs)
        rows.append({"runs": runs, "profiles": profiles, "total": total})
        guard()
    print(json.dumps({"scope": "exact indexed source counts, not a distance certificate",
        "whole_code_certificate": False, "arithmetic": "positive integer recurrence",
        "pairs": 4096, "weight": 256, "rows": rows, "tests": tests,
        "seconds": time.monotonic()-started, "generated_bytes": 0,
        "estimated_peak_bytes": 32*(1 << 20),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
