"""Exact indexed cycle-support counts with O(r) positive terms per profile.

On a cycle of length 2M, let r be its weight and f0,f1 the 11 counts in
the even and odd perfect matchings. Fix t, the weight at even positions.
The four vertices E0,E1,O0,O1 have directed edge multiplicities
 E->O: A=M-r+f0, B=r-t-f0, C=t-f0, D=f0;
 O->E: E=M-r+f1, F=t-f1, G=r-t-f1, H=f1.
Degrees are M-t,t,M-r+t,r-t. Delete zero-degree vertices; let tau be a
directed spanning-tree cofactor of this balanced graph.
When all four vertices remain, tau=C*E*(G+H)+D*G*(E+F).

Choose a last-exit arborescence and order the outgoing labelled edges.
At root v all orders are allowed; elsewhere the tree edge is last.
This counts tau*d_v*product_u(d_u-1)! Euler walks. Divide by parallel-edge
factorials to erase labels. Summing even roots gives sum d_v=M. Thus
 N_t=M*tau*product_u binom(d_u, first outgoing count)/product_u d_u,
Omit zero degrees. These are indexed words, not necklaces: no rotation
correction is needed, even for periodic words. Disconnected graphs give
tau=0. Sum N_t over feasible t; every division is exact, including odd r.

Four initial binomials suffice; recurrences update them along t.
Not distance certificates.48s/20GiB guards; stdout only.
"""

from collections import defaultdict
import hashlib
import importlib.util
from itertools import permutations
import json
from math import comb, prod
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
DEADLINE = float("inf")
OLD_SHA = "25e0d9221a62ffdb79d1d625e6e858ffbb8c0b671078d2bde5d2c6fb27e43ed3"


def guard():
    if time.monotonic() > DEADLINE:
        raise RuntimeError("48-second deadline exhausted; no automatic retry")


def determinant(matrix):
    n = len(matrix)
    total = 0
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        total += (-1 if inversions % 2 else 1)*prod(matrix[i][p[i]] for i in range(n))
    return total


def tree_cofactor(edges, degrees, root=None):
    a, b, c, d, e, f, g, h = edges
    adjacency = ((0, 0, a, b), (0, 0, c, d),
                 (e, f, 0, 0), (g, h, 0, 0))
    active = [i for i, degree in enumerate(degrees) if degree]
    assert all(sum(adjacency[j][i] for j in range(4)) == degrees[i]
               for i in range(4))
    root = active[0] if root is None else root
    assert root in active
    rest = [i for i in active if i != root]
    minor = [[(degrees[i] if i == j else 0)-adjacency[i][j]
              for j in rest] for i in rest]
    result = determinant(minor)
    assert result >= 0
    return result


def profile_count(pairs, weight, first_full, second_full):
    values = (pairs, weight, first_full, second_full)
    if any(not isinstance(x, int) or isinstance(x, bool) for x in values):
        raise TypeError("all count parameters must be integers")
    if pairs < 1 or not 0 <= weight <= 2*pairs:
        raise ValueError("require M>=1 and 0<=r<=2M")
    m, r, f0, f1 = values
    if not max(0, r-m) <= min(f0, f1) <= max(f0, f1) <= min(m, r//2):
        return 0
    low = max(0, r-m, f0, f1)
    high = min(m, r, r-f0, r-f1)
    a, e = m-r+f0, m-r+f1
    bins = [comb(m-low, a), comb(low, f0),
            comb(m-r+low, e), comb(r-low, f1)]
    answer = 0
    for t in range(low, high+1):
        guard()
        b, c, d = r-t-f0, t-f0, f0
        f, g, h = t-f1, r-t-f1, f1
        degrees = (m-t, t, m-r+t, r-t)
        if min(degrees) > 0:
            tau = c*e*(g+h)+d*g*(e+f)
        else:
            tau = tree_cofactor((a, b, c, d, e, f, g, h), degrees)
        numerator = m*tau*prod(bins)
        denominator = prod(x for x in degrees if x)
        term, remainder = divmod(numerator, denominator)
        assert remainder == 0
        answer += term
        if t != high:
            factors = ((m-t-a, m-t), (t+1, t+1-f0),
                       (m-r+t+1, t-f1+1), (r-t-f1, r-t))
            for i, (num, den) in enumerate(factors):
                bins[i], remainder = divmod(bins[i]*num, den)
                assert remainder == 0
    return answer


def run_total(pairs, weight, runs):
    if weight in (0, 2*pairs):
        return int(runs == 0)
    if not 1 <= runs <= min(weight, 2*pairs-weight):
        return 0
    numerator = 2*pairs*comb(weight-1, runs-1)*comb(2*pairs-weight-1, runs-1)
    assert numerator % runs == 0
    return numerator//runs


def checks():
    masks = comparisons = totals = cofactors = 0
    for pairs in (3, 4, 5, 6):
        exact = defaultdict(int)
        for mask in range(1 << (2*pairs)):
            masks += 1
            weight = mask.bit_count()
            f0 = sum(((mask >> (2*i)) & 3) == 3 for i in range(pairs))
            f1 = sum(((mask >> (2*i+1)) & 1)
                     *((mask >> ((2*i+2) % (2*pairs))) & 1) for i in range(pairs))
            exact[weight, f0, f1] += 1
        for weight in range(2*pairs+1):
            by_runs = defaultdict(int)
            for f0 in range(pairs+1):
                for f1 in range(pairs+1):
                    got = profile_count(pairs, weight, f0, f1)
                    assert got == exact[weight, f0, f1]
                    assert got == profile_count(pairs, weight, f1, f0)
                    by_runs[weight-f0-f1] += got
                    comparisons += 1
            assert sum(by_runs.values()) == comb(2*pairs, weight)
            for runs in range(1, min(weight, 2*pairs-weight)+1):
                assert by_runs[runs] == run_total(pairs, weight, runs)
                totals += 1
        for (weight, f0, f1), count in exact.items():
            if not count:
                continue
            for t in range(max(1, weight-pairs+1, f0, f1),
                           min(pairs-1, weight-1, weight-f0, weight-f1)+1):
                a, b, c, d = pairs-weight+f0, weight-t-f0, t-f0, f0
                e, f, g, h = pairs-weight+f1, t-f1, weight-t-f1, f1
                degrees = (pairs-t, t, pairs-weight+t, weight-t)
                edges = (a,b,c,d,e,f,g,h)
                tau = c*e*(g+h)+d*g*(e+f)
                assert all(tree_cofactor(edges, degrees, root) == tau for root in range(4))
                cofactors += 1
        guard()
    path = Path(__file__).with_name("cheap_lane_run_profile_counts.py")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == OLD_SHA
    spec = importlib.util.spec_from_file_location("frozen_run_counts", path)
    old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    old.DEADLINE = DEADLINE
    cross = 0
    for runs in range(1, 7):
        subtotal = 0
        for a in range(runs+1):
            f0, f1 = 128-a, 128-(runs-a)
            got = profile_count(4096, 256, f0, f1)
            assert got == old.profile_count(4096, 256, f0, f1)
            subtotal += got; cross += 1
        assert subtotal == run_total(4096, 256, runs)
    old.parity_histogram.cache_clear()
    return {"exhaustive_masks": masks, "all_weight_profile_checks": comparisons,
            "run_totals": totals, "tree_cofactors": cofactors,
            "frozen_counter_crosschecks": cross}


def main():
    global DEADLINE
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires 20GiB free on the source volume")
    started = time.monotonic(); DEADLINE = started+48
    tests = checks(); rows = []
    for runs in (16, 32, 64, 128, 256):
        a, b = runs//2, runs-runs//2
        f0, f1 = 128-a, 128-b
        rows.append({"runs": runs, "full0": f0, "full1": f1,
                     "indexed_sources": profile_count(4096, 256, f0, f1),
                     "all_profiles_at_runs": run_total(4096, 256, runs)})
    print(json.dumps({"scope": "exact counts, not a distance certificate",
        "whole_code_certificate": False, "arithmetic": "positive exact integers",
        "tests": tests, "sampled_profiles": rows, "generated_bytes": 0,
        "estimated_peak_bytes": 8*(1 << 20), "seconds": time.monotonic()-started,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "frozen_counter_sha256": OLD_SHA}, indent=2))


if __name__ == "__main__":
    main()
