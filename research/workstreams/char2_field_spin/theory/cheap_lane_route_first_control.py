"""Retain regional routes while summing zero/full-state mask paths.

FLOAT control at z=.85, alpha=.05 and consecutive256-group profiles.
Raise each fixed local occupancy cap to alpha FIRST, as in the pinned
fixed-path helper. Form the two-state matrix polynomial with entries
Q,A,E,L, and terminal T,U. Extract the exact regional coefficients of
P**128 and P**127*T; divide by the exact route count OUTSIDE alpha.
The finite256-region chain then sums every zero/full16 packet path.

This deliberately pays separate path concavity at each packet. There
are no partial actual masks, no endpoint multiplicities, and no Jensen
over route occupancies. It is a diagnostic positive submajorant, NOT
an upper bound for all source words or a whole-code certificate.

Two fixed runs: unmarked; quiet mark at the previously found pooled
control's stationary point, with correction for q<=24575. No tuning.
Keep the source comparison and support multiplicity unchanged. Terminal
is output-only and gets no quiet, birth, or feedback factor.

One48-second invocation,20GiB free,128MiB estimated peak, stdout only.
No output files, caches, frozen edits, benchmarks, or automatic retries.
"""

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
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
PINS = {
    "cheap_lane_path_marks.py": "e5bd5beda090ca63e6648a07106181d733327a4d5614a9dfb8a340a98f07c483",
    "cheap_lane_degree_conditioned_tile.py": "264a216fd5beb0c2b2d18826e5fc3f7cd4638561ff2aa37d9862382e36bbd0d2",
}
Z, ALPHA = .85, .05
LN2 = math.log(2)


def load(name):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != PINS[name]:
        raise RuntimeError("frozen dependency hash mismatch: " + name)
    spec = importlib.util.spec_from_file_location("route_first_" + path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def direct_toy(np, degree):
    """Enumerate every route and mask path of three one-slot packets."""
    alpha = .3
    caps = [np.array([[.6, .4], [.2, .5]]), np.array([[.1, .7], [.3, .2]])]
    ends = [np.array([[.8], [.6]]), np.array([[.2], [.7]])]
    p = np.full((1, 2, 2, 2), -math.inf)
    t = np.full((1, 2, 2, 1), -math.inf)
    for r in (0, 1):
        p[0, r] = alpha*np.log(caps[r])
        t[0, r] = alpha*np.log(ends[r])
    actual = degree.convolution(degree.convolution(p, p, 3, 0), t, 3, 0)
    count = 0
    for f in range(4):
        exact = np.zeros((2, 1))
        for route in itertools.product((0, 1), repeat=3):
            if sum(route) != f:
                continue
            for i, j, k in itertools.product((0, 1), repeat=3):
                exact[i, 0] += float(caps[route[0]][i, j]*caps[route[1]][j, k]
                                     *ends[route[2]][k, 0])**alpha
        normalizer = math.comb(3, f)
        assert np.allclose(np.exp(actual[0, f])/normalizer, exact/normalizer,
                           rtol=3e-14, atol=3e-14)
        count += 2
    # Fractionalizing the already normalized route sum is a different bound.
    left = .5*(Fraction(1, 7)**alpha + Fraction(3, 5)**alpha)
    right = float((Fraction(1, 7) + Fraction(3, 5))/2)**alpha
    assert left < right
    return count


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed control accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    frozen = load("cheap_lane_path_marks.py")
    degree = load("cheap_lane_degree_conditioned_tile.py")
    degree.np, degree.lse, degree.DEADLINE = np, logsumexp, started+48
    pilot = degree.load()
    tests = degree.tests(pilot)
    tests["direct_three_packet_route_mask_entries"] = direct_toy(np, degree)
    grouped = frozen.load_grouped()
    singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(started+48)
    tests.update(frozen.self_test(np, paired.guard))
    filtered, _ = paired.imports()
    args = SimpleNamespace(r=256, z=Z, alpha=ALPHA, activity_mark=1.)
    polys, _ = frozen.local_polynomials(grouped, singleton, paired, filtered, args)
    fixed = frozen.straddling(polys, args, paired, {})
    assert abs(fixed["summed_log2_counterterm"]+4598.424024587422) < 1e-7
    initial = np.empty((3, 33, 2, 2))
    for name, i, j in (("Q", 0, 0), ("A", 0, 1), ("E", 1, 0), ("L", 1, 1)):
        initial[:, :, i, j] = polys[name]
    terminal = np.stack((polys["T"], polys["U"]), axis=2)[:, :, :, None]
    source = (12+ALPHA*(math.log2((256/255)**2)
        +256*(1280-32*math.log2(255))-(35*(1 << 21)//100)*math.log2(Z)))
    results = []
    for mark in (1., .002331281099898768):
        a, b = initial.copy(), terminal.copy()
        # This mark is raw: take its alpha power at the same point as the cap.
        a[:, :, 0, 0] += ALPHA*math.log(mark)
        for _ in range(7):
            b, a = degree.convolution(a, b), degree.convolution(a, a)
        regions, tails = [], []
        for f, s in ((128, 0), (127, 2)):
            logv = math.log(math.comb(4096, f)*math.comb(4096-f, s))
            regions.append(a[s, f]-logv)
            tails.append(b[s, f]-logv)
        paths, orientations = pilot.chain(regions, tails)
        correction = -ALPHA*24575*math.log2(mark)
        value = source+correction+paths/LN2
        # The same marked fixed path is included in this all-path sum.
        fixed_marked = fixed["summed_log2_counterterm"]-ALPHA*255*math.log2(mark)
        assert value >= fixed_marked-1e-7
        results.append({"quiet_raw_mark": mark, "log2_value": value,
            "source_bits_per_orientation": source, "quiet_correction_bits": correction,
            "orientation_chain_bits": [x/LN2 for x in orientations],
            "same_marked_fixed_path_log2": fixed_marked,
            "sum_over_paths_cost_bits": value-fixed_marked,
            "region_log2": [(r/LN2).tolist() for r in regions],
            "terminal_log2": [(t/LN2).tolist() for t in tails]})
    degree.guard()
    print(json.dumps({"scope": "FLOAT route-first exact regional0/full16 positive submajorant",
        "whole_code_certificate": False, "full_source_class_bound": False,
        "directed": False, "parameters": {"z": Z, "alpha": ALPHA},
        "results": results, "tests": tests, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pins": PINS}, indent=2))


if __name__ == "__main__":
    main()
