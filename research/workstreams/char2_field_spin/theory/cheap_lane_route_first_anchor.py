"""One preselected route-first anchor comparison, not a parameter search.

Compare z=.85 versus the earlier fixed-path anchor z=.9, alpha=.05,
neutral count marks. All raw caps, exact regional route coefficients,
0/full16 mask paths, source factors and terminal match the pinned control.
No other z/alpha/mark is evaluated. Positive or negative values concern
only this FLOAT restricted operator; no full-source certificate follows.

One48-second invocation,20GiB free,128MiB estimated peak, stdout only.
No output files, caches, frozen edits, benchmarks, or automatic retries.
"""

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
PIN = "fbf6d6b6649e0962feb885db0cdaff5348f7300b11e8ec03116016da7c8053b8"


def main():
    if len(sys.argv) != 1:
        raise ValueError("preselected A/B accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    path = Path(__file__).with_name("cheap_lane_route_first_control.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen route-first control hash mismatch")
    spec = importlib.util.spec_from_file_location("anchor_route_first", path)
    route = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(route)
    frozen = route.load("cheap_lane_path_marks.py")
    degree = route.load("cheap_lane_degree_conditioned_tile.py")
    degree.np, degree.lse, degree.DEADLINE = np, logsumexp, started+48
    pilot = degree.load()
    tests = degree.tests(pilot)
    tests["direct_three_packet_route_mask_entries"] = route.direct_toy(np, degree)
    grouped = frozen.load_grouped()
    singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(started+48)
    filtered, _ = paired.imports()
    results = []
    for z in (.85, .9):
        args = SimpleNamespace(r=256, z=z, alpha=route.ALPHA, activity_mark=1.)
        polys, _ = frozen.local_polynomials(grouped, singleton, paired, filtered, args)
        fixed = frozen.straddling(polys, args, paired, {})
        a = np.empty((3, 33, 2, 2))
        for name, i, j in (("Q", 0, 0), ("A", 0, 1), ("E", 1, 0), ("L", 1, 1)):
            a[:, :, i, j] = polys[name]
        b = np.stack((polys["T"], polys["U"]), axis=2)[:, :, :, None]
        for _ in range(7):
            b, a = degree.convolution(a, b), degree.convolution(a, a)
        regions, tails = [], []
        for f, s in ((128, 0), (127, 2)):
            norm = math.log(math.comb(4096, f)*math.comb(4096-f, s))
            regions.append(a[s, f]-norm)
            tails.append(b[s, f]-norm)
        paths, orientations = pilot.chain(regions, tails)
        source = (12+route.ALPHA*(math.log2((256/255)**2)
            +256*(1280-32*math.log2(255))-(35*(1 << 21)//100)*math.log2(z)))
        score = source+paths/route.LN2
        assert score >= fixed["summed_log2_counterterm"]-1e-7
        if z == .85:
            assert abs(score-537.7089953669492) < 1e-6
        results.append({"z": z, "alpha": route.ALPHA, "log2_all_binary_paths": score,
            "log2_fixed_path": fixed["summed_log2_counterterm"],
            "source_bits_per_orientation": source,
            "orientation_chain_bits": [x/route.LN2 for x in orientations]})
    degree.guard()
    print(json.dumps({"scope": "FLOAT two preselected route-first0/full16 controls",
        "whole_code_certificate": False, "full_source_class_bound": False, "directed": False,
        "results": results, "tests": tests, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "route_first_pin": PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
