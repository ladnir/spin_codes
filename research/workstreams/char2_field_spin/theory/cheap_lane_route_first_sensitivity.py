"""FLOAT count derivatives for the route-first zero/full16 control.

Use the pinned route-first local polynomials, exact regional coefficients,
and finite chain. Seed count derivatives AFTER taking local fractional
powers; these are derivatives with respect to powered marks. Therefore
reported means require NO division by alpha. A raw-mark score derivative
still equals alpha*(mean-threshold)/ln2. No second regional alpha occurs.

At neutral marks, report birth, quiet and active counts, the output floor
4a, and the separate all-quiet path. No search, no parameter changes.
The expression covers only actual packet masks0/full16 and consecutive
256-group supports; it is not a full-source or whole-code certificate.

One48-second invocation,20GiB free,128MiB estimated peak, stdout only.
No generated files, caches, frozen edits, benchmarks, or automatic retries.
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
PINS = {
    "cheap_lane_route_first_control.py": "fbf6d6b6649e0962feb885db0cdaff5348f7300b11e8ec03116016da7c8053b8",
    "cheap_lane_branch_sensitivity.py": "84b6b757a3fc4e85f28c0c3e9f59f6063ae8daf70cc3704aca6feb9ebe27600b",
}


def load(name):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != PINS[name]:
        raise RuntimeError("frozen dependency hash mismatch: " + name)
    spec = importlib.util.spec_from_file_location("route_sensitivity_"+path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed diagnostic accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    route = load("cheap_lane_route_first_control.py")
    dual = load("cheap_lane_branch_sensitivity.py")
    dual.np, dual.lse, dual.DEADLINE = np, logsumexp, started+48
    pool, degree, pilot = dual.load()
    tests = dual.tests(pilot)
    frozen = route.load("cheap_lane_path_marks.py")
    grouped = frozen.load_grouped()
    singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(started+48)
    filtered, _ = paired.imports()
    args = SimpleNamespace(r=256, z=route.Z, alpha=route.ALPHA, activity_mark=1.)
    polys, _ = frozen.local_polynomials(grouped, singleton, paired, filtered, args)
    logs = np.empty((3, 33, 2, 2))
    for name, i, j in (("Q", 0, 0), ("A", 0, 1), ("E", 1, 0), ("L", 1, 1)):
        logs[:, :, i, j] = polys[name]
    terminal = np.stack((polys["T"], polys["U"]), axis=2)[:, :, :, None]
    a, b = dual.seed(logs), dual.seed(terminal, True)
    for _ in range(7):
        b, a = dual.convolution(a, b), dual.convolution(a, a)
    regions, tails = [], []
    for f, s in ((128, 0), (127, 2)):
        norm = math.log(math.comb(4096, f)*math.comb(4096-f, s))
        regions.append((a[0][s, f]-norm, a[1][s, f]))
        tails.append((b[0][s, f]-norm, b[1][s, f]))
    paths, means, _ = dual.chain(regions, tails)
    source = (12+route.ALPHA*(math.log2((256/255)**2)
        +256*(1280-32*math.log2(255))-(35*(1 << 21)//100)*math.log2(route.Z)))
    score = source+paths/route.LN2
    assert abs(score-537.7089953669492) < 1e-6
    identity = float(means[0]+means[1]+means[2]-means[3])
    assert abs(identity-32767) < 1e-6
    q = frozen.power(polys["Q"], 128, np, paired.guard)
    qt = frozen.convolve(frozen.power(polys["Q"], 127, np, paired.guard),
                         polys["T"], 128, np, paired.guard)
    quiet_regions, quiet_tails = [], []
    for f, s in ((128, 0), (127, 2)):
        norm = math.log(math.comb(4096, f)*math.comb(4096-f, s))
        quiet_regions.append(np.array([[q[s, f]-norm]]))
        quiet_tails.append(np.array([[qt[s, f]-norm]]))
    quiet, _ = pilot.chain(quiet_regions, quiet_tails)
    pure_quiet = source+quiet/route.LN2
    assert pure_quiet <= score+1e-7
    slopes = {"birth_le_512": route.ALPHA*(means[0]-512)/route.LN2,
        "quiet_le_24575": route.ALPHA*(means[1]-24575)/route.LN2,
        "floor_4a_le_61503": route.ALPHA*(4*means[2]-61503)/route.LN2}
    dual.guard()
    peak = pool.peak_working_set()
    assert peak < 128*(1 << 20)
    print(json.dumps({"scope": "FLOAT route-first0/full16 sensitivity; neutral marks",
        "whole_code_certificate": False, "full_source_class_bound": False, "directed": False,
        "log2_control": score, "log2_pure_quiet_path": pure_quiet,
        "tilted_means": {"births": float(means[0]), "quiet": float(means[1]),
            "active_incoming": float(means[2]), "terminal_active": float(means[3])},
        "slopes_bits_per_ln_raw_mark": {k: float(v) for k, v in slopes.items()},
        "count_identity": identity, "tests": tests,
        "seconds": time.monotonic()-started, "peak_working_set_bytes": peak,
        "generated_bytes": 0, "pins": PINS,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
