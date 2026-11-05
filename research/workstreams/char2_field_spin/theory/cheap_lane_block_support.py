"""One FLOAT block-support diagnostic; no class bound or certificate.

At r=256, z=.85, alpha=.523, compute the frozen gamma/joint-live minimum
before adding the selected-block cap. NEVER feed the new terminal cap
into the old joint-live calculation: they bound different objects.

For an incoming mask of size h, ceil(h/4)<=b<=min(h,4) blocks are active.
Put ell=sum ell_i and f=sum_active max(0,13-ell_i). Selecting one incoming
coordinate per active block and one outgoing coordinate gives, for j>0,

 K <= [q/(q-b-1)] ((q-1)/q)^b 255^(-(16-j)) g(z)^ell z^f,
 g(z)=z+(1-z)/q, q=256.

The terminal omits feedback and uses q/(q-b). Extinction j=0 is unchanged.
Maximize over compatible b and placements before applying alpha. For each
b, use the b smallest floors. The factors g^ell and z^f count disjoint
positions. Destination multiplicity and activity stay outside alpha.

Pair-only regional coefficients use FLOAT. Repeating them in both shifted
matchings is impossible for proper nonempty supports; the score is artificial.

Use python -B with --run or --self-test. Guards: 20 GiB free, 48 seconds,
arrays below 256 MiB, hash-pinned imports. No frozen source is edited.
"""

import argparse
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
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
PAIRED_SHA256 = "4397d6b39f56786243fd68347ab01cacbd19970068480061b804951976a39969"


def load_paired(deadline):
    path = Path(__file__).with_name("cheap_lane_paired_coefficient.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PAIRED_SHA256:
        raise RuntimeError("frozen paired helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_pair_engine", path)
    paired = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(paired)
    paired.DEADLINE = deadline
    return paired


def floor_test():
    cases = 0
    for shape in ((0, 0, 0, 0), (1, 1, 0, 0), (2, 8, 12, 16),
                  (3, 7, 13, 15), (16, 16, 16, 16)):
        floors = tuple(max(0, 13-x) for x in shape)
        for h in range(17):
            for b in range((h+3)//4, min(h, 4)+1):
                explicit = min(sum(floors[i] for i in active)
                               for active in itertools.combinations(range(4), b))
                assert explicit == sum(sorted(floors)[:b])
                assert b <= h <= 4*b
                cases += 1
    return {"compatible_block_placement_checks": cases}


def local_polynomials(args, paired, filtered, middle):
    np = paired.np
    shapes = np.array(sorted(set(tuple(sorted((a, a, b, b)))
                      for a, b in itertools.product(range(17), repeat=2))), dtype=np.int16)
    index = {tuple(shape): i for i, shape in enumerate(shapes)}
    local = terminal = None
    # Finish old joint caps, using original terminals, before either new cap.
    for gamma in (0, .5, 1):
        candidate, final = filtered.local_tables(shapes, args.z, args.alpha,
                                                  gamma, args.activity_mark)
        for j in range(1, 17):
            factor = min(1.0, 256/255*255.0**(-(16-j)))
            candidate[:, :, :, j] = np.minimum(candidate[:, :, :, j],
                final[:, :, :, 0]*factor**args.alpha*math.comb(16, j))
        local = candidate if local is None else np.minimum(local, candidate)
        terminal = final if terminal is None else np.minimum(terminal, final)
    original, original_terminal = local.copy(), terminal.copy()
    g = args.z+(1-args.z)/256
    for row, shape in enumerate(shapes):
        paired.guard()
        ell = int(shape.sum())
        floors = sorted(max(0, 13-int(x)) for x in shape)
        for h in range(17):
            regular_cap = terminal_cap = 0.0
            for b in range((h+3)//4, min(h, 4)+1):
                f = sum(floors[:b])
                common = (255/256)**b*g**ell*args.z**f
                regular_cap = max(regular_cap, 256/(255-b)*common)
                terminal_cap = max(terminal_cap, 256/(256-b)*common)
            flag, activity = (0, 1.0) if h == 0 else (2, args.activity_mark)
            for j in range(1, 17):
                cap = activity*math.comb(16, j)*(regular_cap*255.0**(-(16-j)))**args.alpha
                local[row, flag, h, j] = min(local[row, flag, h, j], cap)
            cap = activity*terminal_cap**args.alpha
            terminal[row, flag, h, 0] = min(terminal[row, flag, h, 0], cap)
    assert np.array_equal(local[:, :, :, 0], original[:, :, :, 0])
    changed = (local > 0) & (local < original)
    effect = {"improved_local_entries": int(changed.sum()),
              "improved_terminal_entries": int((terminal < original_terminal).sum()),
              "maximum_local_log2_gain": float(np.log2(original[changed]/local[changed]).max())
              if changed.any() else 0.0, "extinction_entries_unchanged": True}
    degree = min(32, args.r//2)
    p = np.zeros((degree+1, 4, 17, 17))
    t = np.zeros((degree+1, 4, 17, 1))
    for a, b in itertools.product(range(17), repeat=2):
        if a+b <= degree:
            coefficient = math.comb(16, a)*math.comb(16, b)
            row = index[tuple(sorted((a, a, b, b)))]
            p[a+b] += coefficient*local[row]
            t[a+b] += coefficient*terminal[row]
    assert not np.count_nonzero(p[:, 0, 1:])
    assert not np.count_nonzero(t[:, 0, 1:])
    return ((middle.to_log(p[:, 0, 0]), middle.to_log(p[:, 2])),
            (middle.to_log(t[:, 0, 0]), middle.to_log(t[:, 2])), effect,
            int(local.nbytes+terminal.nbytes+p.nbytes+t.nbytes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--self-test", action="store_true")
    mode = parser.parse_args()
    free = shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free < 20*(1 << 30):
        raise RuntimeError("requires 20 GiB free on actual source volume")
    start = time.monotonic()
    paired = load_paired(start+48)
    tests = {**floor_test(), **paired.self_test()}
    if mode.self_test:
        print(json.dumps({"tests": tests, "seconds": time.monotonic()-start,
                          "generated_bytes": 0}))
        return
    args = SimpleNamespace(r=256, z=.85, alpha=.523, activity_mark=1.0)
    filtered, middle = paired.imports()
    a, b, effect, array_bytes = local_polynomials(args, paired, filtered, middle)
    for _ in range(7):
        b, a = (paired.positive_convolution(a, b, args.r//2),
                paired.positive_convolution(a, a, args.r//2))
    denominator = math.log(math.comb(4096, args.r//2))
    region = paired.expand_coefficient(a, args.r//2, denominator)
    tail = paired.expand_coefficient(b, args.r//2, denominator)
    result = filtered.bound(args, region, tail, middle)
    regional = {"quiet00": float(region[0, 0, 0]/paired.LN2),
                "nonquiet00": float(region[2, 0, 0]/paired.LN2),
                "nonquiet16to16": float(region[2, 16, 16]/paired.LN2),
                "terminal_quiet0": float(tail[0, 0, 0]/paired.LN2)}
    print(json.dumps({"scope": "FLOAT pair-only REGIONAL majorant diagnostic",
        "whole_code_certificate": False, "actual_support_class_bound": False,
        "repeated_marginal_whole_chain_is_artificial": True,
        "parameters": vars(args), "block_support_effect": effect,
        "regional_log2": regional, "artificial_repeated_marginal_score": result,
        "artificial_restart_floor_log2": result["source_log2_prefactor"]+256*regional["nonquiet00"],
        "tests": tests, "initial_array_bytes": array_bytes,
        "seconds": time.monotonic()-start, "generated_bytes": 0,
        "paired_source_sha256": PAIRED_SHA256}, indent=2))


if __name__ == "__main__":
    main()
