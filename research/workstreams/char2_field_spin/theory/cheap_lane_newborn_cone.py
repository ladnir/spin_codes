"""FLOAT cone retaining the occupied-block count of the previous birth.

Six states: zero, newborn1..4, mature. A birth in b occupied RS blocks
has actual support h>=17-4b. Keep that restriction for the next packet.
The next live transition enters the mature cone, a return uses the
restricted support sum, and the terminal uses restricted output only.
In particular b_previous+b_current<=4 excludes an immediate return.
The label b_previous is determined by the route, not another union factor.
Same fixed p, same actual-mask comparison, same source and route factors.
Apply alpha BEFORE route averaging, after the restricted cone bound.

Only consecutive256-group sources at z=.85,alpha=.05. All masks included,
but neither all source supports nor a directed or whole-code certificate.
One48-second job,20GiB free,128MiB estimated peak, stdout only.
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

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
PIN = "a3eff28456bad943f64932c8f10219ad316d9baede9220c3fcc380d080361abf"
FLOOR_PIN = "2143482d17afe3624bd0f445773c561a47c60cf7ec45ff8bcef57ff3417b3588"


def load_return():
    path = Path(__file__).with_name("cheap_lane_disjoint_return.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen disjoint return hash mismatch")
    spec = importlib.util.spec_from_file_location("newborn_return", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_floor():
    path = Path(__file__).with_name("cheap_lane_return_floor.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != FLOOR_PIN:
        raise RuntimeError("frozen return-floor hash mismatch")
    spec = importlib.util.spec_from_file_location("newborn_floor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def newborn_tables(m, t, p, types, thresholds, coarse, terminal, np):
    count = len(thresholds)
    mature = count+1
    out = np.zeros((len(m), count+2, count+2))
    end = np.zeros((len(m), count+2, 1))
    out[:, 0, 0] = coarse[:, 0, 0]
    out[:, mature, 0] = coarse[:, 1, 0]
    out[:, mature, mature] = coarse[:, 1, 1]
    end[:, 0, 0] = terminal[:, 0, 0]
    end[:, mature, 0] = terminal[:, 1, 0]
    for birth, hmin in enumerate(thresholds, 1):
        selected = types == birth
        out[selected, 0, birth] = coarse[selected, 0, 1]
        if hmin > 1:
            assert np.all(m[selected, 0, 1:hmin] == 0.)
        restricted = p.copy()
        restricted[:hmin-1] = 0.
        action = np.einsum("i,rij->rj", restricted, m[:, 1:, 1:])
        out[:, birth, mature] = np.max(action/p[None, :], axis=1)
        # Both independently bound the restricted actual-mask return.
        out[:, birth, 0] = np.minimum(m[:, 1:, 0]@restricted, coarse[:, 1, 0])
        end[:, birth, 0] = t[:, 1:, 0]@restricted
        assert np.all(out[:, birth, mature] <= coarse[:, 1, 1]*(1+3e-14))
        assert np.all(end[:, birth, 0] <= terminal[:, 1, 0]*(1+3e-14))
    assert np.all(out[types == 0, 0, 1:] == 0.)
    return out, end


def exact_tests(cone):
    # Two nonzero cardinalities, with two actual masks in class1. Birth
    # type1 reaches only class2; type2 can reach both. This tests exact
    # mask enumeration, newborn labels, restrictions, and honest terminal.
    classes = (0, 1, 1, 2)
    d = (1, 2, 1)
    p = (Fraction(1, 7), Fraction(1))
    caps, matrices, terminals, coarse, aged = [], [], [], [], []
    for r in range(3):
        cap = [[Fraction(2+i+2*j+r, 13+i+3*j+r) for j in range(3)] for i in range(3)]
        for j in range(1, 3):
            if j < 3-r:
                cap[0][j] = 0
        for h in range(1, 3):
            if h > r:
                cap[h][0] = 0
        m = [[cap[i][j]*d[j] for j in range(3)] for i in range(3)]
        t = [Fraction(1+i+r, 17+2*i+r) for i in range(3)]
        c, ct = cone.exact_coarse(m, t, p)
        a = [[Fraction(0) for _ in range(4)] for _ in range(4)]
        at = [Fraction(0)]*4
        a[0][0], a[3][0], a[3][3] = c[0][0], c[1][0], c[1][1]
        at[0], at[3] = ct
        if r:
            a[0][r] = c[0][1]
        for born, hmin in ((1, 2), (2, 1)):
            restricted = [p[h-1] if h >= hmin else 0 for h in (1, 2)]
            a[born][3] = max(sum(restricted[h-1]*m[h][j] for h in (1, 2))/p[j-1] for j in (1, 2))
            a[born][0] = sum(restricted[h-1]*m[h][0] for h in (1, 2))
            at[born] = sum(restricted[h-1]*t[h] for h in (1, 2))
        caps.append(cap); matrices.append(m); terminals.append(t)
        coarse.append((c, ct)); aged.append((a, at))
    checked = 0
    for route in itertools.product(range(3), repeat=4):
        for skeleton in itertools.product((0, 1), repeat=3):
            binary = (0,)+skeleton
            labels = [0]
            for k in range(3):
                if not binary[k+1]:
                    labels.append(0)
                elif not binary[k]:
                    labels.append(route[k] if route[k] else 1)
                else:
                    labels.append(3)
            fine = Fraction(0)
            for masks in itertools.product(range(4), repeat=3):
                if tuple(int(x > 0) for x in masks) != skeleton:
                    continue
                states = (0,)+tuple(classes[x] for x in masks)
                term = terminals[route[3]][states[3]]
                for k in range(3):
                    term *= caps[route[k]][states[k]][states[k+1]]
                fine += term
            value = aged[route[3]][1][labels[3]]
            old = coarse[route[3]][1][binary[3]]
            for k in range(3):
                value *= aged[route[k]][0][labels[k]][labels[k+1]]
                old *= coarse[route[k]][0][binary[k]][binary[k+1]]
            assert fine <= value <= old
            checked += 1
    return {"exact_newborn_mask_skeleton_bounds": checked}


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed newborn diagnostic accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    ret = load_return()
    cone = ret.load_cone()
    tests = exact_tests(cone)
    degree, pilot = cone.load_degree(np, logsumexp, started+48)
    singleton, filtered = pilot.imports()
    shapes, index, degrees, weights, records, route_checks = pilot.geometry()
    m, t, _ = pilot.raw_tables(shapes, index, records, cone.Z, singleton, filtered)
    p = np.asarray([math.comb(16, j)*min(1., 256/255*255.**(-(16-j))) for j in range(1, 17)])
    coarse, terminal = cone.cone_caps(m, t, p, np)
    floor = load_floor()
    return_cap, _, _, _, _ = floor.floor_caps(records, cone.Z, ret, np)
    coarse[:, 1, 0] = np.minimum(coarse[:, 1, 0], return_cap)
    old = cone.score_polys(*cone.route_polys(coarse, terminal, degrees, weights, cone.ALPHA, pilot), degree, pilot)
    assert abs(old["log2_bound"]-3262.2132476853913) < 1e-6
    types = np.asarray([sum(x > 0 for x in record[0]) for record in records])
    matrix, end = newborn_tables(m, t, p, types, (13, 9, 5, 1), coarse, terminal, np)
    excluded = 0
    for b in range(1, 5):
        excluded += int(np.count_nonzero((types+b <= 4) & (types > 0)))
        assert np.all(matrix[types+b <= 4, b, 0] == 0.)
    result = cone.score_polys(*cone.route_polys(matrix, end, degrees, weights, cone.ALPHA, pilot), degree, pilot)
    assert result["log2_bound"] <= old["log2_bound"]+1e-7
    degree.guard()
    print(json.dumps({"scope": "FLOAT newborn support cone, consecutive256-group sources",
        "whole_code_certificate": False, "directed": False,
        "all_state_masks_covered": True, "all_source_supports_covered": False,
        "parameters": {"z": cone.Z, "alpha": cone.ALPHA, "path_marks": 1},
        "old_return_floor_control": old, "newborn_cone": result,
        "improvement_bits": old["log2_bound"]-result["log2_bound"],
        "immediate_return_records_excluded": excluded,
        "tests": {**tests, "exact_route_normalizers": route_checks,
            "short_excursion_support_exclusion": "passed", "old_control": "passed"},
        "seconds": time.monotonic()-started, "estimated_peak_bytes": 128*(1 << 20),
        "generated_bytes": 0, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "return_pin": PIN, "floor_pin": FLOOR_PIN}, indent=2))


if __name__ == "__main__":
    main()
