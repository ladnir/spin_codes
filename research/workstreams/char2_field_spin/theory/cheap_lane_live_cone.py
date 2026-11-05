"""Collapse all nonzero masks with ONE fixed positive row-cone vector.

Raw M already contains destination-mask multiplicities. For positive p
on classes1..16 define Q=M00, A=max_j M0j/p_j,
L=max_j (p M_live)_j/p_j, E=p M_live,0; terminal(T0,U)=(t0,p t_live).
Then any row v<=c p remains <=c L p during a live step. An activation
enters this cone with weight A; a return and the honest terminal cost E,U.
Hence every binary skeleton's SUM over all fine masks is bounded before
alpha. Raise these coarse scalar caps to alpha BEFORE averaging routes.
No global route Jensen, no coarse-state binomial factor, no uniform-F
refresh substitution. Local raw caps retain their historical point Jensen.

The same p is used at every packet and occupancy. Its normalization p16=1
does not affect a complete zero-initialized trajectory. Test one analytic
choice from the existing joint-live cap, NOT independent local p choices.
Fixed z=.85,alpha=.05, neutral path marks, consecutive256-group supports.
Exact regional coefficients and both actual profile orientations included.
All state masks are covered for this source class; NOT all source classes.
FLOAT only, no directed or whole-code certificate.

One48-second job,20GiB free,128MiB peak, stdout only. No generated data,
frozen edits, benchmarks, parameter scans, or automatic retries.
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
PIN = "264a216fd5beb0c2b2d18826e5fc3f7cd4638561ff2aa37d9862382e36bbd0d2"
Z, ALPHA = .85, .05
LN2 = math.log(2)


def load_degree(np, lse, deadline):
    path = Path(__file__).with_name("cheap_lane_degree_conditioned_tile.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen degree helper hash mismatch")
    spec = importlib.util.spec_from_file_location("cone_degree", path)
    degree = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(degree)
    degree.np, degree.lse, degree.DEADLINE = np, lse, deadline
    return degree, degree.load()


def exact_coarse(m, terminal, p):
    live = len(p)
    q = m[0][0]
    a = max(m[0][j+1]/p[j] for j in range(live))
    l = max(sum(p[i]*m[i+1][j+1] for i in range(live))/p[j] for j in range(live))
    e = sum(p[i]*m[i+1][0] for i in range(live))
    u = sum(p[i]*terminal[i+1] for i in range(live))
    return ((q, a), (e, l)), (terminal[0], u)


def exact_tests():
    # Classes0,1,2 have destination multiplicities1,2,1. The two actual
    # masks in class1 are independently enumerated inside every skeleton.
    classes = (0, 1, 1, 2)
    d = (1, 2, 1)
    caps = [[[Fraction(1+i+2*j+r, 11+3*i+j+2*r) for j in range(3)]
             for i in range(3)] for r in range(2)]
    matrices = [[[caps[r][i][j]*d[j] for j in range(3)] for i in range(3)] for r in range(2)]
    terminals = [[Fraction(2+i+r, 13+2*i+r) for i in range(3)] for r in range(2)]
    p = (Fraction(1, 7), Fraction(1))
    coarse = [exact_coarse(matrices[r], terminals[r], p) for r in range(2)]
    scaled = [exact_coarse(matrices[r], terminals[r], [3*x for x in p]) for r in range(2)]
    checked = 0
    for route in itertools.product(range(2), repeat=4):
        for skeleton in itertools.product((0, 1), repeat=3):
            binary = (0,)+skeleton
            fine = Fraction(0)
            for masks in itertools.product(range(4), repeat=3):
                if tuple(int(x > 0) for x in masks) != skeleton:
                    continue
                states = (0,)+tuple(classes[x] for x in masks)
                term = terminals[route[3]][states[3]]
                for t in range(3):
                    term *= caps[route[t]][states[t]][states[t+1]]
                fine += term
            bound = coarse[route[3]][1][binary[3]]
            rescaled = scaled[route[3]][1][binary[3]]
            for t in range(3):
                bound *= coarse[route[t]][0][binary[t]][binary[t+1]]
                rescaled *= scaled[route[t]][0][binary[t]][binary[t+1]]
            assert fine <= bound and bound == rescaled
            checked += 1
    return {"exact_binary_skeleton_mask_bounds": checked,
            "exact_cone_scale_invariances": checked}


def cone_caps(m, t, p, np):
    out = np.empty((len(m), 2, 2))
    end = np.empty((len(m), 2, 1))
    out[:, 0, 0] = m[:, 0, 0]
    out[:, 0, 1] = np.max(m[:, 0, 1:]/p[None, :], axis=1)
    action = np.einsum("i,rij->rj", p, m[:, 1:, 1:])
    out[:, 1, 1] = np.max(action/p[None, :], axis=1)
    out[:, 1, 0] = m[:, 1:, 0]@p
    end[:, 0, 0] = t[:, 0, 0]
    end[:, 1, 0] = t[:, 1:, 0]@p
    assert np.all(m[:, 0, 1:] <= out[:, 0, 1, None]*p[None, :]*(1+2e-14))
    assert np.all(action <= out[:, 1, 1, None]*p[None, :]*(1+2e-14))
    assert np.isfinite(out).all() and (out >= 0).all()
    return out, end


def route_polys(m, t, degrees, weights, alpha, pilot):
    np = pilot.np
    p = np.zeros((3, 33, m.shape[1], m.shape[2]))
    end = np.zeros((3, 33, t.shape[1], 1))
    index = (degrees[:, 1].astype(int), degrees[:, 0].astype(int))
    np.add.at(p, index, weights[:, None, None]*m**alpha)
    np.add.at(end, index, weights[:, None, None]*t**alpha)
    return pilot.tolog(p), pilot.tolog(end)


def score_polys(a, b, degree, pilot):
    for _ in range(7):
        b, a = degree.convolution(a, b), degree.convolution(a, a)
    regions, tails = [], []
    for f, s in ((128, 0), (127, 2)):
        normalizer = math.log(math.comb(4096, f)*math.comb(4096-f, s))
        regions.append(a[s, f]-normalizer)
        tails.append(b[s, f]-normalizer)
    paths, orientations = pilot.chain(regions, tails)
    source = 12+ALPHA*(math.log2((256/255)**2)
        +256*(1280-32*math.log2(255))-(35*(1 << 21)//100)*math.log2(Z))
    return {"log2_bound": source+paths/LN2,
            "source_bits_per_orientation": source,
            "orientation_chain_bits": [x/LN2 for x in orientations]}


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed cone diagnostic accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    degree, pilot = load_degree(np, logsumexp, started+48)
    tests = exact_tests()
    tests.update(degree.tests(pilot))
    singleton, filtered = pilot.imports()
    shapes, index, degrees, weights, records, route_checks = pilot.geometry()
    m, t, refinements = pilot.raw_tables(shapes, index, records, Z, singleton, filtered)
    p = np.array([math.comb(16, j)*min(1., 256/255*255.**(-(16-j))) for j in range(1, 17)])
    coarse, terminal = cone_caps(m, t, p, np)
    fullm = m[:, [0, 16], :][:, :, [0, 16]]
    fullt = t[:, [0, 16], :]
    full = score_polys(*route_polys(fullm, fullt, degrees, weights, ALPHA, pilot), degree, pilot)
    assert abs(full["log2_bound"]-537.7089953669492) < 1e-6
    result = score_polys(*route_polys(coarse, terminal, degrees, weights, ALPHA, pilot), degree, pilot)
    assert result["log2_bound"] >= full["log2_bound"]-1e-7
    empty = degrees.sum(axis=1) == 0
    assert np.all(coarse[empty, 0, 1] == 0) and np.all(coarse[empty, 1, 0] == 0)
    full_return = fullm[:, 1, 0] > 0
    degree.guard()
    print(json.dumps({"scope": "FLOAT all-state cone bound for consecutive256-group supports",
        "whole_code_certificate": False, "directed": False,
        "all_state_masks_covered": True, "all_source_supports_covered": False,
        "parameters": {"z": Z, "alpha": ALPHA, "path_marks": 1},
        "cone_p": p.tolist(), "same_p_for_every_packet": True,
        "full_state_only_control": full, "all_state_cone": result,
        "cone_cost_bits": result["log2_bound"]-full["log2_bound"],
        "max_local_ratio_vs_full": {
            "live": float(np.max(coarse[:, 1, 1]/fullm[:, 1, 1])),
            "activation_nonempty": float(np.max(coarse[~empty, 0, 1]/fullm[~empty, 0, 1])),
            "return_when_full_cap_positive": float(np.max(coarse[full_return, 1, 0]/fullm[full_return, 1, 0])),
            "terminal": float(np.max(terminal[:, 1, 0]/fullt[:, 1, 0]))},
        "route_shapes_with_return_opened_by_partial_states": int(np.count_nonzero(
            (~full_return) & (coarse[:, 1, 0] > 0))),
        "tests": {**tests, "exact_route_normalizers": route_checks,
            "empty_input_no_birth_or_return": "passed", "old_full_state_control": "passed",
            "local_cone_inequalities": 2*len(m)*16},
        "pair_only_refinements": refinements, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "degree_pin": PIN}, indent=2))


if __name__ == "__main__":
    main()
