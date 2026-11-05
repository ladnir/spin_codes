"""FLOAT disjoint-pair diagnostic: sum interior state masks before alpha.

For fixed packet occupancies, let K0,K1 be raw per-fixed-destination scalar
caps and D_j=C(16,j). Store M=K D. Sum the interior masks at alpha=1:

    B=(M0 M1) D^-1,       pair_transfer=B^alpha D.

Powers are entrywise. Destination multiplicity is restored only afterward.
The terminal pair uses M0 T1 and has no final destination multiplicity.
Each packet chooses gamma deterministically from0,1/2,1 before its output
zero expansion. Old joint-live caps use their original terminal tables.
The singleton extinction moment is averaged at alpha1 before the pair
power; shared points within a half-packet remain coupled. Conditional on
all occupancies, subsets in different packets are independent. Jensen at
this step is permitted; no already-alpha table is reused as a raw kernel.

Keep quiet and nonquiet paths as separate positive classes. A pair is quiet
exactly when its incoming and intermediate states are zero. For incoming
zero, the nonquiet sum uses only nonzero intermediate masks. For incoming
nonzero, it uses all intermediate masks. No subtraction is used.

Both occupancies remain explicit while powering. Their pair-slot counts
C(16,a)C(16,b) stay outside alpha, with a factor2 only for unequal unordered
halves. A region has64 disjoint pairs. Its last pair contains one regular
packet followed by the output-only terminal. Six positive degree doublings
preserve all128 packets and32767 total updates in the repeated-chain score.
The comparison baseline uses the same raw tables but powers each packet.

Only r256,z.85,alpha.523 is tested. Repeating the fully paired marginal in
both shifted matchings is impossible for a proper nonempty source support.
The score is artificial, NOT a support-class bound or certificate.
Use python -B --self-test or --run. Guards:20GiB free,48seconds,<256MiB.
No frozen edits or generated data. The dominated block/nonzero-moment caps
are omitted; gamma/joint and alpha1 singleton caps are retained.
"""

import argparse
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
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
SINGLETON_SHA256 = "992e6b7bd02e4785b9bfc4ae61437625afabb28a430e085fe32d60a3ab4a3453"


def load_singleton():
    path = Path(__file__).with_name("cheap_lane_singleton_extinction.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != SINGLETON_SHA256:
        raise RuntimeError("frozen singleton helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_singleton", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def raw_tables(args, singleton, paired, filtered):
    np = paired.np
    halves = list(itertools.combinations_with_replacement(range(17), 2))
    shapes = np.array([(a, a, b, b) for a, b in halves], dtype=np.int16)
    raw = terminal = None
    for gamma in (0, .5, 1):
        local, final = filtered.local_tables(shapes, args.z, 1.0, gamma, 1.0)
        for j in range(1, 17):
            factor = min(1.0, 256/255*255.0**(-(16-j)))
            local[:, :, :, j] = np.minimum(local[:, :, :, j],
                final[:, :, :, 0]*factor*math.comb(16, j))
        raw = local if raw is None else np.minimum(raw, local)
        terminal = final if terminal is None else np.minimum(terminal, final)
    m = raw[:, 0]+raw[:, 2]
    t = terminal[:, 0]+terminal[:, 2]
    singleton_half = singleton.half_tables(args.z, 1.0)
    for row, (a, b) in enumerate(halves):
        for h in range(1, 17):
            m[row, h, 0] = min(m[row, h, 0],
                singleton.packet_cap(singleton_half, a, b, h))
    degrees = np.array([a+b for a, b in halves], dtype=np.int16)
    weights = np.array([math.comb(16, a)*math.comb(16, b)*(1 if a == b else 2)
                        for a, b in halves], dtype=float)
    assert len(halves) == 153 and math.isclose(float(weights.sum()), 2.0**32)
    return m, t, degrees, weights


def pair_batch(mi, all_m, all_t, alpha, multiplicities, np):
    """One first occupancy, all second occupancies; positive activity split."""
    quiet = mi[0, 0]*all_m[:, 0, :]
    nonquiet = np.matmul(mi[None, :, :], all_m)
    nonquiet[:, 0, :] = np.matmul(mi[None, 0:1, 1:], all_m[:, 1:, :])[:, 0]
    terminal_quiet = mi[0, 0]*all_t[:, 0, :]
    terminal_nonquiet = np.matmul(mi[None, :, :], all_t)
    terminal_nonquiet[:, 0, :] = np.matmul(mi[None, 0:1, 1:], all_t[:, 1:, :])[:, 0]
    d = multiplicities
    return ((quiet/d)**alpha*d, (nonquiet/d)**alpha*d,
            terminal_quiet**alpha, terminal_nonquiet**alpha)


def assemble(args, m, t, degrees, weights, paired, middle):
    np = paired.np
    d = np.array([math.comb(16, j) for j in range(17)], dtype=float)
    q = np.zeros((65, 17)); n = np.zeros((65, 17, 17))
    tq = np.zeros((65, 1)); tn = np.zeros((65, 17, 1))
    for i in range(len(m)):
        paired.guard()
        bq, bn, btq, btn = pair_batch(m[i], m, t, args.alpha, d, np)
        route = weights[i]*weights
        degree = degrees[i]+degrees
        np.add.at(q, degree, route[:, None]*bq)
        np.add.at(n, degree, route[:, None, None]*bn)
        np.add.at(tq, degree, route[:, None]*btq)
        np.add.at(tn, degree, route[:, None, None]*btn)
    # Same raw table, but move alpha before the interior mask sum.
    bq = np.zeros((33, 17)); bn = np.zeros((33, 17, 17))
    btq = np.zeros((33, 1)); btn = np.zeros((33, 17, 1))
    per_packet = (m/d)**args.alpha*d
    per_terminal = t**args.alpha
    np.add.at(bq, degrees, weights[:, None]*per_packet[:, 0])
    np.add.at(bn[:, 1:], degrees, weights[:, None, None]*per_packet[:, 1:])
    np.add.at(btq, degrees, weights[:, None]*per_terminal[:, 0])
    np.add.at(btn[:, 1:], degrees, weights[:, None, None]*per_terminal[:, 1:])
    arrays = (q, n, tq, tn, bq, bn, btq, btn)
    assert all(np.isfinite(x).all() and (x >= 0).all() for x in arrays)
    return ((middle.to_log(q), middle.to_log(n)),
            (middle.to_log(tq), middle.to_log(tn)),
            (middle.to_log(bq), middle.to_log(bn)),
            (middle.to_log(btq), middle.to_log(btn)),
            sum(x.nbytes for x in arrays))


def toy_test(paired):
    # Categories0,1,2 represent1,2,1 actual masks. No incoming-mask count.
    d = (1, 2, 1)
    masks = (0, 1, 1, 2)
    kernels = []
    terminals = []
    packets = []
    final_packets = []
    for occupancy in range(2):
        k = [[Fraction(1+i+occupancy, 7+2*i+3*j+occupancy)
              for j in range(3)] for i in range(3)]
        t = [[Fraction(2+i+occupancy, 11+2*i+occupancy)] for i in range(3)]
        p = paired.exact_zero(3, 3); end = paired.exact_zero(3, 1)
        for i in range(3):
            flag = int(i > 0)
            p[flag][i] = [k[i][j]*d[j] for j in range(3)]
            end[flag][i] = t[i]
        kernels.append(k); terminals.append(t)
        packets.append(p); final_packets.append(end)
    expected = paired.exact_poly_product(packets, packets, 2)
    expected_t = paired.exact_poly_product(packets, final_packets, 2)
    actual = [paired.exact_zero(3, 3) for _ in range(3)]
    actual_t = [paired.exact_zero(3, 1) for _ in range(3)]
    comparisons = 0
    np = paired.np
    for a, b in itertools.product(range(2), repeat=2):
        ma = np.array([[float(kernels[a][i][j]*d[j]) for j in range(3)] for i in range(3)])
        mb = np.array([[[float(kernels[b][i][j]*d[j]) for j in range(3)] for i in range(3)]])
        tb = np.array([[[float(terminals[b][i][0])] for i in range(3)]])
        q, n, tq, tn = pair_batch(ma, mb, tb, 1.0, np.array(d), np)
        for i in range(3):
            for mid in masks:
                flag = int(i > 0 or mid > 0)
                for end in masks:
                    actual[a+b][flag][i][end] += kernels[a][i][mid]*kernels[b][mid][end]
                actual_t[a+b][flag][i][0] += kernels[a][i][mid]*terminals[b][mid][0]
        single = paired.exact_flag_product(packets[a], packets[b])
        single_t = paired.exact_flag_product(packets[a], final_packets[b])
        for i, j in itertools.product(range(3), repeat=2):
            assert math.isclose(n[0, i, j], float(single[1][i][j]), rel_tol=3e-14, abs_tol=1e-100)
            if i == 0:
                assert math.isclose(q[0, j], float(single[0][i][j]), rel_tol=3e-14)
            comparisons += 1
        for i in range(3):
            assert math.isclose(tn[0, i, 0], float(single_t[1][i][0]), rel_tol=3e-14, abs_tol=1e-100)
        assert math.isclose(tq[0, 0], float(single_t[0][0][0]), rel_tol=3e-14)
        # Subadditivity: grouping a positive interior sum cannot increase it
        # relative to raising each fixed-mask contribution separately.
        q, n, tq, tn = pair_batch(ma, mb, tb, .523, np.array(d), np)
        for i, j in itertools.product(range(3), repeat=2):
            old = sum((float(kernels[a][i][mid]*kernels[b][mid][j]))**.523*d[j]
                      for mid in masks if i > 0 or mid > 0)
            assert n[0, i, j] <= old*(1+2e-14)
    assert actual == expected and actual_t == expected_t
    return {"explicit_mask_categories": [1, 2, 1], "raw_alpha1_comparisons": comparisons,
            "exact_mask_polynomial_convolution": "passed",
            "terminal_and_positive_activity_split": "passed",
            "fractional_grouping_subadditivity": "passed"}


def evaluate(args, a, b, doublings, paired, filtered, middle):
    for _ in range(doublings):
        b, a = (paired.positive_convolution(a, b, args.r//2),
                paired.positive_convolution(a, a, args.r//2))
    denominator = math.log(math.comb(4096, args.r//2))
    region = paired.expand_coefficient(a, args.r//2, denominator)
    tail = paired.expand_coefficient(b, args.r//2, denominator)
    result = filtered.bound(args, region, tail, middle)
    regional = {"quiet00": float(region[0, 0, 0]/paired.LN2),
                "nonquiet00": float(region[2, 0, 0]/paired.LN2),
                "nonquiet16to16": float(region[2, 16, 16]/paired.LN2)}
    return {"regional_log2": regional, "artificial_chain": result,
            "artificial_restart_floor_log2": result["source_log2_prefactor"]+256*regional["nonquiet00"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--self-test", action="store_true")
    modes.add_argument("--run", action="store_true")
    options = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    start = time.monotonic()
    singleton = load_singleton(); block = singleton.load_block()
    paired = block.load_paired(start+48)
    tests = {**toy_test(paired), **singleton.self_test(), **paired.self_test()}
    paired.guard()
    if options.self_test:
        print(json.dumps({"tests": tests, "seconds": time.monotonic()-start, "generated_bytes": 0}))
        return
    args = SimpleNamespace(r=256, z=.85, alpha=.523, activity_mark=1.0)
    filtered, middle = paired.imports()
    m, t, degrees, weights = raw_tables(args, singleton, paired, filtered)
    a, b, baseline_a, baseline_b, array_bytes = assemble(args, m, t, degrees, weights, paired, middle)
    grouped = evaluate(args, a, b, 6, paired, filtered, middle)
    baseline = evaluate(args, baseline_a, baseline_b, 7, paired, filtered, middle)
    print(json.dumps({"scope": "FLOAT artificial repeated-marginal grouped-mask diagnostic",
        "whole_code_certificate": False, "actual_support_class_bound": False,
        "parameters": vars(args), "gamma_minimum": [0, .5, 1],
        "singleton_point_average_alpha": 1, "activity_mark": 1,
        "grouped": grouped, "same_raw_per_packet_baseline": baseline,
        "score_change": grouped["artificial_chain"]["nonquiet_log2_bound"]
                        -baseline["artificial_chain"]["nonquiet_log2_bound"],
        "tests": tests, "polynomial_array_bytes": array_bytes,
        "seconds": time.monotonic()-start, "generated_bytes": 0,
        "singleton_source_sha256": SINGLETON_SHA256}, indent=2))


if __name__ == "__main__":
    main()
