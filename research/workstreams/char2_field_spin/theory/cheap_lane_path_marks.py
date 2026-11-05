"""FLOAT lower counterterms of a SELECTED proof majorant, NOT code probabilities.

Five fixed(z,alpha) anchors; no all-parameter rejection. Actual consecutive
256-group profiles(128full,0single)/(127full,2single) alternate over256 regions;
4096 supports per orientation. The last packet is output-only. Masks0/full16
have multiplicity1. Raw gamma/joint caps include alpha1 singleton extinction
only on pair-only packets. Raise raw fixed-mask caps to alpha before route
averaging; pooled Jensen majorants using these SAME caps dominate the result.

Half weights C16f*C(16-f,s)*Csj/2^s retain the fair singleton orientation.
Normalize regional coefficients by C4096f*C(4096-f,s); keep s<=2 exactly.
Source factor per orientation:4096*(eta*beta^256*z^-L)^alpha,
eta=(256/255)^2,beta=256^160/255^32,L=floor(.35*2^21).

ClassIII(births<=K,quiet<Q): corrected raw marks0<v<1,0<u<=1 cannot
decrease a retained path. ClassII(births>K) similarly forbids a decrease
under u>1. Positive terms therefore reject those marks at that anchor only.
The q95 one-excursion path belongs to the PREVIOUS44-window certificate.
Additional separated/straddling paths explicitly escape that class.
One48s stdout run,20GiB free,<128MiB. Exact Fraction and path checks first.
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
GROUPED_SHA256 = "4526522c47c9920482f8990bb104243eada12318774ee2d0dbb68fda3704ab36"
INF = float("-inf")


def load_grouped():
    path = Path(__file__).with_name("cheap_lane_grouped_masks.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != GROUPED_SHA256:
        raise RuntimeError("frozen grouped helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_path_grouped", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def half_terms(length=16, max_single=2):
    for f in range(length+1):
        for s in range(min(max_single, length-f)+1):
            for j in range(s+1):
                weight = Fraction(math.comb(length, f)*math.comb(length-f, s)
                                  *math.comb(s, j), 2**s)
                yield f, s, (f+j, f+s-j), weight


def packet_terms(length=16):
    terms = {}
    halves = list(half_terms(length))
    for f0, s0, o0, w0 in halves:
        for f1, s1, o1, w1 in halves:
            if s0+s1 > 2:
                continue
            key = (f0+f1, s0+s1, tuple(sorted(o0+o1)))
            terms[key] = terms.get(key, Fraction(0))+w0*w1
    return terms


def convolve(a, b, limit, np, guard):
    """Positive log convolution; axis0 singletons, axis1 full pairs."""
    width = min(limit+1, a.shape[1]+b.shape[1]-1)
    result = np.full((3, width), INF)
    for s in range(3):
        for left_s in range(s+1):
            av, bv = a[left_s], b[s-left_s]
            if not np.isfinite(av).any() or not np.isfinite(bv).any():
                continue
            for f in range(width):
                guard()
                lo, hi = max(0, f-len(bv)+1), min(f, len(av)-1)
                if lo <= hi:
                    value = np.logaddexp.reduce(av[lo:hi+1]+bv[f-hi:f-lo+1][::-1])
                    result[s, f] = np.logaddexp(result[s, f], value)
    return result


def power(a, exponent, np, guard):
    result = np.full((3, 1), INF)
    result[0, 0] = 0.0
    while exponent:
        if exponent & 1:
            result = convolve(result, a, 128, np, guard)
        exponent //= 2
        if exponent:
            a = convolve(a, a, 128, np, guard)
    return result


def exact_convolve(a, b, limit):
    out = {}
    for (sa, fa), x in a.items():
        for (sb, fb), y in b.items():
            if sa+sb <= 2 and fa+fb <= limit:
                key = sa+sb, fa+fb
                out[key] = out.get(key, Fraction(0))+x*y
    return out


def path_counts(m, ell, quiet, terminal):
    state = 0
    birth_count = quiet_count = active_count = steps = 0
    sequence = []
    for i in range(m):
        sequence.append((0, 16))
        if i == 0:
            sequence.extend([(16, 16)]*ell)
        sequence.append((16, 0))
    sequence.extend([(0, 0)]*quiet)
    assert len(sequence) == 128
    for index, (h, j) in enumerate(sequence):
        assert h == state
        active_count += int(h > 0)
        if terminal and index == 127:
            assert h == 0
            break
        steps += 1
        birth_count += int(h == 0 and j > 0)
        quiet_count += int(h == 0 and j == 0)
        state = j
    assert state == 0
    return birth_count, quiet_count, active_count, steps


def uncovered_paths():
    A, E, L, Q = (0,16), (16,0), (16,16), (0,0)
    separated = [A]+[L]*29+[E]+[Q]*94+[A,E,Q]
    first = [Q]*95+[A]+[L]*32
    middle = [L]*14+[E]+[Q]*95+[A]+[L]*17
    final = [L]*14+[E]+[Q]*95+[A]+[L]*16+[(16,None)]
    return separated, first, middle, final


def check_paths():
    separated, first, middle, final = uncovered_paths()
    def scan(path, incoming):
        assert len(path)==128
        birth=quiet=active=updates=0; live=[]
        for t,(h,j) in enumerate(path):
            assert h==incoming
            active+=int(h>0)
            if j is None:
                assert t==127
                continue
            updates+=1; birth+=int(h==0 and j>0); quiet+=int(h==j==0)
            if h or j: live.append(t)
            incoming=j
        fits=not live or live[-1]-live[0]+1<=44
        return birth,quiet,active,updates,incoming,fits
    a=scan(separated,0)
    assert a==(2,95,31,128,0,False)
    end=separated[:-1]+[(0,None)]
    assert scan(end,0)==(2,94,31,127,0,False)
    f,m,t=scan(first,0),scan(middle,16),scan(final,16)
    assert f==(1,95,32,128,16,True)
    assert m==(1,95,32,128,16,False)
    assert t==(1,95,32,127,16,False)
    assert tuple(f[i]+254*m[i]+t[i] for i in range(4))==(256,24320,8192,32767)
    return 5


def self_test(np, guard):
    geometry_checks = 0
    for length in (1, 2, 3):
        terms = packet_terms(length)
        totals = {}
        for (f, s, _), value in terms.items():
            totals[f, s] = totals.get((f, s), Fraction(0))+value
        for f in range(2*length+1):
            for s in range(min(2, 2*length-f)+1):
                assert totals[f, s] == math.comb(2*length, f)*math.comb(2*length-f, s)
                geometry_checks += 1
    # Independent direct slot enumeration includes the fair singleton flip.
    direct = {}
    for symbols in itertools.product(range(4), repeat=4):
        f = symbols.count(1); s = symbols.count(2)+symbols.count(3)
        if s > 2:
            continue
        occ = []
        for pair in (symbols[:2], symbols[2:]):
            occ.extend((sum(x in (1, 2) for x in pair), sum(x in (1, 3) for x in pair)))
        key = f, s, tuple(sorted(occ))
        direct[key] = direct.get(key, Fraction(0))+Fraction(1, 2**s)
    assert direct == packet_terms(2)
    toy = {(0, 0): Fraction(1, 2), (0, 1): Fraction(2, 3),
           (1, 0): Fraction(3, 7), (2, 1): Fraction(1, 5)}
    logged = np.full((3, 2), INF)
    for key, value in toy.items():
        logged[key] = math.log(value)
    exact = {(0, 0): Fraction(1)}
    checks = 0
    for exponent in range(1, 5):
        exact = exact_convolve(exact, toy, 128)
        result = power(logged, exponent, np, guard)
        for key, value in exact.items():
            assert math.isclose(float(result[key]), math.log(value), abs_tol=3e-14)
            checks += 1
    count_checks = 0
    for m, q in ((1, 63), (1, 79), (1, 95), (2, 84), (3, 84), (5, 84)):
        ell = 128-2*m-q
        reg, end = path_counts(m, ell, q, False), path_counts(m, ell, q, True)
        assert reg == (m, q, m+ell, 128)
        assert end == (m, q-1, m+ell, 127)
        assert 255*reg[3]+end[3] == 32767
        count_checks += 1
    marker_checks = 0
    for q in (63, 79, 95):
        q0 = 256*q-1
        threshold = 256*(q+1)
        for mark in (Fraction(1, 4), Fraction(3, 4), Fraction(99, 100)):
            assert mark**(q0-threshold+1) >= 1
            marker_checks += 1
    for m, threshold in ((2, 256), (3, 512), (5, 1024)):
        for mark in (Fraction(101, 100), Fraction(2), Fraction(4)):
            assert mark**(256*m-threshold-1) >= 1
            marker_checks += 1
    for alpha in (.125, .523, 1.0):
        values, weights = (Fraction(1, 7), Fraction(3, 5)), (Fraction(1, 4), Fraction(3, 4))
        left = sum(float(w)*float(x)**alpha for w, x in zip(weights, values))
        right = float(sum(w*x for w, x in zip(weights, values)))**alpha
        assert left <= right+1e-14
    return {"geometry_coefficients": geometry_checks, "direct_oriented_slot_assignments": 4**4,
            "exact_fraction_convolutions": checks, "terminal_count_checks": count_checks,
            "exact_marker_inequalities": marker_checks, "pooled_Jensen_checks": 3,
            "uncovered_concrete_paths":check_paths()}


def local_polynomials(grouped, singleton, paired, filtered, args):
    np = paired.np
    terms = packet_terms()
    shape_list = sorted({key[2] for key in terms})
    index = {shape: i for i, shape in enumerate(shape_list)}
    shapes = np.array(shape_list, dtype=np.int16)
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
    # Match the frozen pair-only cap exactly on the s=0 route slice.
    old_m, old_t, _, _ = grouped.raw_tables(args, singleton, paired, filtered)
    pair_shapes = list(itertools.combinations_with_replacement(range(17), 2))
    old_index = {tuple(sorted((a, a, b, b))): i for i, (a, b) in enumerate(pair_shapes)}
    values = {name: np.zeros((3, 33)) for name in ("A", "L", "E", "Q", "T", "U")}
    for (f, s, shape), weight in terms.items():
        row = index[shape]
        rm, rt = (old_m[old_index[shape]], old_t[old_index[shape]]) if s == 0 else (m[row], t[row])
        caps = (rm[0,16], rm[16,16], rm[16,0], rm[0,0], rt[0,0], rt[16,0])
        for name, cap in zip(values, caps):
            values[name][s, f] += float(weight)*float(cap)**args.alpha
    logs = {}
    for name, value in values.items():
        assert np.isfinite(value).all() and (value >= 0).all()
        out = np.full(value.shape, INF)
        positive = value > 0
        out[positive] = np.log(value[positive])
        logs[name] = out
    return logs, {"occupancy_shapes": len(shapes), "positive_route_terms": len(terms),
                  "local_array_bytes": raw.nbytes+terminal.nbytes+m.nbytes+t.nbytes}


def contribution(polys, m, ell, quiet, args, paired, cache):
    np, guard = paired.np, paired.guard
    def raised(name, exponent):
        key = name, exponent
        if key not in cache:
            cache[key] = power(polys[name], exponent, np, guard)
        return cache[key]
    common = convolve(raised("A", m), raised("E", m), 128, np, guard)
    common = convolve(common, raised("L", ell), 128, np, guard)
    regular = convolve(common, raised("Q", quiet), 128, np, guard)
    end = convolve(common, raised("Q", quiet-1), 128, np, guard)
    end = convolve(end, polys["T"], 128, np, guard)
    profiles = ((128, 0), (127, 2))
    reg, tail = [], []
    for f, s in profiles:
        normalizer = math.log(math.comb(4096, f))+math.log(math.comb(4096-f, s))
        reg.append(float(regular[s, f]-normalizer)/math.log(2))
        tail.append(float(end[s, f]-normalizer)/math.log(2))
    limit = 35*(1 << 21)//100
    eta = (256/255)**2
    logbeta = 160*8-32*math.log2(255)
    source = 12+args.alpha*(math.log2(eta)+256*logbeta-limit*math.log2(args.z))
    # Matching index is region&1. Region255 has the honest final packet.
    scores = [source+128*reg[0]+127*reg[1]+tail[1],
              source+128*reg[1]+127*reg[0]+tail[0]]
    total = float(np.logaddexp(scores[0]*math.log(2), scores[1]*math.log(2))/math.log(2))
    return {"births": 256*m, "quiet_regular_steps": 256*quiet-1,
            "actual_path_output_floor": max(0, 52*256*(m+ell)-256*256),
            "packets_per_region": {"birth": m, "full_live": ell, "return": m, "quiet": quiet},
            "regional_log2_coefficients": reg, "terminal_log2_coefficients": tail,
            "per_orientation_source_log2": source, "orientation_log2_counterterms": scores,
            "summed_log2_counterterm": total}


def straddling(polys,args,paired,cache):
    np,guard=paired.np,paired.guard
    def build(factors):
        result=np.full((3,1),INF); result[0,0]=0
        for name,n in factors:
            if (name,n) not in cache:
                cache[name,n]=power(polys[name],n,np,guard)
            result=convolve(result,cache[name,n],128,np,guard)
        return [float(result[s,f]-math.log(math.comb(4096,f)*math.comb(4096-f,s)))/math.log(2)
                for f,s in ((128,0),(127,2))]
    first=build((("Q",95),("A",1),("L",32)))
    mid=build((("L",31),("E",1),("Q",95),("A",1)))
    last=build((("L",30),("E",1),("Q",95),("A",1),("U",1)))
    limit=35*(1<<21)//100
    source=12+args.alpha*(math.log2((256/255)**2)+256*(1280-32*math.log2(255))-limit*math.log2(args.z))
    scores=[source+first[i]+127*sum(mid)+last[1-i] for i in (0,1)]
    total=float(np.logaddexp(*(x*math.log(2) for x in scores))/math.log(2))
    return {"class":"III","path":"boundary_straddling","births":256,
        "quiet_regular_steps":24320,"quiet_threshold_Q":24576,"birth_thresholds_K":[256,512,1024],
        "outside_previous_44window_class":True,"actual_path_output_floor":360448,
        "regional_first_log2":first,"regional_middle_log2":mid,"terminal_log2":last,
        "orientation_log2_counterterms":scores,"summed_log2_counterterm":total}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--self-test", action="store_true")
    mode.add_argument("--run", action="store_true")
    options = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    start = time.monotonic()
    grouped = load_grouped(); singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(start+48)
    tests = self_test(paired.np, paired.guard)
    if options.self_test:
        print(json.dumps({"tests": tests, "seconds": time.monotonic()-start, "generated_bytes": 0}))
        return
    filtered, _ = paired.imports()
    points = []
    for z, alpha in ((.85, .523), (.7, .25), (.5, .125), (.99, .05), (.9, .05)):
        args = SimpleNamespace(r=256, z=z, alpha=alpha, activity_mark=1.0)
        polys, stats = local_polynomials(grouped, singleton, paired, filtered, args)
        cache = {}; results = []
        for threshold in (16384, 20480, 24576):
            quiet = threshold//256-1
            item = contribution(polys, 1, 126-quiet, quiet, args, paired, cache)
            item.update({"class": "III", "quiet_threshold_Q": threshold,
                         "birth_thresholds_K": [256, 512, 1024],
                         "rejected_raw_markers_if_positive": "0<v<1; optional0<u<=1"})
            results.append(item)
        for m, threshold in ((2, 256), (3, 512), (5, 1024)):
            item = contribution(polys, m, 44-2*m, 84, args, paired, cache)
            item.update({"class": "II", "birth_threshold_K": threshold,
                         "rejected_raw_markers_if_positive": "u>1"})
            results.append(item)
        item=contribution(polys,2,29,95,args,paired,cache)
        item.update({"class":"III","path":"two_separated_excursions",
            "quiet_threshold_Q":24576,"birth_thresholds_K":[512,1024],
            "outside_previous_44window_class":True})
        results.extend((item,straddling(polys,args,paired,cache)))
        points.append({"parameters": vars(args), "results": results, "local_tables": stats})
    paired.guard()
    print(json.dumps({"scope": "FLOAT LOWER counterterms of one selected positive proof majorant",
        "actual_probability_bound": False, "whole_code_certificate": False,
        "all_z_alpha_choices_rejected": False,
        "actual_matching_profiles": [[128, 0], [127, 2]],
        "source_supports_per_orientation": 4096, "points": points,
        "tests": tests, "seconds": time.monotonic()-start,
        "generated_bytes": 0, "estimated_peak_bytes": 128*(1 << 20),
        "grouped_source_sha256": GROUPED_SHA256}, indent=2))


if __name__ == "__main__":
    main()
