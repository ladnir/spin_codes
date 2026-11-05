"""FLOAT binary-pair pooling at fixed actual r1024 matching profiles.

Start with the frozen raw two-state cone and the genuine one/two-input
return refinement. A packet has32 route slots and raw polynomials M,T.
For two packets form H=M*M and J=M*T by positive log-domain convolution.
For total route degree(f,s), let D=C(64,f)C(64-f,s). Conditional route
Jensen bounds the grouped fractional entry by D^(1-alpha)*H[f,s]^alpha;
the terminal coefficient uses J and the same64 route slots.

The old pair coefficients instead multiply the two packet polynomials
whose local caps were already raised to alpha. Both bound the same
conditional-degree object. Take their entrywise minimum BEFORE regional
composition. No destination-mask factor is added: the binary intermediate
state is summed in H, whereas the old candidate sums it after alpha.

There are64 pairs per region. The last region uses63 regular pairs and
one terminal pair M*T, preserving32767 updates and one output-only packet.
Keep the own-f rational markers, exact source counts and source factor
unchanged. Replay the old pair candidate against the frozen single-packet
engine, then the pure grouped candidate and the per-degree minimum.

All outputs are selected-profile FLOAT diagnostics, not certificates.
No optimizer, parameter scan, generated data, or frozen source edits.
Default modulus512; optional256.180s/256MiB/20GiB, stdout only.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement, product
import json
import math
from math import comb
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
RETURN_PIN = "62eff6f2f831647696f451fdc3e736e63de0484298781725b2ac3fe89f132a3d"
WEIGHT, SECONDS, MEMORY = 1024, 180, 256*(1 << 20)
PROFILES = ((384,500),(384,384))
LN2 = math.log(2)


def load_return():
    path = Path(__file__).with_name("cheap_lane_exact_small_return.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != RETURN_PIN:
        raise RuntimeError("frozen exact-small-return hash mismatch")
    spec = importlib.util.spec_from_file_location("binary_pair_small_return", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def as_logs(coefficients, np):
    if np.any(coefficients < 0) or not np.isfinite(coefficients).all():
        raise ArithmeticError("invalid positive local coefficient")
    logs = np.full(coefficients.shape, -math.inf)
    positive = coefficients > 0
    logs[positive] = np.log(coefficients[positive])
    return logs


def pair_convolution(logs, np, guard):
    """Six positive channels of M^2 and M*T; no FFT or raw underflow."""
    n = logs.shape[0]-1
    out = np.full((2*n+1,2*n+1,6), -math.inf)
    for f in range(n+1):
        guard()
        for s in range(n+1-f):
            a = logs[f,s]
            terms = np.empty_like(logs)
            terms[:,:,0] = np.logaddexp(a[0]+logs[:,:,0], a[1]+logs[:,:,2])
            terms[:,:,1] = np.logaddexp(a[0]+logs[:,:,1], a[1]+logs[:,:,3])
            terms[:,:,2] = np.logaddexp(a[2]+logs[:,:,0], a[3]+logs[:,:,2])
            terms[:,:,3] = np.logaddexp(a[2]+logs[:,:,1], a[3]+logs[:,:,3])
            terms[:,:,4] = np.logaddexp(a[0]+logs[:,:,4], a[1]+logs[:,:,5])
            terms[:,:,5] = np.logaddexp(a[2]+logs[:,:,4], a[3]+logs[:,:,5])
            target = out[f:f+n+1,s:s+n+1]
            np.logaddexp(target, terms, out=target)
    return out


def degree_power(raw_pair, alpha, slots, np):
    out = np.full(raw_pair.shape, -math.inf)
    for f in range(slots+1):
        for s in range(slots+1-f):
            mass = comb(slots,f)*comb(slots-f,s)
            out[f,s] = alpha*raw_pair[f,s]+(1-alpha)*math.log(mass)
    return out


def route_mass_checks(guard):
    """Exact positive convolution of the two32-slot multinomials."""
    half = [(f,s,comb(32,f)*comb(32-f,s)) for f in range(33) for s in range(33-f)]
    totals = [[0]*65 for _ in range(65)]
    for i,(f,s,a) in enumerate(half):
        if i%16 == 0: guard()
        for j in range(i,len(half)):
            ff,ss,b = half[j]
            totals[f+ff][s+ss] += a*b*(1 if i == j else 2)
    checks = 0
    for f in range(65):
        for s in range(65-f):
            assert totals[f][s] == comb(64,f)*comb(64-f,s)
            checks += 1
    return checks


def exact_tests(np, guard):
    """Small rational route/mask enumeration, alpha1 and ordered terminal."""
    categories = ((0,0),(1,0),(0,1))
    packet = (((F(1,3),F(2,5)),(F(3,7),F(5,11))),
              ((F(7,13),F(11,17)),(F(13,19),F(17,23))),
              ((F(19,29),F(23,31)),(F(29,37),F(31,41))))
    terminal = ((F(37,43),F(41,47)),(F(43,53),F(47,59)),(F(53,61),F(59,67)))
    coeff = np.zeros((2,2,6))
    for (f,s),m,t in zip(categories,packet,terminal):
        coeff[f,s] = [float(x) for row in m for x in row]+list(map(float,t))
    paired = pair_convolution(as_logs(coeff,np),np,guard)
    exact_pair = defaultdict(lambda: [F(0)]*6)
    for ca,cb in product(range(3),repeat=2):
        f = categories[ca][0]+categories[cb][0]
        s = categories[ca][1]+categories[cb][1]
        for h in (0,1):
            for j in (0,1):
                exact_pair[f,s][2*h+j] += sum(packet[ca][h][k]*packet[cb][k][j] for k in (0,1))
            exact_pair[f,s][4+h] += sum(packet[ca][h][k]*terminal[cb][k] for k in (0,1))
    checks = 0
    for (f,s),row in exact_pair.items():
        for k,value in enumerate(row):
            assert abs(math.exp(paired[f,s,k])-float(value)) < 2e-12
            checks += 1
    powered_one = degree_power(paired,1.0,2,np)
    assert np.array_equal(powered_one,paired); checks += 1
    # Two pairs with a terminal second pair equal the noncommuting P^3*T.
    direct = defaultdict(lambda: [F(0)]*2)
    for routes in product(range(3),repeat=4):
        f = sum(categories[x][0] for x in routes)
        s = sum(categories[x][1] for x in routes)
        for h in (0,1):
            for hidden in product((0,1),repeat=3):
                states = (h,)+hidden
                mass = terminal[routes[3]][states[3]]
                for i in range(3): mass *= packet[routes[i]][states[i]][states[i+1]]
                direct[f,s][h] += mass
    by_pairs = defaultdict(lambda: [F(0)]*2)
    for (f,s),left in exact_pair.items():
        for (ff,ss),right in exact_pair.items():
            for h in (0,1):
                by_pairs[f+ff,s+ss][h] += sum(left[2*h+k]*right[4+k] for k in (0,1))
    assert direct == by_pairs; checks += sum(len(x) for x in direct.values())
    # Alpha1 grouped and old powers agree for every coefficient and endpoint.
    for alpha in (F(1),F(1,2)):
        oldpair = pair_convolution(as_logs(coeff**float(alpha),np),np,guard)
        grouped = degree_power(paired,float(alpha),2,np)
        lower = np.minimum(oldpair,grouped)
        assert np.all(lower <= oldpair) and np.all(lower <= grouped)
        if alpha == 1: assert np.array_equal(oldpair,grouped)
        checks += 3
    return checks+route_mass_checks(guard)


class PairReplay:
    """Selected normalized regions only are cached; FFT tables are discarded."""
    def __init__(self, logs, modulus, calibration, cells, fft, screen, np, guard):
        self.modulus,self.calibration,self.cells,self.fft = modulus,calibration,cells,fft
        self.screen,self.np,self.guard = screen,np,guard
        indices = [(f,s) for f in range(65) for s in range(65-f)]
        self.fd = np.asarray([f for f,s in indices])
        self.sd = np.asarray([s for f,s in indices])
        self.logs = np.asarray([logs[f,s] for f,s in indices])
        self.positive = np.isfinite(self.logs); self.cache = {}

    def marginal(self,full,markers):
        key = (full,markers)
        if key in self.cache: return self.cache[key]
        self.guard(); np,screen,fft = self.np,self.screen,self.fft
        lu,lv = map(self.calibration.logq,markers); single = WEIGHT-2*full
        logs = self.logs+(self.fd*lu+self.sd*lv)[:,None]
        channel_logs = []
        for j in range(6):
            peak = float(np.max(logs[:,j]))
            channel_logs.append(peak+math.log(float(np.sum(np.exp(logs[:,j]-peak)))))
        def normalized(values):
            scale = max(values); entries = tuple(math.exp(x-scale) for x in values)
            if not all(x > 0 for x in entries): raise ArithmeticError("ordinary pair channel underflow")
            return screen.normalized_matrix(entries,scale)
        p,t = normalized(channel_logs[:4]),normalized(channel_logs[4:])
        ordinary = (screen.matrix_power(p,64),screen.matrix_vector(screen.matrix_power(p,63),t))
        denominator = math.log(comb(4096,full)*comb(4096-full,single))+full*lu+single*lv
        result = {"ordinary": {"status":"resolved","region":(*ordinary,denominator)}}
        try:
            shift = float(np.max(logs)); mass = np.exp(logs-shift)
            if np.any(self.positive & (mass == 0)):
                raise ArithmeticError("positive paired FFT input underflow")
            m = self.modulus; degree = np.zeros((m,m,6))
            np.add.at(degree,(self.fd%m,self.sd%m),mass)
            spectra = fft.regional_spectra(degree[:,:,:4].reshape(m,m,2,2),
                degree[:,:,4:].reshape(m,m,2,1),np,self.guard,doublings=6)
            chosen,diagnostics = [],[]
            for (raw,scale),expected in zip(spectra,ordinary):
                scale += 64*shift
                fft.scaled_check(tuple(float(x.real) for x in raw[0,0].ravel()),scale,expected)
                residue,_,_,diagnostic = fft.recover((raw,scale),np)
                thresholds = self.cells.selection_thresholds(fft,residue,raw,np)
                values,minimum = self.cells.select_entries(residue,(full%m,single%m),thresholds)
                chosen.append(screen.normalized_matrix(values,scale))
                diagnostics.append({**diagnostic,"selected_minimum_over_DC":minimum})
            result["residue"] = {"status":"resolved","region":(*chosen,denominator),
                                 "FFT_controls":diagnostics}
        except ArithmeticError as error:
            result["residue"] = {"status":"unresolved","reason":str(error)}
        self.cache[key] = result
        return result

    def profile(self,fulls,markers,source_count):
        entries = [self.marginal(f,m) for f,m in zip(fulls,markers)]
        report = {}
        for name in ("ordinary","residue"):
            selected = [entry[name] for entry in entries]
            if any(x["status"] != "resolved" for x in selected):
                report[name] = {"status":"unresolved","reasons":[x.get("reason") for x in selected]}
            else:
                conditional = self.screen.paired_chains(selected[0]["region"],selected[1]["region"])[0]
                report[name] = {"status":"resolved","log2_bound_FLOAT":(source_count+conditional)/LN2,
                    "log2_conditional_chain_FLOAT":conditional/LN2,
                    "FFT_controls":[x.get("FFT_controls") for x in selected] if name=="residue" else None}
        return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--modulus",type=int,choices=(256,512),default=512)
    args = parser.parse_args()
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator)) > 64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+SECONDS
    small = load_return(); calibration = small.load_calibration(); broad = calibration.load_broad()
    neighbor = broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"pair_neighbor")
    cells = neighbor.load_cells(); fft = cells.load_fft(); screen = fft.load_screen()
    old = screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"pair_old")
    gradient = old.load_gradient(); anchor = gradient.load_anchor(); peak = anchor.peak_memory_reader()
    def guard():
        if time.monotonic() > deadline: raise TimeoutError("180-second pair-pooling deadline; no retry")
        if peak() > MEMORY: raise MemoryError("256MiB pair-pooling budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)

    import numpy as np
    from flint import ctx
    ctx.prec = 256
    tests = exact_tests(np,guard)+fft.exact_toy(np,guard)+cells.exact_scale_tests()
    a,c,s,inverses,weights = small.actual_matrices()
    enumerators = small.pair_enumerators(c,s,inverses,guard)
    field_tests = small.exact_tests(c,s,inverses,weights,enumerators,guard)
    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z,core.ALPHA = args.z,args.alpha
    class LocalVerifier(core.Verifier):
        def guard(self): guard()
    verifier = LocalVerifier(2,deadline)
    _,quiet,old_d,_,_ = core.rational_caps()
    new_d,block_report = small.refined_blocks(args.z,quiet,old_d,enumerators)
    shapes = list(combinations_with_replacement(range(17),4))
    original = verifier.local_caps(shapes); raw,powered = {},{}
    exponent = verifier.exact(args.alpha)
    for occ,row in original.items():
        guard()
        adapter.exact_point(row[2],small.return_cap(occ,quiet,old_d,core.KAPPA),verifier.fmpq)
        refined = list(row)
        refined[2] = verifier.exact(small.return_cap(occ,quiet,new_d,core.KAPPA))
        raw[occ] = tuple(refined)
        powered[occ] = tuple(verifier.positive(x)**exponent for x in refined)
        for x in raw[occ]+powered[occ]:
            value = float(x.mid())
            if not math.isfinite(value) or value < 0 or (x > 0 and value == 0):
                raise ArithmeticError("positive local cap not representable in FLOAT")
    del original
    raw_coeff = old.float_coefficients(raw,np,guard)
    single_coeff = old.float_coefficients(powered,np,guard)
    del raw,powered
    raw_pair = pair_convolution(as_logs(raw_coeff,np),np,guard)
    old_pair = pair_convolution(as_logs(single_coeff,np),np,guard)
    grouped = degree_power(raw_pair,float(args.alpha),64,np)
    selected = np.minimum(old_pair,grouped)
    assert np.all(selected <= old_pair) and np.all(selected <= grouped)
    finite = np.isfinite(old_pair)
    improved = finite & (grouped < old_pair-1e-12)
    controls = {"improved_degree_entries":int(np.count_nonzero(improved)),
        "finite_old_degree_entries":int(np.count_nonzero(finite)),
        "largest_local_gain_bits_FLOAT":float(np.max(old_pair[finite]-selected[finite])/LN2),
        "raw_pair_route_slots":64,"regular_pairs_per_region":64,"terminal_regular_pairs":63}
    progress("local_pair_coefficients_ready",tests=tests,field_tests=field_tests,controls=controls)
    replays = {name:PairReplay(table,args.modulus,calibration,cells,fft,screen,np,guard)
               for name,table in (("old_pair",old_pair),("grouped_pair",grouped),("minimum_pair",selected))}
    single = calibration.Replay(single_coeff,WEIGHT,cells,fft,screen,np,guard)
    counts = screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"pair_counts")
    counts.DEADLINE = deadline
    h,c_h,eta = anchor.source_factor(WEIGHT)
    logsource = float((exponent*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()
                                -core.LIMIT*verifier.exact(args.z).log())).mid())
    cf,cs = F(1562175,990029),F(1140210,880327)
    for fulls in PROFILES:
        guard(); number = counts.profile_count(4096,WEIGHT,*fulls)
        assert number > 0
        markers = tuple(broad.marker_curve(WEIGHT,f,cf,cs) for f in fulls)
        source_count = logsource+math.log(number)
        reference_entries = [single.marginal(f,m,args.modulus) for f,m in zip(fulls,markers)]
        reference = {"status":"unresolved"}
        if all(x["status"]=="resolved" for x in reference_entries):
            chain = screen.paired_chains(reference_entries[0]["region"],reference_entries[1]["region"])[0]
            reference = {"status":"resolved","log2_bound_FLOAT":(source_count+chain)/LN2}
        reports = {name:replay.profile(fulls,markers,source_count) for name,replay in replays.items()}
        old_result = reports["old_pair"]["residue"]
        if reference["status"]==old_result["status"]=="resolved":
            assert abs(reference["log2_bound_FLOAT"]-old_result["log2_bound_FLOAT"]) < 2e-5
        for kind in ("ordinary","residue"):
            before,after = reports["old_pair"][kind],reports["minimum_pair"][kind]
            if before["status"]==after["status"]=="resolved":
                assert after["log2_bound_FLOAT"] <= before["log2_bound_FLOAT"]+2e-5
        progress("selected_profile_complete",fulls=fulls,markers=[[str(x) for x in m] for m in markers],
            log2_exact_source_count=math.log(number)/LN2,
            exact_source_count_sha256=hashlib.sha256(str(number).encode()).hexdigest(),
            frozen_refined_single_packet_control=reference,comparison=reports)
    guard()
    print(json.dumps({"scope":"two actual ordered r1024 profiles; FLOAT pair-pooling A/B only",
        "whole_code_certificate":False,"complete_weight_bound":False,"directed":False,
        "z":str(args.z),"alpha":str(args.alpha),"weight":WEIGHT,"profiles":PROFILES,
        "modulus":args.modulus,"source_factor":{"h":h,"c_h":str(c_h),"eta":str(eta),
            "powered_source_log2_FLOAT":logsource/LN2},"finite_updates":32767,
        "terminal_packets":1,"terminal_factor":"M*T, not T*M or M^2",
        "controls":controls,"tests":tests,"field_tests":field_tests,"refined_blocks":block_report,
        "return_pin":RETURN_PIN,"calibration_pin":small.CALIBRATION_PIN,
        "core_pin":anchor.CORE_PIN,"count_pin":screen.COUNT_PIN,"FFT_pin":cells.FFT_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__ == "__main__":
    main()
