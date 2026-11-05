"""Compose the exact ell1/2/3 return caps with binary-pair corrections.

This selected-profile FLOAT replay pins both independent helpers. Build
the genuine raw ell3 return cap, then raise every local cap to alpha once.
Pass these maps to corrected_pair_logs, with raw supplied for the optional
degree-Jensen minimum. Labels distinguish old, correction-only, Jensen,
and correction-plus-Jensen stages. The frozen single-packet ell3 engine
is an independent normalization/terminal control for the old-pair stage.

Only actual ordered r1024 profiles(384,500) and(384,384) are evaluated.
Their exact counts, own-f markers, source factor, alpha, and output-only
terminal are shared across all stages. A negative FLOAT result is not a
certificate, a complete r1024 bound, or a whole-code proof.
No optimizer or automatic scans.180s/256MiB/20GiB; stdout only.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement
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
CORRECTION_PIN = "804eeb5b9e7df90cba0ffa57651b0e411784e5bb775e6e3d62b532e994728cf1"
TRIPLE_PIN = "d03f60009fcc9f0f757a19681e9647609aad7f501493b1d89980fca8bf1da6b5"
WEIGHT, SECONDS, MEMORY = 1024, 180, 256*(1 << 20)
PROFILES = ((384,500),(384,384))
LN2 = math.log(2)


def load(name,digest,module_name):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("frozen helper hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(module_name,path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def build_caps(z,alpha,triple,anchor,guard,deadline):
    """Return genuine ell1/2/3 raw+powered maps and their exact provenance.

    z and alpha are Fractions. The caller owns its resource guard/deadline
    and chooses the subsequent route algorithm. No route matrix, source
    factor, profile count, or fractional power beyond the single local
    power is computed here. The imported verifier is private to this call.
    """
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError("require0<z<1 and0<alpha<=1")
    small=triple.load_small()
    a,c,s,inverses,weights=small.actual_matrices()
    pairs=small.pair_enumerators(c,s,inverses,guard)
    tests=small.exact_tests(c,s,inverses,weights,pairs,guard)
    records,triple_tests=triple.triples(c,s,inverses,small,guard)
    tests["triple_tests"]=triple_tests
    adapter=anchor.load_adapter(); core=adapter.load_verifier()
    core.Z,core.ALPHA=z,alpha
    class LocalVerifier(core.Verifier):
        def guard(self): guard()
    verifier=LocalVerifier(2,deadline)
    _,quiet,old_d,_,_=core.rational_caps()
    refined_d,block_report=small.refined_blocks(z,quiet,old_d,pairs)
    assert quiet[3] == F(1,256**3)
    values=[(points,small.polynomial(histogram,z)/256**3) for points,histogram in records]
    maximum=max(value for _,value in values)
    before=refined_d[3]; refined_d=list(refined_d)
    refined_d[3]=min(before,maximum)
    block_report.update({"ell1_2_baseline_D3":str(before),"refined_D3":str(refined_d[3]),
        "exact_uniform_D3":str(maximum),
        "D3_maximizing_triples":[points for points,value in values if value==maximum]})
    assert all(0 <= x <= y for x,y in zip(refined_d,old_d))
    original=verifier.local_caps(list(combinations_with_replacement(range(17),4)))
    raw,powered={},{}; exponent=verifier.exact(alpha)
    for occ,row in original.items():
        guard()
        adapter.exact_point(row[2],small.return_cap(occ,quiet,old_d,core.KAPPA),verifier.fmpq)
        refined=list(row)
        refined[2]=verifier.exact(small.return_cap(occ,quiet,refined_d,core.KAPPA))
        raw[occ]=tuple(refined)
        powered[occ]=tuple(verifier.positive(x)**exponent for x in refined)
        for value in raw[occ]+powered[occ]:
            approximate=float(value.mid())
            if not math.isfinite(approximate) or approximate<0 or (value>0 and approximate==0):
                raise ArithmeticError("positive local cap not representable in FLOAT")
    del original
    return {"raw":raw,"powered":powered,"verifier":verifier,"core":core,
        "small":small,"tests":tests,"block_report":block_report}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--modulus",type=int,choices=(256,512),default=512)
    parser.add_argument("--rho-radius",type=int,choices=(4,8,12,16,20,24,28,32),default=16)
    args=parser.parse_args()
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    correction=load("cheap_lane_pair_concavity.py",CORRECTION_PIN,"triple_concavity")
    triple=load("cheap_lane_exact_triple_return.py",TRIPLE_PIN,"triple_return_caps")
    pair=correction.load_pair(); small=triple.load_small()
    calibration=small.load_calibration(); broad=calibration.load_broad()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"combined_neighbor")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"combined_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("180-second composition deadline; no retry")
        if peak()>MEMORY: raise MemoryError("256MiB composition budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    tests=correction.exact_tests(pair,np,guard)+pair.exact_tests(np,guard)
    caps=build_caps(args.z,args.alpha,triple,anchor,guard,deadline)
    verifier,core=caps["verifier"],caps["core"]
    progress("exact_return_caps_ready",tests=tests,field_tests=caps["tests"])
    exponents=tuple(range(-args.rho_radius,args.rho_radius+1,4))
    tables,single,controls=correction.corrected_pair_logs(caps["powered"],float(args.alpha),
        exponents,pair,old,np,guard,raw=caps["raw"])
    del caps["raw"],caps["powered"]
    progress("corrected_pair_coefficients_ready",controls=controls)
    replays={name:pair.PairReplay(logs,args.modulus,calibration,cells,fft,screen,np,guard)
             for name,logs in tables.items()}
    reference=calibration.Replay(single,WEIGHT,cells,fft,screen,np,guard)
    counts=screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"combined_counts")
    counts.DEADLINE=deadline
    h,c_h,eta=anchor.source_factor(WEIGHT); exponent=verifier.exact(args.alpha)
    logsource=float((exponent*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()
                              -core.LIMIT*verifier.exact(args.z).log())).mid())
    cf,cs=F(1562175,990029),F(1140210,880327)
    results=[]
    for fulls in PROFILES:
        guard(); number=counts.profile_count(4096,WEIGHT,*fulls); assert number>0
        markers=tuple(broad.marker_curve(WEIGHT,f,cf,cs) for f in fulls)
        source_count=logsource+math.log(number)
        reports={"ell1_2_3_"+name:replay.profile(fulls,markers,source_count)
                 for name,replay in replays.items()}
        refs=[reference.marginal(f,m,args.modulus) for f,m in zip(fulls,markers)]
        control={"status":"unresolved"}
        if all(x["status"]=="resolved" for x in refs):
            score=(source_count+screen.paired_chains(refs[0]["region"],refs[1]["region"])[0])/LN2
            control={"status":"resolved","log2_bound_FLOAT":score}
            old_result=reports["ell1_2_3_old_pair"]["residue"]
            if old_result["status"]=="resolved": assert abs(score-old_result["log2_bound_FLOAT"])<2e-5
        for kind in ("ordinary","residue"):
            before=reports["ell1_2_3_old_pair"][kind]
            for label in ("concavity_pair","old_Jensen_minimum","concavity_Jensen_minimum"):
                after=reports["ell1_2_3_"+label][kind]
                if before["status"]==after["status"]=="resolved":
                    assert after["log2_bound_FLOAT"]<=before["log2_bound_FLOAT"]+2e-5
        row={"fulls":fulls,"markers":[[str(x) for x in m] for m in markers],
            "exact_count_sha256":hashlib.sha256(str(number).encode()).hexdigest(),
            "log2_exact_count":math.log(number)/LN2,
            "ell1_2_3_single_packet_control":control,"comparison":reports}
        results.append(row); progress("selected_profile_complete",**row)
    guard()
    print(json.dumps({"scope":"two selected ordered r1024 profiles; exact ell3 caps plus FLOAT pair correction",
        "directed":False,"whole_code_certificate":False,"complete_weight_bound":False,
        "weight":WEIGHT,"z":str(args.z),"alpha":str(args.alpha),"modulus":args.modulus,
        "updates":32767,"terminal_outputs":1,"terminal_pair":"M*T",
        "source_factor":{"h":h,"c_h":str(c_h),"eta":str(eta),"log2_powered_source_FLOAT":logsource/LN2},
        "controls":controls,"tests":tests,"field_tests":caps["tests"],"block_caps":caps["block_report"],
        "results":results,"correction_pin":CORRECTION_PIN,"triple_pin":TRIPLE_PIN,
        "pair_pin":correction.PAIR_PIN,"count_pin":screen.COUNT_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__ == "__main__":
    main()
