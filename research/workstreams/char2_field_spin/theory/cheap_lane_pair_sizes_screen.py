"""Complete-size FLOAT paired-bank diagnostic at r2048,4096,6144.

Build the frozen genuine ell1/2/3-return and concavity/Jensen pair caps.
Six cyclic doublings give H64 and H63J (J=M*T), with64 common local
scales restored and4096 route slots in either regional normalizer.
Every ordered matching profile is retained in an exact indexed-run-mass
envelope. No BEST queries, profile omissions, or Perron approximation.

For r>4096 enumerate holes k=8192-r. Reverse the COMPLETE local degree
table: hole coefficient(a,s)=actual coefficient(64-a-s,s), identically
for H and J. This changes no state, time, cap, or refresh distribution.
The source factor remains Omega_r, never Omega_k. The target normalizer
is M4096(a,k-2a), equal to the complementary actual multinomial.

The inherited fifteen-center curve is not optimized. Default centers
exclude endpoints. Explicit endpoint centers use a half-slot in place
of a zero numerator/denominator, as a positive witness, not a limit.
Unresolved residues use an explicitly reported ordinary Cauchy fallback;
no coefficient is clipped. Resolved and fallback conditional leaders
are reported separately, without pretending they have exact counts.

This is FLOAT reconnaissance, not a certificate. One weight per run;
180/300s,256/512MiB,20GiB free, stdout only, no generated files.
"""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from math import comb
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
BANK_PIN="2b5099dc1794a97a2c72c4d09584a62ee248b4982b9889bc791e01e5bd685742"
BROAD_PIN="df955f7f8458b0e22314c82311b2b713c2142cd9c24b1d188441edf1a208c33f"
DEFAULT_Z={2048:F(1,4),4096:F(1,16),6144:F(1,64)}
ALPHA=F(119,1280)
CF,CS=F(1562175,990029),F(1140210,880327)
LN2=math.log(2)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def markers(weight,center):
    """Positive fixed witnesses; only exact zeros receive half-slot smoothing."""
    empty=4096-weight+center
    assert empty>=0 and 0<=2*center<=weight
    den=F(empty) if empty else F(1,2)
    full=F(center) if center else F(1,2)
    single=F(weight-2*center) if 2*center<weight else F(1,2)
    return CF*full/den,CS*single/den


def reverse_degrees(table,slots):
    """Pure coefficient relabeling; leaves channels and temporal order intact."""
    return {(slots-f-s,s):row for (f,s),row in table.items()}


def complement_tests(broad,guard):
    checks=words=0
    for n in (6,8,10):
        by_weight=defaultdict(Counter); m=n//2
        for mask in range(1<<n):
            guard(); bits=[(mask>>i)&1 for i in range(n)]
            r=sum(bits)
            f0=sum(bits[2*i]*bits[2*i+1] for i in range(m))
            f1=sum(bits[2*i+1]*bits[(2*i+2)%n] for i in range(m))
            c=sum(bits[i] and not bits[(i-1)%n] for i in range(n))
            holes=[1-x for x in bits]
            a0=sum(holes[2*i]*holes[2*i+1] for i in range(m))
            a1=sum(holes[2*i+1]*holes[(2*i+2)%n] for i in range(m))
            assert (a0,a1)==(m-r+f0,m-r+f1)
            assert c==sum(holes[i] and not holes[(i-1)%n] for i in range(n))
            assert r-2*f0==n-r-2*a0
            s=r-2*f0
            assert comb(m,f0)*comb(m-f0,s)==comb(m,a0)*comb(m-a0,s)
            by_weight[r][c]+=1; words+=1; checks+=4
        for r in range(1,n):
            expected=dict(broad.run_masses(n,r))
            assert dict(by_weight[r])==expected
            assert expected==dict(broad.run_masses(n,n-r))
            assert sum(expected.values())==comb(n,r)
            checks+=3
    # Exact noncommuting MM/MT products commute with degree reversal.
    def mul(a,b):
        # Each coefficient is a rectangular tuple of rows.
        result={}
        for (f,s),x in a.items():
            for (g,t),y in b.items():
                key=(f+g,s+t)
                if key not in result: result[key]=[[F(0) for _ in y[0]] for _ in x]
                for i in range(len(x)):
                    for j in range(len(y[0])):
                        result[key][i][j]+=sum(x[i][k]*y[k][j] for k in range(len(y)))
        return {key:tuple(map(tuple,row)) for key,row in result.items()}
    for slots in (1,2):
        m={}; t={}
        for f in range(slots+1):
            for s in range(slots-f+1):
                mass=comb(slots,f)*comb(slots-f,s)
                m[f,s]=((F(mass*(1+f),7),F(mass*(2+s),11)),
                        (F(mass*(3+s),13),F(mass*(5+f),17)))
                t[f,s]=((F(mass*(7+f),19),),(F(mass*(11+s),23),))
                assert mass==comb(slots,slots-f-s)*comb(f+s,s)
                checks+=1
        rm,rt=reverse_degrees(m,slots),reverse_degrees(t,slots)
        assert reverse_degrees(rm,slots)==m
        assert reverse_degrees(mul(m,m),2*slots)==mul(rm,rm)
        assert reverse_degrees(mul(m,t),2*slots)==mul(rm,rt)
        checks+=3
    return {"exact_complement_checks":checks,"exhaustive_cycle_words":words,
            "noncommuting_regular_terminal_reversal":True}


def marginals(log_coefficients,actual_weight,centers,modulus,ordinary_fallback,
              broad,cells,fft,screen,np,guard,progress):
    route_weight=min(actual_weight,8192-actual_weight); holes=actual_weight>4096
    assignment=[min(centers,key=lambda c:(abs(f-c),c)) for f in range(route_weight//2+1)]
    indices=[(f,s) for f in range(65) for s in range(65-f)]
    fd=np.asarray([f for f,s in indices]); sd=np.asarray([s for f,s in indices])
    # Reverse before any residue folding: this includes every local degree.
    local={(f,s):log_coefficients[f,s] for f,s in indices}
    if holes: local=reverse_degrees(local,64)
    loglocal=np.asarray([local[f,s] for f,s in indices]); del local
    if np.isnan(loglocal).any() or np.isposinf(loglocal).any():
        raise ArithmeticError("invalid pair log coefficient")
    positive=np.isfinite(loglocal); result=[None]*len(assignment)
    reports=[]; fallbacks=[]
    stats={"max_imaginary_over_DC":0.0,"max_negative_over_DC":0.0,"min_selected_over_DC":1.0}
    for center in centers:
        guard(); targets=[f for f,c in enumerate(assignment) if c==center]
        if not targets: continue
        u,v=markers(route_weight,center); lu,lv=broad.logq(u),broad.logq(v)
        logs=loglocal+(fd*lu+sd*lv)[:,None]
        shift=float(np.max(logs)); mass=np.exp(logs-shift)
        channel_logs=[]
        for j in range(6):
            peak=float(np.max(logs[:,j]))
            channel_logs.append(peak+math.log(float(np.sum(np.exp(logs[:,j]-peak)))))
        def from_logs(values):
            scale=max(values); entries=tuple(math.exp(x-scale) for x in values)
            if not all(math.isfinite(x) and x>0 for x in entries):
                raise ArithmeticError("ordinary positive pair channel underflow")
            return screen.normalized_matrix(entries,scale)
        op,ot=from_logs(channel_logs[:4]),from_logs(channel_logs[4:])
        ordinary=(screen.matrix_power(op,64),screen.matrix_vector(screen.matrix_power(op,63),ot))
        tables=[]; cell_error=None
        spectra=degree=p=t=residue=raw=None
        try:
            if np.any(positive & (mass==0)):
                raise ArithmeticError("positive paired route term underflowed in FFT input")
            degree=np.zeros((modulus,modulus,6))
            np.add.at(degree,(fd%modulus,sd%modulus),mass)
            p=degree[:,:,:4].reshape(modulus,modulus,2,2)
            t=degree[:,:,4:].reshape(modulus,modulus,2,1)
            spectra=fft.regional_spectra(p,t,np,guard,doublings=6)
            for (raw,scale),expected in zip(spectra,ordinary):
                scale+=64*shift
                fft.scaled_check(tuple(float(x.real) for x in raw[0,0].ravel()),scale,expected)
                residue,_,_,diagnostic=fft.recover((raw,scale),np)
                tables.append((residue,scale,cells.selection_thresholds(fft,residue,raw,np)))
                for key in ("max_imaginary_over_DC","max_negative_over_DC"):
                    stats[key]=max(stats[key],diagnostic[key])
        except ArithmeticError as error:
            if not ordinary_fallback: raise ArithmeticError(f"center{center}: {error}") from error
            cell_error=str(error)
        for f in targets:
            guard(); s=route_weight-2*f; error=cell_error; projected=[]
            if error is None:
                try:
                    for residue,scale,thresholds in tables:
                        values,minimum=cells.select_entries(residue,(f%modulus,s%modulus),thresholds)
                        projected.append(screen.normalized_matrix(values,scale))
                        stats["min_selected_over_DC"]=min(stats["min_selected_over_DC"],minimum)
                except ArithmeticError as failure:
                    if not ordinary_fallback: raise ArithmeticError(f"f{f},s{s}: {failure}") from failure
                    error=str(failure)
            if error is not None:
                projected=ordinary
                fallbacks.append({"coordinate_full":f,"actual_full":f+actual_weight-4096 if holes else f,
                    "single":s,"center":center,"reason":error})
            denominator=math.log(comb(4096,f)*comb(4096-f,s))+f*lu+s*lv
            result[f]=(projected[0],projected[1],denominator,u,v)
        report={"center":center,"u":str(u),"v":str(v),"first_coordinate_full":targets[0],
            "last_coordinate_full":targets[-1],"marginals":len(targets),
            "ordinary_fallbacks":sum(x["center"]==center for x in fallbacks)}
        reports.append(report); progress("pair_size_marker_cell_complete",**report)
        tables.clear(); del spectra,degree,p,t,residue,raw
    assert all(x is not None for x in result)
    return result,reports,stats,fallbacks


def complete_join(actual_weight,regions,fallbacks,broad,screen,np,log_source,guard,progress):
    """Exact run masses times maxima over all ordered finite chains."""
    k=min(actual_weight,8192-actual_weight); maximum=k//2
    regular=np.asarray([x[0][0] for x in regions]); terminal=np.asarray([x[1][0] for x in regions])
    rl=np.asarray([x[0][1]-x[2] for x in regions]); tl=np.asarray([x[1][1]-x[2] for x in regions])
    fallback=np.zeros(maximum+1,dtype=bool)
    for row in fallbacks: fallback[row["coordinate_full"]]=True
    offset=actual_weight-4096 if actual_weight>4096 else 0
    total=-math.inf; categorized=[-math.inf,-math.inf]; ordered=source_total=0
    category_counts=[0,0]; leaders=[[],[]]; largest=[]
    for c,number in broad.run_masses(8192,k):
        guard(); first=np.arange(max(0,k-c-maximum),min(maximum,k-c)+1,dtype=int)
        second=k-c-first
        chain=broad.batched_chains(first,second,regular,rl,terminal,tl,np)
        if not np.isfinite(chain).all(): raise ArithmeticError("nonfinite complete conditional chain")
        flags=fallback[first]|fallback[second]
        base=log_source+math.log(number); best=int(np.argmax(chain))
        term=base+float(chain[best]); total=screen.logadd(total,term)
        row={"runs":c,"log2_run_envelope_FLOAT":term/LN2,
             "actual_full":[int(first[best])+offset,int(second[best])+offset],
             "coordinate_full":[int(first[best]),int(second[best])],
             "uses_ordinary_fallback":bool(flags[best]),"ordered_profiles":len(first)}
        largest.append(row); largest=sorted(largest,key=lambda x:x["log2_run_envelope_FLOAT"],reverse=True)[:8]
        for category in (0,1):
            positions=np.flatnonzero(flags==bool(category)); category_counts[category]+=len(positions)
            if len(positions):
                p=int(positions[int(np.argmax(chain[positions]))]); upper=base+float(chain[p])
                categorized[category]=screen.logadd(categorized[category],upper)
                leaders[category].append({"runs":c,"actual_full":[int(first[p])+offset,int(second[p])+offset],
                    "coordinate_full":[int(first[p]),int(second[p])],
                    "single":[k-2*int(first[p]),k-2*int(second[p])],
                    "log2_full_run_mass_charged_FLOAT":upper/LN2,
                    "log2_conditional_chain_FLOAT":float(chain[p])/LN2})
                leaders[category]=sorted(leaders[category],key=lambda x:x["log2_full_run_mass_charged_FLOAT"],reverse=True)[:6]
        ordered+=len(first); source_total+=number
        if c%max(32,k//16)==0 or c==k:
            progress("pair_size_run_envelope",through_runs=c,ordered_profiles=ordered)
    assert source_total==comb(8192,actual_weight)==comb(8192,k)
    assert ordered==(maximum+1)**2-1 and sum(category_counts)==ordered
    return {"actual_weight":actual_weight,"coordinate_weight":k,"complete_source_weight":True,
        "ordered_profiles":ordered,"log2_run_envelope_FLOAT":total/LN2,"exact_count_queries":0,
        "largest_run_bounds":largest,"source_log2_FLOAT":log_source/LN2,
        "method_categories":{name:{"ordered_profiles":category_counts[i],
            "log2_separate_run_envelope_FLOAT":categorized[i]/LN2 if math.isfinite(categorized[i]) else None,
            "leaders":leaders[i]} for i,name in enumerate(("resolved_residue","uses_ordinary_fallback"))},
        "category_scope":"Each category charges the whole run mass; these are separate diagnostic upper bounds, not a partition of source counts.",
        "checks":{"all_indexed_run_masses_sum_to_binomial":True,"all_ordered_chains_retained":True,
            "both_terminal_orientations":True,"updates":32767,"terminal_slots":32,
            "source_factor_applied_once_for_actual_weight":True,"counts_outside_alpha":True}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight",type=int,choices=(2048,4096,6144),default=2048)
    parser.add_argument("--z",type=F)
    parser.add_argument("--alpha",type=F,default=ALPHA)
    parser.add_argument("--modulus",type=int,choices=(128,256,512),default=256)
    parser.add_argument("--centers")
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--memory-mib",type=int,choices=(256,512),default=256)
    parser.add_argument("--ordinary-fallback",action=argparse.BooleanOptionalAction,default=True)
    args=parser.parse_args(); r=args.weight; k=min(r,8192-r)
    z=DEFAULT_Z[r] if args.z is None else args.z
    if not 0<z<1 or not 0<args.alpha<=1: parser.error("require0<z<1,0<alpha<=1")
    if max(abs(n).bit_length() for x in (z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    bank=load("cheap_lane_pair_bank_screen.py",BANK_PIN,"sizes_frozen_bank")
    broad=load("cheap_lane_broad_peeling_screen.py",BROAD_PIN,"sizes_frozen_broad")
    try: centers=tuple(broad.default_centers(k)) if args.centers is None else tuple(sorted(map(int,args.centers.split(","))))
    except ValueError: parser.error("invalid comma-separated centers")
    if not 1<=len(centers)<=24 or len(set(centers))!=len(centers) or not 0<=centers[0]<=centers[-1]<=k//2:
        parser.error("choose1..24 distinct centers in0..floor(coordinate_weight/2)")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"sizes_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"sizes_triple")
    pair=correction.load_pair()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"sizes_neighbor")
    peel=broad.load("cheap_lane_profile_peeling.py",broad.PEEL_PIN,"sizes_peel_tests")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"sizes_old_coefficients")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("size-screen deadline; no complete bound or automatic retry")
        if peak()>args.memory_mib*(1<<20): raise MemoryError("size-screen memory budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,"peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    tests=broad.exact_tests(np,screen,peel)
    tests.update(complement_tests(broad,guard))
    tests["pair_and_correction_checks"]=correction.exact_tests(pair,np,guard)+pair.exact_tests(np,guard)
    tests["FFT_exact_toy_coefficients"]=fft.exact_toy(np,guard)
    tests["scale_checks"]=cells.exact_scale_tests()
    caps=composition.build_caps(z,args.alpha,triple,anchor,guard,deadline)
    verifier,core=caps["verifier"],caps["core"]
    tables,single,controls=correction.corrected_pair_logs(caps["powered"],float(args.alpha),
        tuple(range(-16,17,4)),pair,old,np,guard,raw=caps["raw"])
    del caps["raw"],caps["powered"],single
    selected=tables["concavity_Jensen_minimum"]; del tables
    progress("pair_size_coefficients_ready",actual_weight=r,coordinate_weight=k,holes=r>4096,
        z=str(z),alpha=str(args.alpha),centers=centers,tests=tests,controls=controls)
    regions,bank_report,stats,fallbacks=marginals(selected,r,centers,args.modulus,
        args.ordinary_fallback,broad,cells,fft,screen,np,guard,progress)
    del selected
    missing_group_limit,c_h,eta=anchor.source_factor(r); exponent=verifier.exact(args.alpha)
    log_source=float((exponent*(verifier.exact(eta).log()+r*verifier.exact(core.BETA).log()
        -core.LIMIT*verifier.exact(z).log())).mid())
    progress("pair_size_marginals_ready",marginals=len(regions),ordinary_fallbacks=len(fallbacks))
    report=complete_join(r,regions,fallbacks,broad,screen,np,log_source,guard,progress)
    guard()
    print(json.dumps({"scope":"complete selected actual source weight, FLOAT paired-bank run-mass envelope",
        "directed":False,"negative_is_not_a_certificate":True,"whole_code_certificate":False,
        "other_weights_covered":False,"parameter_optimization":False,
        "parameters":{"actual_weight":r,"coordinate_weight":k,"hole_coordinates":r>4096,
            "z":str(z),"alpha":str(args.alpha),"modulus":args.modulus,"centers":centers,
            "curve_full":str(CF),"curve_single":str(CS),"endpoint_zero_replacement":"1/2 slot",
            "seconds":args.seconds,"memory_budget_mib":args.memory_mib},
        "result":report,"bank":bank_report,"FFT_controls":stats,
        "ordinary_fallback":{"explicitly_enabled":args.ordinary_fallback,"marginals_used":len(fallbacks),
            "first16":fallbacks[:16],"all_coordinate_full_indices":[x["coordinate_full"] for x in fallbacks],
            "method":"ordinary positive Cauchy numerator, no coefficient clipping"},
        "joint_missing_for_actual_weight":{"h":missing_group_limit,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "controls":controls,"tests":tests,"field_tests":caps["tests"],"block_caps":caps["block_report"],
        "bank_pin":BANK_PIN,"broad_pin":BROAD_PIN,"composition_pin":bank.COMPOSITION_PIN,
        "correction_pin":composition.CORRECTION_PIN,"triple_pin":composition.TRIPLE_PIN,
        "pair_pin":correction.PAIR_PIN,"FFT_pin":cells.FFT_PIN,"core_pin":anchor.CORE_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
