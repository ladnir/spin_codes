"""Complete r1024 FLOAT profile scan with exact-small-return pair caps.

Pin the ell1/2/3 composition and use its concavity/Jensen minimum. Local
coefficients have64 route slots. Six positive cyclic matrix doublings
produce H^64 and H^63*J, where H represents two regular packets and J
represents M*T. Both regional objects contain64 paired factors, so their
common local scale is restored64 times. The route denominator remains
C(4096,f)C(4096-f,s)u^f v^s, with s=1024-2f.

Every f=0..512 is assigned to the nearest fixed bank center (ties lower).
The default bank is the frozen scaled bank plus0 and512. Each target
retains its own denominator; markers are never optimized. An independent
frozen PairReplay at f384 checks one shared-center marginal.

The frozen broad join evaluates ALL263168 ordered profiles with the
honest finite terminal. Exact indexed counts are queried lazily; after a
query/time budget expires, the full residual run mass remains under its
largest remaining conditional cap. No source profile is omitted.
Instrumentation separately reports queried weighted leaders whose
marginals are resolved and those using an ordinary fallback.

Ordinary fallback is explicitly ENABLED by default. It can be disabled
with --no-ordinary-fallback. Unresolved residues are not clipped or called
coefficient estimates. This is a complete FLOAT upper-envelope diagnostic,
not a directed certificate or a claim about any other source weight.
225/300s;256MiB;20GiB free; stdout only; no generated files or frozen edits.
"""

import argparse
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

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
COMPOSITION_PIN = "57c3abe184ce6ee643779c48106b79f1d39c85bfb4175e28404ec86b88bc7169"
WEIGHT, MEMORY, LN2 = 1024, 256*(1 << 20), math.log(2)
CURVE_FULL, CURVE_SINGLE = F(1562175,990029), F(1140210,880327)


def load_composition():
    path = Path(__file__).with_name("cheap_lane_pair_concavity_triple.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != COMPOSITION_PIN:
        raise RuntimeError("frozen pair-concavity-triple hash mismatch")
    spec = importlib.util.spec_from_file_location("pair_bank_composition",path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def marginals(log_coefficients,centers,modulus,ordinary_fallback,
              broad,cells,fft,screen,np,guard,progress):
    """Bulk paired adaptation of frozen broad.marginals, without alpha."""
    assignment = [min(centers,key=lambda c:(abs(f-c),c)) for f in range(513)]
    indices = [(f,s) for f in range(65) for s in range(65-f)]
    fd = np.asarray([f for f,s in indices]); sd = np.asarray([s for f,s in indices])
    loglocal = np.asarray([log_coefficients[f,s] for f,s in indices])
    if np.isnan(loglocal).any() or np.isposinf(loglocal).any():
        raise ArithmeticError("invalid positive pair log coefficient")
    positive = np.isfinite(loglocal)
    result = [None]*513; reports=[]; fallbacks=[]
    stats = {"max_imaginary_over_DC":0.0,"max_negative_over_DC":0.0,
             "min_selected_over_DC":1.0}
    for center in centers:
        guard()
        targets = [f for f,c in enumerate(assignment) if c==center]
        if not targets: continue
        u,v = broad.marker_curve(WEIGHT,center,CURVE_FULL,CURVE_SINGLE)
        lu,lv = broad.logq(u),broad.logq(v)
        logs = loglocal+(fd*lu+sd*lv)[:,None]
        shift = float(np.max(logs)); mass = np.exp(logs-shift)
        channel_logs=[]
        for j in range(6):
            peak = float(np.max(logs[:,j]))
            channel_logs.append(peak+math.log(float(np.sum(np.exp(logs[:,j]-peak)))))
        def from_logs(values):
            scale=max(values); entries=tuple(math.exp(x-scale) for x in values)
            if not all(x>0 for x in entries):
                raise ArithmeticError("ordinary paired channel underflow")
            return screen.normalized_matrix(entries,scale)
        op,ot = from_logs(channel_logs[:4]),from_logs(channel_logs[4:])
        ordinary = (screen.matrix_power(op,64),screen.matrix_vector(screen.matrix_power(op,63),ot))
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
                scale += 64*shift
                fft.scaled_check(tuple(float(x.real) for x in raw[0,0].ravel()),scale,expected)
                residue,_,_,diagnostic=fft.recover((raw,scale),np)
                thresholds=cells.selection_thresholds(fft,residue,raw,np)
                tables.append((residue,scale,thresholds))
                for key in ("max_imaginary_over_DC","max_negative_over_DC"):
                    stats[key]=max(stats[key],diagnostic[key])
        except ArithmeticError as error:
            if not ordinary_fallback: raise ArithmeticError(f"center{center}: {error}") from error
            cell_error=str(error)
        for f in targets:
            guard(); s=WEIGHT-2*f; error=cell_error; projected=[]
            if error is None:
                try:
                    for residue,scale,thresholds in tables:
                        values,minimum=cells.select_entries(residue,(f%modulus,s%modulus),thresholds)
                        projected.append(screen.normalized_matrix(values,scale))
                        stats["min_selected_over_DC"]=min(stats["min_selected_over_DC"],minimum)
                except ArithmeticError as failure:
                    if not ordinary_fallback:
                        raise ArithmeticError(f"center{center},f{f},s{s}: {failure}") from failure
                    error=str(failure)
            if error is not None:
                projected=ordinary
                fallbacks.append({"full":f,"single":s,"center":center,"reason":error})
            denominator=math.log(comb(4096,f)*comb(4096-f,s))+f*lu+s*lv
            result[f]=(projected[0],projected[1],denominator,u,v)
        report={"center":center,"u":str(u),"v":str(v),"first_full":targets[0],
                "last_full":targets[-1],"marginals":len(targets),
                "ordinary_fallbacks":sum(x["center"]==center for x in fallbacks)}
        reports.append(report); progress("pair_marker_cell_complete",**report)
        tables.clear()
        del spectra,degree,p,t,residue,raw
    assert all(x is not None for x in result)
    return result,reports,stats,fallbacks,assignment


def join_with_diagnostics(regions,fallbacks,counts,broad,screen,np,log_source,
                          queries,count_seconds,deadline,guard,progress):
    """Only instrumentation: the frozen join's returned bound is unchanged."""
    fallback_set={x["full"] for x in fallbacks}
    original_chain=broad.batched_chains; original_count=counts.profile_count
    current={}; leaders={"resolved_marginals":[],"uses_ordinary_fallback":[]}
    seen=set(); queries_observed=0
    def capture_chains(f0,f1,*args,**kwargs):
        chain=original_chain(f0,f1,*args,**kwargs)
        current.clear()
        current.update({(int(a),int(b)):float(value) for a,b,value in zip(f0,f1,chain)})
        return chain
    def capture_count(pairs,r,f0,f1):
        nonlocal queries_observed
        value=original_count(pairs,r,f0,f1)
        queries_observed+=1
        if value:
            for a,b in {(f0,f1),(f1,f0)}:
                key=(a,b)
                if key in seen: continue
                seen.add(key)
                conditional=current[key]
                category="uses_ordinary_fallback" if a in fallback_set or b in fallback_set else "resolved_marginals"
                row={"full":[a,b],"single":[r-2*a,r-2*b],"runs":r-a-b,
                     "log2_weighted_profile_FLOAT":(log_source+math.log(value)+conditional)/LN2,
                     "log2_exact_count":math.log(value)/LN2,
                     "log2_conditional_chain_FLOAT":conditional/LN2,
                     "fallback_full":[f for f in (a,b) if f in fallback_set]}
                leaders[category].append(row)
                leaders[category]=sorted(leaders[category],key=lambda x:x["log2_weighted_profile_FLOAT"],reverse=True)[:6]
        return value
    broad.batched_chains=capture_chains; counts.profile_count=capture_count
    try:
        result=broad.join(WEIGHT,regions,counts,screen,np,log_source,queries,count_seconds,
                          deadline,guard,progress)
    finally:
        broad.batched_chains=original_chain; counts.profile_count=original_count
    assert queries_observed==result["completed_unique_exact_queries"]
    for row in result["largest_run_bounds"]:
        profile=row["residual_maximizing_full"]
        row["residual_maximizer_uses_fallback"]=(any(f in fallback_set for f in profile)
                                                if profile is not None else None)
    return result,{"queried_weighted_leaders_by_numerical_method":leaders,
        "query_results_observed":queries_observed,
        "classification_scope":"queried exact-count profiles only; all unqueried sources remain in the reported tails",
        "no_additional_count_queries":True,"bound_algorithm_unchanged":True}


def marginal_control(logs,regions,assignment,fallbacks,modulus,pair,calibration,
                     broad,cells,fft,screen,np,guard):
    full=384; center=assignment[full]
    markers=broad.marker_curve(WEIGHT,center,CURVE_FULL,CURVE_SINGLE)
    replay=pair.PairReplay(logs,modulus,calibration,cells,fft,screen,np,guard)
    reference=replay.marginal(full,markers)
    used_fallback=any(x["full"]==full for x in fallbacks)
    kind="ordinary" if used_fallback else "residue"
    if reference[kind]["status"]!="resolved":
        raise ArithmeticError("independent shared-center marginal control unresolved")
    expected=reference[kind]["region"]; actual=regions[full]
    differences=[]
    for (a,la),(b,lb) in zip(actual[:2],expected[:2]):
        for x,y in zip(a,b):
            differences.append(abs(math.log(x)+la-math.log(y)-lb))
    assert abs(actual[2]-expected[2])<1e-10
    assert max(differences)<2e-8
    left=screen.paired_chains(actual,actual)[0]
    right=screen.paired_chains(expected,expected)[0]
    assert abs(left-right)/LN2<2e-5
    return {"full":full,"single":WEIGHT-2*full,"shared_center":center,
        "markers":list(map(str,markers)),"method":kind,
        "maximum_entry_log_difference_FLOAT":max(differences),
        "finite_ordered_chain_log2_difference_FLOAT":abs(left-right)/LN2,
        "normalizer_matches":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--modulus",type=int,choices=(256,512),default=512)
    parser.add_argument("--centers")
    parser.add_argument("--seconds",type=int,choices=(225,300),default=300)
    parser.add_argument("--count-seconds",type=float,default=30)
    parser.add_argument("--max-queries",type=int,default=20000)
    parser.add_argument("--ordinary-fallback",action=argparse.BooleanOptionalAction,default=True)
    args=parser.parse_args()
    if not 0<args.z<1 or not 0<args.alpha<=1:
        parser.error("require0<z<1 and0<alpha<=1")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64 bits")
    if not 0<=args.count_seconds<=min(120,args.seconds-25) or not 0<=args.max_queries<=100000:
        parser.error("count-seconds0..120; max-queries0..100000")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    composition=load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"bank_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"bank_triple")
    pair=correction.load_pair(); small=triple.load_small()
    calibration=small.load_calibration(); broad=calibration.load_broad()
    try:
        centers=(tuple(sorted(set(broad.default_centers(WEIGHT))|{0,512})) if args.centers is None
                 else tuple(sorted(map(int,args.centers.split(",")))))
    except ValueError: parser.error("invalid comma-separated center list")
    if not 1<=len(centers)<=24 or len(set(centers))!=len(centers) or not 0<=centers[0]<=centers[-1]<=512:
        parser.error("choose1..24 distinct centers in0..512")
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"pair_bank_neighbor")
    peel=broad.load("cheap_lane_profile_peeling.py",broad.PEEL_PIN,"pair_bank_peel_tests")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"pair_bank_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("pair-bank deadline; no complete bound or retry")
        if peak()>MEMORY: raise MemoryError("256MiB pair-bank budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    counts=screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"pair_bank_counts")
    counts.DEADLINE=deadline
    tests=broad.exact_tests(np,screen,peel)
    tests["count_small_checks"]=counts.checks()
    tests["pair_and_correction_checks"]=correction.exact_tests(pair,np,guard)+pair.exact_tests(np,guard)
    tests["FFT_exact_toy_coefficients"]=fft.exact_toy(np,guard)
    tests["scale_checks"]=cells.exact_scale_tests()
    caps=composition.build_caps(args.z,args.alpha,triple,anchor,guard,deadline)
    verifier,core=caps["verifier"],caps["core"]
    tables,single,controls=correction.corrected_pair_logs(caps["powered"],float(args.alpha),
        tuple(range(-16,17,4)),pair,old,np,guard,raw=caps["raw"])
    del caps["raw"],caps["powered"],single
    selected=tables["concavity_Jensen_minimum"]
    del tables
    progress("combined_pair_coefficients_ready",controls=controls,tests=tests,field_tests=caps["tests"],
        centers=centers,ordinary_fallback_enabled=args.ordinary_fallback)
    regions,bank,stats,fallbacks,assignment=marginals(selected,centers,args.modulus,
        args.ordinary_fallback,broad,cells,fft,screen,np,guard,progress)
    control=marginal_control(selected,regions,assignment,fallbacks,args.modulus,pair,
        calibration,broad,cells,fft,screen,np,guard)
    del selected
    progress("complete_pair_marginals_ready",marginals=len(regions),fallbacks=len(fallbacks),control=control)
    h,c_h,eta=anchor.source_factor(WEIGHT); exponent=verifier.exact(args.alpha)
    log_source=float((exponent*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()
                               -core.LIMIT*verifier.exact(args.z).log())).mid())
    report,diagnostics=join_with_diagnostics(regions,fallbacks,counts,broad,screen,np,log_source,
        args.max_queries,args.count_seconds,deadline,guard,progress)
    guard()
    print(json.dumps({"scope":"complete r1024 FLOAT upper-envelope scan with combined pair caps",
        "directed":False,"negative_is_not_a_certificate":True,"whole_code_certificate":False,
        "other_weights_covered":False,"fresh_marker_or_parameter_calibration":False,
        "parameters":{"weight":WEIGHT,"z":str(args.z),"alpha":str(args.alpha),
            "modulus":args.modulus,"centers":centers,"curve_full":str(CURVE_FULL),
            "curve_single":str(CURVE_SINGLE),"rho_exponents":list(range(-16,17,4)),
            "seconds":args.seconds,"count_seconds":args.count_seconds,"max_queries":args.max_queries},
        "joint_missing":{"h":h,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "result":report,"diagnostic_separation":diagnostics,"bank":bank,"FFT_controls":stats,
        "ordinary_fallback":{"explicitly_enabled":args.ordinary_fallback,
            "policy":"enabled by this diagnostic's documented default; use --no-ordinary-fallback to reject",
            "marginals_used":len(fallbacks),"first16":fallbacks[:16],
            "all_fallback_full_indices":[x["full"] for x in fallbacks],
            "method":"ordinary positive Cauchy numerator, never coefficient clipping"},
        "frozen_shared_center_control":control,"controls":controls,"tests":tests,
        "field_tests":caps["tests"],"block_caps":caps["block_report"],
        "composition_pin":COMPOSITION_PIN,"correction_pin":composition.CORRECTION_PIN,
        "triple_pin":composition.TRIPLE_PIN,"pair_pin":correction.PAIR_PIN,
        "broad_pin":calibration.BROAD_PIN,"count_pin":screen.COUNT_PIN,"FFT_pin":cells.FFT_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
