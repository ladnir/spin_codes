"""Directed normalized ACTUAL-size bands with proposal-only marker fitting.

For641<=r<=4096, use actual full/single counts (a,r-2a), not holes.
The normalized coordinate is x=2048*a/r. The frozen nine-plane geometry
is valid unchanged: empty counts4096-r+a are nonnegative, and a+b<=r-1
excludes the zero-run profile. Actual local H/J coefficients are never
reversed, either for fitting or for the directed ordinary numerators.

Build one coherent directed local table at the selected exact z,alpha.
A bounded FLOAT procedure proposes31 positive rational marker pairs.
Replay those witnesses against the preserved Arb coefficients. The
directed log denominator contains -alpha*r*log(beta)/256; the256 region
factors restore beta^(alpha*r). Each actual integer size, including odd
sizes, contributes eta_r^alpha times its indexed support count. Only the
count stays outside alpha.

No FLOAT matrix, fitted objective, endpoint proxy, new coefficient power,
marker union factor, or full-field distribution enters the proof. The
anchor may lie outside the requested band: it selects deterministic proof
witnesses only. A completed positive upper bound need not be useful.

Parent execution only:300s/512MiB/20GiB free,-B,stdout,no generated files.
Numerical receipts are indexed separately; source alone is not evidence.
No provider import or numerical work runs on module import.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from math import factorial
from pathlib import Path
import shutil
import sys
import time
from types import MappingProxyType

sys.dont_write_bytecode=True
SUPPORT_PIN="95397e9c6a3fd51272d79b36d911fa0ef144927a479f5be44a7a3226dd2d9354"
MEMORY=512*(1<<20)
REFERENCE_WEIGHT=2048


def load_support():
    path=Path(__file__).with_name("cheap_lane_calibrated_hole_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SUPPORT_PIN:
        raise RuntimeError("frozen calibrated-band support mismatch")
    spec=importlib.util.spec_from_file_location("actual_calibrated_support",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def rational64(text):
    if not isinstance(text,str) or len(text)>128 or "e" in text.lower():
        raise argparse.ArgumentTypeError("use a bounded fraction or decimal without exponent notation")
    try: value=F(text)
    except (ValueError,ZeroDivisionError) as exc:
        raise argparse.ArgumentTypeError("expected an exact rational") from exc
    if max(abs(value.numerator).bit_length(),value.denominator.bit_length())>64:
        raise argparse.ArgumentTypeError("rational witnesses limited to64 bits")
    return value


def fit_target(actual,center):
    """Positive interior proxies apply only to the FLOAT proposal objective."""
    if type(actual) is not int or not 641<=actual<=4096:
        raise ValueError("actual fitting anchor must be641..4096")
    if type(center) is not int or not 0<=center<=1024:
        raise ValueError("normalized center must be0..1024")
    full=F(actual*center,REFERENCE_WEIGHT)
    if center==0: full=F(1,2)
    elif center==1024: full=F(actual-1,2)
    single=actual-2*full
    if min(full,single,F(4096-actual)+full)<=0:
        raise ArithmeticError("fitting proxy must have positive full/single/empty counts")
    return full,single


def fit_actual_markers(rational_degrees,actual,fitter,guard,progress=None):
    """Use pinned low-level utilities, not the hole fitter or a fake anchor."""
    fit_target(actual,0)
    template,control,moments,sizes,pins=fitter._dependencies(guard)
    import numpy as np
    table,fingerprint,positives,zero_count=fitter._prepare(rational_degrees,guard)
    logs=np.full((65,65,6),-np.inf,dtype=float)
    for (f,s),row in table.items():
        guard(); logs[f,s]=tuple(fitter._log_fraction(value) for value in row)
        if row!=rational_degrees[f,s]:
            raise ArithmeticError("actual degree labels changed in proposal conversion")
    if np.isnan(logs).any() or np.isposinf(logs).any():
        raise ArithmeticError("invalid actual-degree proposal log table")
    del table
    objective=control.Objective(logs,0.0,moments,np,guard)
    bank=[]; reports=[]; selected=[]; failures=[]; evaluations=0
    for center,lo,hi,_,_ in template:
        guard(); full,single=fit_target(actual,center)
        objective.targets=256*np.asarray([float(full),float(single)])
        seed=tuple(fitter._positive_rational(value) for value in sizes.markers(actual,full))
        marker=seed; history=[]; rounded=None
        status="positive rational seed retained"
        try:
            best,history,status=control.coordinate_fit(objective,
                np.asarray([fitter._log_fraction(value) for value in seed]),np,guard)
            marker=tuple(fitter._positive_rational(F.from_float(math.exp(float(x)))) for x in best["x"])
            rounded=objective.evaluate(np.asarray([fitter._log_fraction(value) for value in marker]))
        except (ArithmeticError,ValueError) as error:
            marker=seed; rounded=None
            status="FLOAT fit failed; positive rational seed retained"
            failures.append({"center":center,"reason":type(error).__name__+": "+str(error)})
        if len(history)>fitter.MAX_PROPOSALS:
            raise ArithmeticError("frozen coordinate fit exceeded12 proposals")
        if not all(value>0 and value.denominator<=fitter.DENOMINATOR_LIMIT for value in marker):
            raise ArithmeticError("invalid exported rational witness")
        bank.append((center,lo,hi,*marker)); evaluations+=len(history)
        report={"center":center,"status":status,"evaluations":len(history)}
        if rounded is not None and history:
            report["ordinary_fit_gain_bits_FLOAT"]=history[0]["log2_ordinary_FLOAT"]-rounded["objective"]/math.log(2)
        reports.append(report)
        if center in (0,256,512,768,1024):
            detail={"center":center,"fit_targets":[str(full),str(single)]}
            if rounded is not None:
                detail["rounded_region_means_FLOAT"]=[float(x/256) for x in rounded["means"]]
                detail["rounded_region_gradient_FLOAT"]=[float(x/256) for x in rounded["gradient"]]
            selected.append(detail)
        if progress: progress("actual_marker_fit",center=center,status=status,evaluations=len(history))
    if tuple(row[:3] for row in bank)!=tuple(row[:3] for row in template):
        raise ArithmeticError("fitting changed reference coverage")
    return tuple(bank),{
        "scope":"FLOAT-only proposals from unchanged actual-degree rational coefficients",
        "certificate":False,"actual_anchor":actual,"coordinate_weight":actual,
        "reference_weight":REFERENCE_WEIGHT,"reference_centers":31,"degree_reversals":0,
        "local_cap_builds":0,"local_degree_rows":2145,"input_utility_sha256":fingerprint,
        "positive_coefficients_per_channel":positives,"zero_coefficients":zero_count,
        "frozen_coordinate_evaluation_limit":fitter.MAX_PROPOSALS,
        "successful_proposal_evaluations_recorded":evaluations,
        "one_rounded_witness_check_after_each_successful_fit":True,
        "marker_denominator_limit":fitter.DENOMINATOR_LIMIT,"positive_floor":"1/1048576",
        "objective_numerator":"e0 H^16383 J","objective_targets":"256*(actual_full,actual_single)",
        "marker_independent_objective_constant":0,"source_or_count_factor_applied":False,
        "endpoint_proxies_only_for_fitting":True,"fits":reports,
        "selected_fit_diagnostics":selected,"fit_failures":failures,"pins":pins}


def actual_denominators(v,beta,alpha,band,guard):
    """Reuse cached directed Gamma, explicitly replacing the hole source tilt."""
    class ActualLogDenominators(band.LogDenominators):
        def value(self,r,a,logu,logw):
            r,a=F(r),F(a); single=r-2*a; empty=4096-r+a
            if not 1<=r<=4096 or min(a,single,empty)<0:
                raise ValueError("invalid exact actual regional counts")
            self.calls+=1
            return (self.logfactorial-self.lgamma(a+1)-self.lgamma(single+1)
                -self.lgamma(empty+1)+self.v.exact(a)*logu+self.v.exact(single)*logw
                -self.alpha*self.v.exact(r/256)*self.logbeta)
    return ActualLogDenominators(v,beta,alpha,guard)


def actual_coverage(normalized,cells,first,last,guard):
    # The frozen routine is coordinate-neutral; only its label says "hole".
    report=dict(normalized.coverage(cells,first,last,guard))
    report["actual_sizes_checked"]=report.pop("hole_sizes_checked")
    report["coordinate_weight_is_actual_size"]=True
    return report


def exact_tests(template,fitter,normalized,geometry,guard):
    checks=0
    for r in (641,642,2047,2048,3072,4095,4096):
        for center,_,_,_,_ in template:
            guard(); full,single=fit_target(r,center)
            assert 2*full+single==r and min(full,single,4096-r+full)>0; checks+=1
            if center not in (0,1024):
                assert full==F(r*center,2048); checks+=1
            else:
                assert full==(F(1,2) if center==0 else F(r-1,2)); checks+=1
    for bad in (640,4097,True,F(3072)):
        try: fit_target(bad,512)
        except ValueError: checks+=1
        else: raise AssertionError("invalid actual fitting anchor accepted")
    cells=normalized.normalized_cells(template)
    controls=[]
    for first,last in ((641,643),(2047,2049),(3071,3073),(4095,4096)):
        controls.append(actual_coverage(normalized,cells,first,last,guard)); checks+=1
    for r,a,b in ((5,1,2),(641,0,320),(2047,511,1023),(4096,0,2048)):
        x,y=F(2048*a,r),F(2048*b,r)
        assert geometry.vertices(r,r,x,x,y,y,guard)==((F(r),F(a),F(b)),); checks+=1
    assert geometry.vertices(4096,4096,1024,1024,1024,1024,guard)==(); checks+=1
    toy={(f,s):tuple(F(1+f+3*s+c,17+c) for c in range(6))
         for f in range(3) for s in range(3-f)}
    original=dict(toy); copied,_,_,_=fitter._prepare(MappingProxyType(toy),guard,capacity=2)
    assert copied==original and toy==original; checks+=1
    u,w=F(2,3),F(5,7)
    for c in range(6):
        assert sum(row[c]*u**f*w**s for (f,s),row in copied.items())==sum(
            toy[f,s][c]*u**f*w**s for f in range(3) for s in range(3-f)); checks+=1
    return {"actual_coordinate_exact_checks":checks,"anchor4096_supported":True,
        "actual_degree_copy_not_reversed":True,"zero_run_excluded":True,
        "odd_size_coverage_controls":controls,"geometry":geometry.exact_tests(guard)}


def directed_tests(v,den,extraction,bank,alpha,beta,join,band,guard):
    checks=0
    for n in (0,1,2,3,7,16,32,4096):
        guard(); assert den.lgamma(F(n+1)).overlaps(v.exact(F(factorial(n))).log()); checks+=1
    marker=(bank[15][3],bank[15][4]); logs=tuple(v.exact(x).log() for x in marker)
    for r,a in ((5,1),(6,3),(641,0),(2047,511),(3072,768),(4096,0),(4096,2048)):
        exact=extraction.target_denominator(r,a,marker)
        direct=v.exact(exact).log()-v.exact(alpha)*v.exact(F(r,256))*v.exact(beta).log()
        evaluated=den.value(r,a,*logs)
        assert evaluated.overlaps(direct); checks+=1
        low=v.exact(band.lower_fraction(evaluated))
        assert low<=direct or low.overlaps(direct); checks+=1
    markers=(marker,(bank[16][3],bank[16][4]))
    regions=tuple((tuple(v.exact(F(n,17)) for n in (2+t,3,5,7+t)),
                   tuple(v.exact(F(n,19)) for n in (11+t,13))) for t in (0,1))
    marker_logs=tuple(tuple(v.exact(x).log() for x in m) for m in markers)
    point=(F(5),F(1),F(2))
    low,vertex,bounds=band.evaluate_cell(v,den,(point,),marker_logs[0],marker_logs[1],
        regions[0],regions[1],join)
    assert vertex==point; checks+=1
    divided=[]
    for region,a,m in zip(regions,(1,2),markers):
        divisor=v.exact(extraction.target_denominator(5,a,m))
        divided.append(tuple(tuple(x/divisor for x in row) for row in region))
    exact_log_den=den.value(5,1,*marker_logs[0])+den.value(5,2,*marker_logs[1])
    for (_,bound),chain in zip(bounds,join.two_orientations(*divided)):
        direct=chain.log()+v.exact(alpha)*5*v.exact(beta).log()
        assert bound>=direct or bound.overlaps(direct); checks+=1
        assert (bound-128*(exact_log_den-v.exact(low))).overlaps(direct); checks+=1
    return {"actual_directed_checks":checks,"source_exponent_is5_not8187":True,
        "singleton_cell_both_terminal_orders":True,"source_tilt_has256_factors":True,
        "fractional_vertices_not_rounded":True,"lower_log_endpoint_direction_checked":True}


def cell_envelope(v,band,geometry,den,numerators,marker_logs,cells,first,last,join,guard,progress=None):
    """Same frozen cell evaluator, with actual coordinates and truthful labels."""
    maximum=None; worst=None; worst_upper=None; controls=[]
    vertex_hash=hashlib.sha256(); vertex_count=nonempty=empty=ordered=0
    logtwo=v.arb(2).log()
    for i,left in enumerate(cells):
        for j in range(i,len(cells)):
            guard(); right=cells[j]
            vertices=geometry.vertices(first,last,left[1],left[2],right[1],right[2],guard)
            vertex_hash.update((str(i)+","+str(j)+":"+repr(vertices)+"\n").encode("ascii"))
            if not vertices: empty+=1; continue
            nonempty+=1; vertex_count+=len(vertices); ordered+=1 if i==j else 2
            low,point,bounds=band.evaluate_cell(v,den,vertices,marker_logs[i],marker_logs[j],
                numerators[i],numerators[j],join,diagonal=i==j)
            pair_upper=None; winner=None
            for direction,logbound in bounds:
                maximum=logbound if maximum is None else maximum.max(logbound)
                if pair_upper is None or logbound.upper()>pair_upper.upper():
                    pair_upper,winner=logbound,direction
                if worst_upper is None or logbound.upper()>worst_upper:
                    worst_upper=logbound.upper(); r,a,b=point
                    worst={"centers":[left[0],right[0]],"order":direction,
                        "normalized_first_cut":list(map(str,left[1:3])),
                        "normalized_second_cut":list(map(str,right[1:3])),
                        "minimum_vertex_actual_full_full":list(map(str,point)),
                        "actual_size_at_minimum_vertex":str(r),
                        "normalized_vertex":[str(2048*a/r),str(2048*b/r)],
                        "lower_log_denominator_sum":str(v.exact(low)),
                        "log2_conditional_cell_upper_interval":str(logbound/logtwo)}
            controls.append((pair_upper.upper(),{"centers":[left[0],right[0]],"vertices":len(vertices),
                "larger_order":winner,"log2_cell_upper_interval":str(pair_upper/logtwo)}))
        if progress: progress("actual_band_geometry",row=i,center=left[0],vertices=vertex_count)
    if maximum is None: raise ArithmeticError("no nonempty actual normalized cells")
    return maximum,worst,{"nonempty_unordered_cells":nonempty,"empty_unordered_cells":empty,
        "ordered_cells":ordered,"rational_vertices":vertex_count,"vertices_sha256":vertex_hash.hexdigest(),
        "loggamma_cache_entries":len(den.gamma),"log_denominator_evaluations":den.calls,
        "fractional_vertices_preserved":True,"coordinate_weight_is_actual_size":True,
        "largest_cell_controls":[row for _,row in sorted(controls,key=lambda x:x[0],reverse=True)[:8]]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor",type=int,default=3072)
    parser.add_argument("--first",type=int,default=2944)
    parser.add_argument("--last",type=int,default=3200)
    parser.add_argument("--z",type=rational64,default=F(1,8))
    parser.add_argument("--alpha",type=rational64,default=F(119,1280))
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    args=parser.parse_args()
    if not 641<=args.first<=args.last<=4096 or not 641<=args.anchor<=4096:
        raise ValueError("actual anchor and inclusive band must lie in641..4096")
    if not 0<args.z<1 or not 0<args.alpha<=1:
        raise ValueError("require0<z<1 and0<alpha<=1")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    support=load_support(); normalized=support.load_normalized(); band=normalized.load_band()
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"actual_fitting_utilities")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"actual_normalized_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"actual_reference_cuts")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"actual_ordinary")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"actual_source_interface")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"actual_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"actual_denominator")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"actual_terminal")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("actual band deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("actual band512MiB budget exceeded")
    def progress(phase,**values):
        selected={key:values[key] for key in ("center","status","evaluations","row","vertices") if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**selected}),flush=True)
    tests={"wrapper":support.wrapper_tests(reference.BANK,guard),
        "fitting_utilities":fitter.exact_tests(guard),
        "actual":exact_tests(reference.BANK,fitter,normalized,geometry,guard),
        "ordinary":ordinary.exact_tests(guard)}
    built=local.build_local_pairs(args.z,args.alpha,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    degrees=built["rational_coefficients"]
    before=support.coefficient_fingerprint(degrees,guard)
    assert before==built["report"]["rational_coefficient_sha256"]
    arb_identity=dict(built["coefficients"])
    proposed,fit_report=fit_actual_markers(MappingProxyType(degrees),args.anchor,fitter,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    after=support.coefficient_fingerprint(degrees,guard)
    if before!=after or len(arb_identity)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed the coherent local table")
    assert all(row is built["coefficients"][key] for key,row in arb_identity.items())
    cells=normalized.normalized_cells(fitted)
    cover=actual_coverage(normalized,cells,args.first,args.last,guard)
    den=actual_denominators(v,core.BETA,args.alpha,band,guard)
    tests["directed_actual_bank"]=directed_tests(v,den,extraction,fitted,args.alpha,core.BETA,join,band,guard)
    numerators,marker_logs=normalized.numerator_tables(v,built["coefficients"],cells,ordinary,guard,progress)
    maximum,worst,geometry_report=cell_envelope(v,band,geometry,den,numerators,marker_logs,cells,
        args.first,args.last,join,guard,progress)
    weighted,source_report=band.outer_source_sum(v,anchor,args.alpha,args.first,args.last,guard)
    if not weighted>0: raise ArithmeticError("nonpositive weighted source sum")
    log_count=weighted.log(); log_output=-v.exact(args.alpha)*core.LIMIT*v.exact(args.z).log()
    log_total=log_count+log_output+maximum
    guard()
    print(json.dumps({"scope":"all indexed supports in the inclusive actual-size band and all histories at35%",
        "directed_arithmetic":True,"whole_code_certificate":False,"all_integer_sizes_in_band":True,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"fit_anchor_actual":args.anchor,
            "anchor_inside_band":args.first<=args.anchor<=args.last,"z":str(args.z),"alpha":str(args.alpha),
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,
            "memory_budget_bytes":MEMORY},
        "log2_union_bound_interval":str(log_total/logtwo),
        "band_below_2_pow_minus_40":bool(log_total<-40*logtwo),
        "band_below_2_pow_minus_60":bool(log_total<-60*logtwo),
        "log2_eta_weighted_support_count_interval":str(log_count/logtwo),
        "log2_output_correction_interval":str(log_output/logtwo),
        "worst_cell":worst,"geometry":geometry_report,"coverage":cover,"source_counts":source_report,
        "cell_coordinate":"x=2048*actual_full/actual_weight","continuous_cut_rule":"reference hi+1/2; endpoints0,1024",
        "positive_rational_bank":public_bank,"positive_rational_bank_sha256":bank_hash,
        "fit":{"status":"FLOAT proposals only; rational witnesses independently replayed with directed arithmetic",
            "not_a_certificate":True,"report":fit_report},
        "local_table_builds":1,"alpha_explicitly_rebuilt":True,"no_coefficient_power_after_fitting":True,
        "local_coefficient_sha256_before_fit":before,"local_coefficient_sha256_after_fit":after,
        "arb_coefficients_preserved_by_identity":True,"actual_degree_fit_input":True,
        "fit_degree_reversals":0,"directed_degree_reversals":0,"source_beta_only_in_denominator_tilt":True,
        "output_factor_once":True,"boundary_replacement":False,"extra_profile_or_marker_factor":False,
        "tests":tests,"local_coefficients":built["report"],"elapsed_seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"generated_bytes":0,"nonnegative_intersections":v.intersections,
        "pins":{"support":SUPPORT_PIN,"normalized":support.NORMALIZED_PIN,"fitting_utilities":support.FITTER_PIN,
            "band":normalized.BAND_PIN,"geometry":normalized.GEOMETRY_PIN,"reference_cuts":band.HOLE_PIN,
            "ordinary":reference.ORDINARY_PIN,"local":large.LOCAL_PIN,"anchor":large.ANCHOR_PIN,
            "join":large.JOIN_PIN,"normalizer":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
