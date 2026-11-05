"""Directed calibrated bands with cell-specific run-count upper masses.

Use actual coordinates for641..4096, or hole coordinates for4097..7551.
The anchor and complete inclusive band must lie on the same side. Build
the coherent local pair table once at the selected z,alpha. A bounded
FLOAT fitter proposes positive rational witnesses; proof numerators use
only the preserved directed coefficients and those exact witnesses.

Each normalized cell pair keeps its source-tilted conditional bound.
Exact run-range masses upper-bound the number of indexed supports in that
cell; these are not exact profile-cell counts. Run intervals may overlap.
The new complete sum charges both ordered terminal chains off diagonal.
The old complete global maximum uses every nonempty geometry cell.
Take the minimum of these two complete upper bounds, then apply the
single global output factor. No interval intersection or extra alpha is
used. beta is already in the256 denominator tilts; eta_r^alpha stays with
each actual size's count. Positive results do not certify the target.

Parent execution only:300s default/600s explicit,512MiB,20GiB free,-B,
stdout only. No work or provider imports run on import. Numerical receipts
are indexed separately; this source alone is not numerical evidence.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time
from types import MappingProxyType

sys.dont_write_bytecode=True
ACTUAL_PIN="ca72174f0fc19848e02de159093f476eb61369629d8dc6f1271332171e7ae635"
COUNTS_PIN="2468da0c91e6e44c7b707a1158252bb695818980de33eb7f270d620748b4f028"
MEMORY=512*(1<<20)


def load_actual():
    path=Path(__file__).with_name("cheap_lane_calibrated_actual_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=ACTUAL_PIN:
        raise RuntimeError("frozen calibrated actual-band mismatch")
    spec=importlib.util.spec_from_file_location("run_band_actual",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def coordinate_mode(first,last,anchor,requested=None):
    if any(type(x) is not int for x in (first,last,anchor)):
        raise TypeError("band and anchor must be exact integer sizes")
    if 641<=first<=last<=4096 and 641<=anchor<=4096:
        mode="actual"
    elif 4097<=first<=last<=7551 and 4097<=anchor<=7551:
        mode="hole"
    else:
        raise ValueError("require a same-side band and anchor in641..4096 or4097..7551")
    if requested is not None and requested!=mode:
        raise ValueError("explicit coordinate convention does not match actual sizes")
    return mode


def wrapper_tests(guard):
    checks=0
    for first,last,anchor,mode in ((641,4096,3072,"actual"),(4096,4096,641,"actual"),
            (4097,7551,4352,"hole"),(4352,4479,4608,"hole")):
        guard(); assert coordinate_mode(first,last,anchor)==mode; checks+=1
        assert coordinate_mode(first,last,anchor,mode)==mode; checks+=1
    for first,last,anchor in ((640,641,641),(4096,4097,4096),(4097,4097,4096),
            (641,642,4097),(7551,7552,7551),(3073,3072,3072)):
        try: coordinate_mode(first,last,anchor)
        except ValueError: checks+=1
        else: raise AssertionError("invalid same-side domain accepted")
    for requested in ("actual","other"):
        try: coordinate_mode(4352,4479,4352,requested)
        except ValueError: checks+=1
        else: raise AssertionError("mismatched coordinate convention accepted")
    try: coordinate_mode(True,4096,3072)
    except TypeError: checks+=1
    else: raise AssertionError("noninteger band accepted")
    return {"same_side_domain_checks":checks,"anchor_may_lie_outside_band":True,
        "crossing4096_forbidden":True}


def cell_bounds(v,band,geometry,den,numerators,marker_logs,cells,klo,khi,mode,join,guard,progress=None):
    """Evaluate each unordered geometry once and retain both ordered logs."""
    result={}; controls=[]; worst=None; worst_upper=None
    digest=hashlib.sha256(); vertices_total=empty=ordered=0
    logtwo=v.arb(2).log()
    for i,left in enumerate(cells):
        for j in range(i,len(cells)):
            guard(); right=cells[j]
            points=geometry.vertices(klo,khi,left[1],left[2],right[1],right[2],guard)
            digest.update((str(i)+","+str(j)+":"+repr(points)+"\n").encode("ascii"))
            if not points: empty+=1; continue
            vertices_total+=len(points); ordered+=1 if i==j else 2
            low,point,bounds=band.evaluate_cell(v,den,points,marker_logs[i],marker_logs[j],
                numerators[i],numerators[j],join,diagonal=i==j)
            result[i,j]=bounds
            expected=("forward",) if i==j else ("forward","reverse")
            if tuple(name for name,_ in bounds)!=expected:
                raise ArithmeticError("incomplete ordered terminal bounds")
            for direction,value in bounds:
                control={"cell_indices":[i,j],"centers":[left[0],right[0]],"order":direction,
                    "log2_conditional_cell_upper_interval":str(value/logtwo)}
                controls.append((value.upper(),control))
                if worst_upper is None or value.upper()>worst_upper:
                    worst_upper=value.upper(); k,a,b=point
                    worst={**control,"coordinate":mode,
                        "minimum_vertex_coordinate_full_full":list(map(str,point)),
                        "actual_size_at_minimum_vertex":str(k if mode=="actual" else 8192-k),
                        "normalized_vertex":[str(2048*a/k),str(2048*b/k)],
                        "first_normalized_cut":list(map(str,left[1:3])),
                        "second_normalized_cut":list(map(str,right[1:3])),
                        "lower_log_denominator_sum":str(v.exact(low))}
        if progress: progress("run_band_geometry",row=i,center=left[0],vertices=vertices_total)
    if not result: raise ArithmeticError("no nonempty normalized cell geometry")
    return result,{"nonempty_unordered_cells":len(result),"empty_unordered_cells":empty,
        "ordered_cells":ordered,"rational_vertices":vertices_total,"vertices_sha256":digest.hexdigest(),
        "worst_conditional_cell":worst,"fractional_vertices_preserved":True,
        "loggamma_cache_entries":len(den.gamma),"log_denominator_evaluations":den.calls,
        "largest_conditional_cells":[row for _,row in sorted(controls,key=lambda x:x[0],reverse=True)[:8]]}


def main():
    actual=load_actual()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor",type=int,default=4352)
    parser.add_argument("--first",type=int,default=4352)
    parser.add_argument("--last",type=int,default=4479)
    parser.add_argument("--coordinate",choices=("actual","hole"))
    parser.add_argument("--z",type=actual.rational64,default=F(7,128))
    parser.add_argument("--alpha",type=actual.rational64,default=F(119,1280))
    parser.add_argument("--seconds",type=int,choices=(180,300,600),default=300)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    args=parser.parse_args()
    mode=coordinate_mode(args.first,args.last,args.anchor,args.coordinate)
    if not 0<args.z<1 or not 0<args.alpha<=1:
        raise ValueError("require0<z<1 and0<alpha<=1")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    support=actual.load_support(); normalized=support.load_normalized(); band=normalized.load_band()
    counts_module=band.load("cheap_lane_cell_run_counts.py",COUNTS_PIN,"run_band_exact_counts")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"run_band_fitting_utilities")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"run_band_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"run_band_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"run_band_ordinary")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"run_band_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"run_band_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"run_band_denominator")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"run_band_terminal")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("calibrated run-band deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("calibrated run-band512MiB budget exceeded")
    def progress(phase,**values):
        selected={key:values[key] for key in ("center","status","evaluations","row","vertices") if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**selected}),flush=True)
    tests={"wrapper":wrapper_tests(guard),"bank":support.wrapper_tests(reference.BANK,guard),
        "fitting_utilities":fitter.exact_tests(guard),"ordinary":ordinary.exact_tests(guard),
        "run_counts":counts_module.exact_tests(join,guard)}
    built=local.build_local_pairs(args.z,args.alpha,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    degrees=built["rational_coefficients"]; arb_identity=dict(built["coefficients"])
    before=support.coefficient_fingerprint(degrees,guard)
    assert before==built["report"]["rational_coefficient_sha256"]
    if mode=="actual":
        proposed,fit_report=actual.fit_actual_markers(MappingProxyType(degrees),args.anchor,fitter,guard,progress)
    else:
        proposed,fit_report=fitter.fit_markers(MappingProxyType(degrees),args.anchor,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    after=support.coefficient_fingerprint(degrees,guard)
    if before!=after or len(arb_identity)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed the shared coherent local table")
    assert all(row is built["coefficients"][key] for key,row in arb_identity.items())
    cells=normalized.normalized_cells(fitted)
    if mode=="actual":
        klo,khi=args.first,args.last
        coefficients=built["coefficients"]
        den=actual.actual_denominators(v,core.BETA,args.alpha,band,guard)
        cover=actual.actual_coverage(normalized,cells,klo,khi,guard)
        tests["actual_coordinates"]=actual.exact_tests(fitted,fitter,normalized,geometry,guard)
        tests["directed_bank"]=actual.directed_tests(v,den,extraction,fitted,args.alpha,core.BETA,join,band,guard)
    else:
        klo,khi=8192-args.last,8192-args.first
        coefficients={(64-f-s,s):row for (f,s),row in built["coefficients"].items()}
        if len(coefficients)!=2145 or any(row is not built["coefficients"][64-f-s,s]
                for (f,s),row in coefficients.items()):
            raise ArithmeticError("complete hole reversal changed local channels")
        den=band.LogDenominators(v,core.BETA,args.alpha,guard)
        cover=normalized.coverage(cells,klo,khi,guard)
        tests["normalized_coordinates"]=normalized.normalized_tests(cells,fitted,geometry,guard)
        tests["directed_bank"]=band.directed_tests(v,den,extraction,fitted,args.alpha,core.BETA,join,guard)
    numerators,marker_logs=normalized.numerator_tables(v,coefficients,cells,ordinary,guard,progress)
    bounds,geometry_report=cell_bounds(v,band,geometry,den,numerators,marker_logs,cells,
        klo,khi,mode,join,guard,progress)
    counts=counts_module.source_counts(v,anchor,join,cells,args.alpha,args.first,args.last,mode,guard)
    old_weight,old_source_report=band.outer_source_sum(v,anchor,args.alpha,args.first,args.last,guard)
    if not counts["full_weight"].overlaps(old_weight):
        raise ArithmeticError("run-count full source sum disagrees with frozen envelope")
    tests["full_source_sum_matches_frozen"]=True
    joined=counts_module.combine(v,counts,bounds,guard)
    old_max=None
    for row in bounds.values():
        for _,value in row:
            old_max=value if old_max is None else old_max.max(value)
    if not joined["old_log_bound"].overlaps(old_weight.log()+old_max):
        raise ArithmeticError("old complete max envelope changed in the count combiner")
    tests["old_complete_max_matches_frozen_formula"]=True
    log_output=-v.exact(args.alpha)*core.LIMIT*v.exact(args.z).log()
    old_log=log_output+joined["old_log_bound"]
    run_log=log_output+joined["run_log_bound"]
    final_log=log_output+joined["minimum_log_bound"]
    guard()
    print(json.dumps({"scope":"complete inclusive actual-size band, all indexed supports and histories at35%",
        "directed_arithmetic":True,"whole_code_certificate":False,"all_integer_sizes_in_band":True,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"fit_anchor_actual":args.anchor,
            "anchor_inside_band":args.first<=args.anchor<=args.last,"coordinate":mode,
            "coordinate_first":klo,"coordinate_last":khi,"z":str(args.z),"alpha":str(args.alpha),
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,
            "memory_budget_bytes":MEMORY},
        "log2_old_global_max_interval":str(old_log/logtwo),
        "log2_cell_run_sum_interval":str(run_log/logtwo),
        "log2_union_bound_interval":str(final_log/logtwo),
        "band_below_2_pow_minus_40":bool(final_log<-40*logtwo),
        "log2_old_eta_weighted_support_count_interval":str(old_weight.log()/logtwo),
        "log2_output_correction_interval":str(log_output/logtwo),
        "run_count_join":joined["report"],"cell_run_counts":counts["report"],
        "old_source_counts_control":old_source_report,"geometry":geometry_report,"coverage":cover,
        "cell_coordinate":"x=2048*coordinate_full/coordinate_weight",
        "cell_count_scope":"exact run-range masses upper-bound cell profile counts; cell ranges may overlap",
        "positive_rational_bank":public_bank,"positive_rational_bank_sha256":bank_hash,
        "fit":{"status":"FLOAT-only proposals independently replayed as positive rational witnesses",
            "not_a_certificate":True,"report":fit_report},
        "local_table_builds":1,"alpha_explicitly_rebuilt":True,"no_coefficient_power_after_fitting":True,
        "local_coefficient_sha256_before_fit":before,"local_coefficient_sha256_after_fit":after,
        "arb_coefficients_preserved_by_identity":True,"directed_degree_reversals":int(mode=="hole"),
        "source_beta_only_in_denominator_tilt":True,"source_eta_powered_once":True,"counts_outside_alpha":True,
        "output_factor_once":True,"both_offdiagonal_terminal_orders_charged":True,
        "completed_upper_bounds_minimized":True,"interval_intersections_used_for_bound":False,
        "boundary_replacement":False,"tests":tests,"local_coefficients":built["report"],
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "nonnegative_intersections":v.intersections,
        "pins":{"actual":ACTUAL_PIN,"counts":COUNTS_PIN,"support":actual.SUPPORT_PIN,
            "normalized":support.NORMALIZED_PIN,"fitter":support.FITTER_PIN,"band":normalized.BAND_PIN,
            "geometry":normalized.GEOMETRY_PIN,"reference":band.HOLE_PIN,"ordinary":reference.ORDINARY_PIN,
            "local":large.LOCAL_PIN,"anchor":large.ANCHOR_PIN,"join":large.JOIN_PIN,"normalizer":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
