"""Complete actual-r6144 FLOAT screen using calibrated hole coordinates.

Build the same actual local caps as cheap_lane_ordered_large_screen:
composition -> residual-rank returns -> ordered129 pair correction.
For actual source weight r>4096, let k=8192-r be the hole weight.
Calibration uses the COMPLETE reversed local table
    H_h[f_h,s]=H[64-f_h-s,s], J_h[f_h,s]=J[64-f_h-s,s].
This is a coefficient relabeling, not evaluation of caps at hole inputs.
The repeated-profile fit targets are256*(f_h,s), with2*f_h+s=k.

Pass the ORIGINAL local logs and actual r to frozen sizes.marginals;
that function performs its own single reversal before residue folding.
Every target uses its own integer hole degrees and denominator
C(4096,f_h)*C(4096-f_h,s)*u^f_h*v^s. The multinomial equals its actual
counterpart at f=r-4096+f_h. Regular and terminal regions remain H^64
and H^63 J; the source factor is Omega_r, never Omega_k.

Fit a fixed31-center grid scaled from the original weight2048 grid.
Only fit endpoints use half-slot proxies; no actual target is changed.
No full-degree-zero boundary is replaced: hole target(0,k) is not an
extremal regional degree when k<4096. All unresolved residues retain
explicit ordinary Cauchy fallback, unless the caller disables it.

Default z1/64,alpha119/1280 preserves the historical r6144 control.
Rational overrides select one witness, not a proof-parameter search.
The optional weight must be even and4098..7550. Odd sizes are rejected
because the frozen complete join checks the even-weight profile count.
Only the requested actual size is covered, not all CLI-eligible sizes.

Numerical output is FLOAT reconnaissance, never a directed certificate.
Import performs no job. Parent-only execution:300s,512MiB,20GiB free,
stdout only, zero generated files. Frozen sources remain unchanged.
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

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
LARGE_PIN="1a0655fde36dbb4f090a143998833bf8ccaf8c25f029df42b1b8837a1db670ed"
MEMORY=512*(1<<20)
LN2=math.log(2)
BASE_CENTERS=(0,3,11,25,44,69,98,131,169,211,256,304,354,406,458,512,
    566,618,670,720,768,813,855,893,926,955,980,999,1013,1021,1024)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def centers_for_holes(weight):
    """Nearest integer scaling, with half ties rounded upward."""
    if type(weight) is not int or weight%2 or not 642<=weight<=4094:
        raise ValueError("even hole weight642..4094 required")
    maximum=weight//2
    centers=tuple((center*maximum+512)//1024 for center in BASE_CENTERS)
    if len(set(centers))!=31 or centers!=tuple(sorted(centers)):
        raise ArithmeticError("scaled center grid must retain31 distinct centers")
    assert centers[0]==0 and centers[-1]==maximum
    return centers


def fit_target(weight,center):
    """Positive endpoint proxies affect only marker proposals."""
    if not 0<=center<=weight//2: raise ValueError("invalid hole center")
    full=F(1,2) if center==0 else F(weight-1,2) if 2*center==weight else F(center)
    single=weight-2*full
    assert full>0 and single>0
    return full,single


def reversed_fit_logs(actual,sizes,np):
    """Reverse every local degree for calibration, preserving all six channels."""
    if actual.shape!=(65,65,6): raise ValueError("expected six pair channels at total degree<=64")
    original={(f,s):actual[f,s] for f in range(65) for s in range(65-f)}
    reversed_table=sizes.reverse_degrees(original,64)
    restored=sizes.reverse_degrees(reversed_table,64)
    result=np.full_like(actual,-np.inf)
    for (f,s),row in reversed_table.items():
        result[f,s]=row
        if not np.array_equal(restored[f,s],original[f,s]):
            raise ArithmeticError("local reversal changed a channel or degree")
    if np.isnan(result).any() or np.isposinf(result).any():
        raise ArithmeticError("invalid reversed local pair log coefficient")
    return result


def exact_hole_tests(large,sizes,broad,guard):
    """Exact complement/count tests and actual-size normalizer controls."""
    tests=sizes.complement_tests(broad,guard)
    checks=0
    assert BASE_CENTERS==large.SMALL_CENTERS
    for actual in (4098,6144,7550):
        guard(); holes=8192-actual; centers=centers_for_holes(holes)
        assert len(centers)==31 and centers[0]==0 and centers[-1]==holes//2
        for center in centers:
            f,s=fit_target(holes,center)
            assert 2*f+s==holes and f>0 and s>0
            checks+=1
        for fh in (0,holes//4,holes//2):
            s=holes-2*fh; fa=actual-4096+fh
            assert 2*fa+s==actual
            assert comb(4096,fa)*comb(4096-fa,s)==comb(4096,fh)*comb(4096-fh,s)
            assert 0<=fh<=4096 and 0<=s<=4096-fh
            checks+=3
        assert holes<4096
        checks+=1
    assert centers_for_holes(2048)==BASE_CENTERS
    for invalid in (641,643,4095,4096):
        try: centers_for_holes(invalid)
        except ValueError: checks+=1
        else: raise AssertionError("unsupported hole grid was accepted")
    # A non-extremal singleton target can combine different local degrees.
    # Consequently H_(0,64)^64 cannot replace the hole f0 coefficient.
    coefficient=F(2)*F(5)+F(3)*F(3)+F(5)*F(2)
    assert coefficient==29 and coefficient!=F(3)**2
    checks+=1
    tests.update({"actual_size_normalizer_and_grid_checks":checks,
        "all_actual_degrees_remain_integer":True,"no_extremal_hole_endpoint_substitution":True,
        "even_weight_only_due_to_frozen_join":True})
    return tests


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight",type=int,default=6144)
    parser.add_argument("--z",type=F,default=F(1,64))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--modulus",type=int,choices=(256,512),default=256)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--ordinary-fallback",action=argparse.BooleanOptionalAction,default=True)
    args=parser.parse_args(); actual=args.weight; holes=8192-actual
    if actual%2 or not 4098<=actual<=7550:
        parser.error("only even actual weights4098..7550: frozen complete_join has an even-only count check")
    if not 0<args.z<1 or not 0<args.alpha<=1: parser.error("require0<z<1,0<alpha<=1")
    if max(abs(n).bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64bits")
    centers=centers_for_holes(holes)
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    large=load("cheap_lane_ordered_large_screen.py",LARGE_PIN,"holes_large_interface")
    sizes=large.load("cheap_lane_pair_sizes_screen.py",large.SIZES_PIN,"holes_sizes")
    rho=large.load("cheap_lane_rho_bucket_control.py",large.RHO_PIN,"holes_coarse")
    ordered=large.load("cheap_lane_rho_ordered_buckets.py",large.ORDERED_PIN,"holes_ordered")
    control=rho.load("cheap_lane_pair_profile_calibration.py",rho.CONTROL_PIN,"holes_control")
    rank=rho.load("cheap_lane_residual_rank_return.py",rho.RANK_PIN,"holes_rank")
    calibration=control.load("cheap_lane_broad_profile_calibration.py",control.CALIBRATION_PIN,"holes_moments")
    bank=sizes.load("cheap_lane_pair_bank_screen.py",sizes.BANK_PIN,"holes_bank")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"holes_broad")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"holes_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"holes_triple")
    pair=correction.load_pair()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"holes_neighbor")
    peel=broad.load("cheap_lane_profile_peeling.py",broad.PEEL_PIN,"holes_peel_tests")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"holes_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("ordered hole-screen deadline; no incomplete cover")
        if peak()>MEMORY: raise MemoryError("512MiB ordered hole-screen budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    tests=exact_hole_tests(large,sizes,broad,guard)
    tests["complete_join_toys"]=broad.exact_tests(np,screen,peel)
    tests["coarse_correction"]=rho.exact_tests(pair,correction,np,guard)
    tests["ordered_correction"]=ordered.exact_tests(ordered.load_fine(),pair,correction,np,guard)
    tests["residual_rank"]=rank.exact_tests(guard)
    tests["positive_moments"]=calibration.exact_tests(np,screen,old)
    tests["FFT_exact_toys"]=fft.exact_toy(np,guard)
    tests["scale_checks"]=cells.exact_scale_tests()
    caps=rank.augment_caps(composition.build_caps(args.z,args.alpha,triple,anchor,guard,deadline),
        args.z,args.alpha,guard)
    verifier,core=caps["verifier"],caps["core"]
    tables,local_report=ordered.ordered_bucket_pair_logs(caps["powered"],float(args.alpha),
        pair,old,correction,np,guard,caps["raw"],progress,coarse_builder=rho.bucket_pair_logs)
    actual_logs=tables["ordered129_bucket_sum"]
    del tables,caps["raw"],caps["powered"]
    fit_logs=reversed_fit_logs(actual_logs,sizes,np)
    objective=control.Objective(fit_logs,0.0,calibration,np,guard)
    witnesses={}; fit_reports=[]
    for center in centers:
        guard(); fstar,sstar=fit_target(holes,center)
        objective.targets=256*np.asarray([float(fstar),float(sstar)])
        start=sizes.markers(holes,fstar)
        best,history,status=control.coordinate_fit(objective,np.asarray(list(map(broad.logq,start))),np,guard)
        marker=tuple(max(F(1,1<<30),F(math.exp(float(x))).limit_denominator(1<<20)) for x in best["x"])
        rounded=objective.evaluate(np.asarray(list(map(broad.logq,marker))))
        witnesses[center]=marker
        report={"hole_center":center,"actual_full_at_center":actual-4096+center,
            "fit_hole_targets":[str(fstar),str(sstar)],"markers":list(map(str,marker)),
            "evaluations":len(history),"status":status,
            "ordinary_objective_gain_bits_FLOAT":history[0]["log2_ordinary_FLOAT"]-rounded["objective"]/LN2,
            "rounded_region_hole_means_FLOAT":[float(x/256) for x in rounded["means"]],
            "rounded_per_region_log_marker_gradient_bits_FLOAT":[float(x/(256*LN2)) for x in rounded["gradient"]]}
        fit_reports.append(report); progress("ordered_hole_center_fitted",**report)
    del objective,fit_logs
    original_markers=sizes.markers
    def calibrated_markers(coordinate_weight,center):
        assert coordinate_weight==holes and center in witnesses
        return witnesses[center]
    sizes.markers=calibrated_markers
    try:
        # This call intentionally receives ORIGINAL actual-degree logs.
        # Frozen marginals reverses them exactly once using actual>4096.
        regions,bank_report,stats,fallbacks=sizes.marginals(actual_logs,actual,centers,args.modulus,
            args.ordinary_fallback,broad,cells,fft,screen,np,guard,progress)
    finally:
        sizes.markers=original_markers
    del actual_logs
    assert len(regions)==holes//2+1
    for f,region in enumerate(regions):
        center=min(centers,key=lambda c:(abs(c-f),c))
        assert region[3:]==witnesses[center]
    missing,c_h,eta=anchor.source_factor(actual)
    assert missing==actual//256 and c_h==1-F(actual,256*(missing+1))
    if actual==6144: assert missing==24 and c_h==F(1,25)
    log_source=float((verifier.exact(args.alpha)*(verifier.exact(eta).log()
        +actual*verifier.exact(core.BETA).log()-core.LIMIT*verifier.exact(args.z).log())).mid())
    progress("ordered_hole_marginals_ready",actual_weight=actual,hole_weight=holes,
        marginals=len(regions),ordinary_fallbacks=len(fallbacks),boundary_replaced=False)
    result=sizes.complete_join(actual,regions,fallbacks,broad,screen,np,log_source,guard,progress)
    assert result["ordered_profiles"]==(holes//2+1)**2-1
    assert result["actual_weight"]==actual and result["coordinate_weight"]==holes
    guard()
    print(json.dumps({"scope":"complete requested actual source size; FLOAT current-operator hole-coordinate run-mass envelope",
        "directed":False,"negative_is_not_a_certificate":True,"whole_code_certificate":False,
        "other_weights_covered":False,"odd_weights_supported":False,
        "odd_weight_restriction":"frozen complete_join verifies even-weight ordered-profile cardinality",
        "proof_parameter_search":False,"marker_fit_only":True,
        "parameters":{"actual_weight":actual,"hole_weight":holes,"hole_coordinates":True,
            "z":str(args.z),"alpha":str(args.alpha),"paired_modulus":args.modulus,
            "hole_centers":centers,"max_fit_evaluations":control.MAX_PROPOSALS,
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "pair_bound":"current residual-rank returns and ordered129; N4/coarse/Jensen ceilings retained",
        "proof_interface":"CHEAP_LANE_RANK_INTERFACE.md: H_h(u,v)=u^64 H(1/u,v/u); exact correction between different ratio buckets",
        "reversal":{"complete_local_table":True,"caps_evaluated_at_actual_inputs":True,
            "calibration_reversals":1,"regional_reversals":1,"regional_input_logs":"original actual-degree table",
            "actual_full_from_holes":"r-4096+f_h","local_capacity":64,"regional_capacity":4096,
            "same_MM_MT_order":True},
        "normalization":"C(4096,f_h) C(4096-f_h,k-2f_h) u^f_h v^(k-2f_h); own integer target",
        "fits":fit_reports,"bank":bank_report,"result":result,"boundary_replaced":False,
        "FFT_controls":stats,"ordinary_fallback":{"explicitly_enabled":args.ordinary_fallback,
            "marginals_used":len(fallbacks),"first16":fallbacks[:16],
            "all_hole_full_indices":[row["coordinate_full"] for row in fallbacks],
            "all_actual_full_indices":[row["actual_full"] for row in fallbacks],
            "method":"ordinary positive Cauchy numerator; no clipping or omitted profiles"},
        "fit_scope":"ordinary repeated-hole-profile proposals; no projected stationarity claim",
        "endpoint_proxy_changes_actual_degrees":False,"source_factor_applied_once":True,
        "source_factor_weight":actual,"source_factor_uses_hole_weight":False,
        "joint_missing_for_actual_weight":{"h":missing,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "local_ordered_report":local_report,"tests":tests,"field_tests":caps["tests"],
        "residual_rank_receipt":caps["residual_rank_receipt"],
        "pins":{"large":LARGE_PIN,"sizes":large.SIZES_PIN,"rho":large.RHO_PIN,
            "ordered":large.ORDERED_PIN,"rank":rho.RANK_PIN,"control":rho.CONTROL_PIN,
            "moments":control.CALIBRATION_PIN,"composition":bank.COMPOSITION_PIN,
            "correction":composition.CORRECTION_PIN,"triple":composition.TRIPLE_PIN,
            "pair":correction.PAIR_PIN,"FFT":cells.FFT_PIN,"core":anchor.CORE_PIN},
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
