"""Complete r2048 FLOAT rho-bucket screen with one fixed calibrated bank.

Keep z23/100,alpha13/128, residual-rank returns and the frozen dense33
disjoint-left bucket pair cap. Build the caps ONCE. A fixed31-center
variance-shaped grid is calibrated with the frozen twelve-evaluation
coordinate fitter; no grid or proof-parameter search is performed.

The fit is a repeated-profile finite-chain witness objective. Endpoint
centers0/1024 use strictly interior fit targets f*=1/2,2047/2 with
s*=2048-2f*. These are proposal targets ONLY. Every actual extracted
profile has integer f,s, its own exact multinomial denominator and the
usual positive marker powers. No interpolated matrix or probability is
used. All rounded rational witnesses are positive and replayed directly.

Every one of1025 marginals and every compatible ordered profile is
retained, including explicitly reported ordinary Cauchy fallbacks.
The frozen complete join charges exact indexed run masses outside alpha
and applies Omega2048 once. H64/H63J preserves32767 updates and terminal.
This is FLOAT reconnaissance, not a certificate or whole-code result.
One180/300s job,512MiB,20GiB free, stdout only; no generated files.
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
SIZES_PIN="7a7024b78912abb78eaaef9e90f766eefdd7e7f807d8dfe6143f38c86c233afd"
RHO_PIN="57391fb3ab3ada8078c17a98eb921dbb5777452422c429960e3fa47aa07624b8"
WEIGHT=2048
Z,ALPHA=F(23,100),F(13,128)
CENTERS=(0,3,11,25,44,69,98,131,169,211,256,304,354,406,458,512,
         566,618,670,720,768,813,855,893,926,955,980,999,1013,1021,1024)
MEMORY=512*(1<<20)
LN2=math.log(2)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def fit_target(center):
    f=F(1,2) if center==0 else F(2047,2) if center==1024 else F(center)
    return f,WEIGHT-2*f


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modulus",type=int,choices=(256,512),default=512)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--ordinary-fallback",action=argparse.BooleanOptionalAction,default=True)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    sizes=load("cheap_lane_pair_sizes_screen.py",SIZES_PIN,"rho_fitbank_sizes")
    rho=load("cheap_lane_rho_bucket_control.py",RHO_PIN,"rho_fitbank_buckets")
    control=rho.load("cheap_lane_pair_profile_calibration.py",rho.CONTROL_PIN,"rho_fitbank_control")
    rank=rho.load("cheap_lane_residual_rank_return.py",rho.RANK_PIN,"rho_fitbank_rank")
    calibration=control.load("cheap_lane_broad_profile_calibration.py",control.CALIBRATION_PIN,"rho_fitbank_moments")
    bank=sizes.load("cheap_lane_pair_bank_screen.py",sizes.BANK_PIN,"rho_fitbank_bank")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"rho_fitbank_broad")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"rho_fitbank_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"rho_fitbank_triple")
    pair=correction.load_pair()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"rho_fitbank_neighbor")
    peel=broad.load("cheap_lane_profile_peeling.py",broad.PEEL_PIN,"rho_fitbank_peel")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"rho_fitbank_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("calibrated bank deadline; no automatic retry")
        if peak()>MEMORY: raise MemoryError("512MiB calibrated bank budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    assert len(CENTERS)==31 and tuple(sorted(set(CENTERS)))==CENTERS
    assert fit_target(0)==(F(1,2),F(2047)) and fit_target(1024)==(F(2047,2),F(1))
    tests=broad.exact_tests(np,screen,peel)
    tests["rho"]=rho.exact_tests(pair,correction,np,guard)
    tests["residual_rank"]=rank.exact_tests(guard)
    tests["positive_moments"]=calibration.exact_tests(np,screen,old)
    tests["FFT_exact_toys"]=fft.exact_toy(np,guard)
    tests["scale_checks"]=cells.exact_scale_tests()
    tests["fixed_grid_and_positive_endpoint_fit_targets"]=True
    caps=rank.augment_caps(composition.build_caps(Z,ALPHA,triple,anchor,guard,deadline),Z,ALPHA,guard)
    verifier,core=caps["verifier"],caps["core"]
    tables,local_report=rho.bucket_pair_logs(caps["powered"],float(ALPHA),pair,old,
        correction,np,guard,caps["raw"],progress)
    selected=tables["dense33_bucket_sum"]
    del tables,caps["raw"],caps["powered"]
    # Fixed scalar normalizers/source counts do not affect a marker fit.
    # Reuse the six-moment geometry, replacing ONLY the requested fit targets.
    objective=control.Objective(selected,0.0,calibration,np,guard)
    witnesses={}; fit_reports=[]
    for center in CENTERS:
        guard(); fstar,sstar=fit_target(center)
        objective.targets=256*np.asarray([float(fstar),float(sstar)])
        start=sizes.markers(WEIGHT,fstar)
        best,history,status=control.coordinate_fit(objective,np.asarray(list(map(broad.logq,start))),np,guard)
        markers=tuple(max(F(1,1<<30),F(math.exp(float(x))).limit_denominator(1<<20)) for x in best["x"])
        rounded=objective.evaluate(np.asarray(list(map(broad.logq,markers))))
        witnesses[center]=markers
        row={"center":center,"fit_targets":[str(fstar),str(sstar)],"markers":list(map(str,markers)),
            "evaluations":len(history),"status":status,
            "ordinary_objective_gain_bits_FLOAT":history[0]["log2_ordinary_FLOAT"]-rounded["objective"]/LN2,
            "rounded_region_means_FLOAT":[float(x/256) for x in rounded["means"]],
            "rounded_per_region_log_marker_gradient_bits_FLOAT":[float(x/(256*LN2)) for x in rounded["gradient"]]}
        fit_reports.append(row); progress("rho_bank_center_fitted",**row)
    del objective
    original_markers=sizes.markers
    def calibrated_markers(weight,center):
        assert weight==WEIGHT and center in witnesses
        return witnesses[center]
    sizes.markers=calibrated_markers
    try:
        regions,bank_report,stats,fallbacks=sizes.marginals(selected,WEIGHT,CENTERS,args.modulus,
            args.ordinary_fallback,broad,cells,fft,screen,np,guard,progress)
    finally:
        sizes.markers=original_markers
    del selected
    assert len(regions)==1025
    for f,region in enumerate(regions):
        center=min(CENTERS,key=lambda c:(abs(c-f),c))
        assert region[3:]==witnesses[center]
    missing,c_h,eta=anchor.source_factor(WEIGHT)
    log_source=float((verifier.exact(ALPHA)*(verifier.exact(eta).log()
        +WEIGHT*verifier.exact(core.BETA).log()-core.LIMIT*verifier.exact(Z).log())).mid())
    progress("rho_calibrated_marginals_ready",marginals=len(regions),ordinary_fallbacks=len(fallbacks))
    anchor_controls=[]; fallback_full={row["coordinate_full"] for row in fallbacks}
    for f in (768,898):
        c=WEIGHT-2*f
        numerator=8192*comb(WEIGHT-1,c-1)*comb(8192-WEIGHT-1,c-1)
        count,remainder=divmod(numerator,c); assert remainder==0
        forward,reverse=screen.paired_chains(regions[f],regions[f]); assert abs(forward-reverse)<1e-7
        anchor_controls.append({"actual_full":[f,f],"single":[c,c],"runs":c,
            "markers":list(map(str,regions[f][3:])),"uses_ordinary_fallback":f in fallback_full,
            "log2_selected_full_run_mass_charged_FLOAT":(log_source+math.log(count)+forward)/LN2})
    result=sizes.complete_join(WEIGHT,regions,fallbacks,broad,screen,np,log_source,guard,progress)
    guard()
    print(json.dumps({"scope":"complete actual r2048 FLOAT rho-bucket calibrated-bank run-mass envelope",
        "directed":False,"negative_is_not_a_certificate":True,"whole_code_certificate":False,
        "other_weights_covered":False,"proof_parameter_search":False,"marker_fit_only":True,
        "parameters":{"actual_weight":WEIGHT,"z":str(Z),"alpha":str(ALPHA),
            "paired_modulus":args.modulus,"centers":CENTERS,"max_proposals_per_center":control.MAX_PROPOSALS,
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "pair_bound":"same dense33 bucket table at every center; old/dense/Jensen ceilings retained",
        "fit_scope":"repeated-profile ordinary finite-chain witnesses; not stationary residue or per-entry optima",
        "endpoint_proxies_affect_extraction":False,"fits":fit_reports,"bank":bank_report,
        "result":result,"selected_anchor_controls":anchor_controls,
        "selected_anchor_scope":"diagnostic only; not added to the complete envelope",
        "FFT_controls":stats,"ordinary_fallback":{"explicitly_enabled":args.ordinary_fallback,
            "marginals_used":len(fallbacks),"first16":fallbacks[:16],
            "all_coordinate_full_indices":[row["coordinate_full"] for row in fallbacks],
            "method":"ordinary positive Cauchy numerator, no clipping or omitted profiles"},
        "joint_missing_for_actual_weight":{"h":missing,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "local_rho_report":local_report,"tests":tests,"field_tests":caps["tests"],
        "residual_rank_receipt":caps["residual_rank_receipt"],
        "sizes_pin":SIZES_PIN,"rho_pin":RHO_PIN,"rank_pin":rho.RANK_PIN,"control_pin":rho.CONTROL_PIN,
        "moment_pin":control.CALIBRATION_PIN,"bank_pin":sizes.BANK_PIN,"broad_pin":sizes.BROAD_PIN,
        "composition_pin":bank.COMPOSITION_PIN,"correction_pin":composition.CORRECTION_PIN,
        "triple_pin":composition.TRIPLE_PIN,"pair_pin":correction.PAIR_PIN,"FFT_pin":cells.FFT_PIN,
        "core_pin":anchor.CORE_PIN,"seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "generated_bytes":0,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
