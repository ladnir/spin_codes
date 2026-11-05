"""Selected r2048 equal-profile FLOAT calibration with frozen rho buckets.

Default full898, singles252, run count252. Every comparison charges the
SAME exact indexed run mass, not an exact selected-profile count. This is
not a cover of that run, other profiles, or the whole source weight.

Build residual-rank return caps and the dense33 disjoint-left bucket
correction once, coherently at the chosen rational z/alpha. Compare the
nearest original fifteen-bank marker, the target's own curve, and at most
twelve ordinary-objective proposals, then replay positive rational markers
at m256/m512. No projected stationarity or numerical certificate is claimed.
The fixed finite chain is H16383J/d256, preserving32767 updates and the
output-only terminal. Source is applied once; run counts stay outside alpha.

Only private imported control FULL/SINGLE/RUNS fields are set, before
constructing Objective or calling replay. The frozen control.main (with
its historical center704 assertion) is never called. No frozen edits.
180s/512MiB/20GiB, stdout only, no generated files or numerical work on import.
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
RHO_PIN="57391fb3ab3ada8078c17a98eb921dbb5777452422c429960e3fa47aa07624b8"
WEIGHT=2048
SECONDS,MEMORY,LN2=180,512*(1<<20),math.log(2)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full",type=int,default=898)
    parser.add_argument("--z",type=F,default=F(23,100))
    parser.add_argument("--alpha",type=F,default=F(13,128))
    args=parser.parse_args()
    if not 1<=args.full<=1023: parser.error("require full in1..1023")
    if not 0<args.z<1 or not 0<args.alpha<=1: parser.error("require0<z<1,0<alpha<=1")
    if max(abs(n).bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    rho=load("cheap_lane_rho_bucket_control.py",RHO_PIN,"rho_profile_buckets")
    control=rho.load("cheap_lane_pair_profile_calibration.py",rho.CONTROL_PIN,"rho_profile_control")
    rank=rho.load("cheap_lane_residual_rank_return.py",rho.RANK_PIN,"rho_profile_rank")
    sizes=control.load("cheap_lane_pair_sizes_screen.py",control.SIZES_PIN,"rho_profile_sizes")
    calibration=control.load("cheap_lane_broad_profile_calibration.py",control.CALIBRATION_PIN,"rho_profile_moments")
    bank=sizes.load("cheap_lane_pair_bank_screen.py",sizes.BANK_PIN,"rho_profile_bank")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"rho_profile_broad")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"rho_profile_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"rho_profile_triple")
    pair=correction.load_pair()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"rho_profile_neighbor")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"rho_profile_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("180-second rho-profile deadline; no automatic retry")
        if peak()>MEMORY: raise MemoryError("512MiB rho-profile budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    f=args.full; s=WEIGHT-2*f; c=s
    assert control.WEIGHT==WEIGHT and 1<=c<=min(WEIGHT,8192-WEIGHT)
    control.FULL,control.SINGLE,control.RUNS=f,s,c
    tests={"rho":rho.exact_tests(pair,correction,np,guard),
        "residual_rank":rank.exact_tests(guard),
        "positive_moments":calibration.exact_tests(np,screen,old),
        "FFT_exact_toys":fft.exact_toy(np,guard),"scale_checks":cells.exact_scale_tests()}
    caps=rank.augment_caps(composition.build_caps(args.z,args.alpha,triple,anchor,guard,deadline),
        args.z,args.alpha,guard)
    verifier,core=caps["verifier"],caps["core"]
    tables,local_report=rho.bucket_pair_logs(caps["powered"],float(args.alpha),pair,old,
        correction,np,guard,caps["raw"],progress)
    selected=tables["dense33_bucket_sum"]
    del tables,caps["raw"],caps["powered"]
    missing,c_h,eta=anchor.source_factor(WEIGHT)
    log_source=float((verifier.exact(args.alpha)*(verifier.exact(eta).log()
        +WEIGHT*verifier.exact(core.BETA).log()-core.LIMIT*verifier.exact(args.z).log())).mid())
    numerator=8192*comb(WEIGHT-1,c-1)*comb(8192-WEIGHT-1,c-1)
    run_mass,remainder=divmod(numerator,c); assert remainder==0
    source_count=log_source+math.log(run_mass)
    constant=source_count-256*math.log(comb(4096,f)*comb(4096-f,s))
    objective=control.Objective(selected,constant,calibration,np,guard)
    assert tuple(objective.targets)==(256*f,256*s)
    centers=tuple(broad.default_centers(WEIGHT))
    center=min(centers,key=lambda x:(abs(x-f),x))
    inherited=sizes.markers(WEIGHT,center); own=sizes.markers(WEIGHT,f)
    best,history,status=control.coordinate_fit(objective,np.asarray(list(map(broad.logq,own))),np,guard)
    fitted=tuple(max(F(1,1<<30),F(math.exp(float(x))).limit_denominator(1<<20)) for x in best["x"])
    progress("rho_profile_markers_fixed",full=f,single=s,runs=c,nearest_original_center=center,
        proposals=len(history),status=status,fitted_markers=list(map(str,fitted)))
    results=[]
    for label,markers in (("inherited_original_bank",inherited),("own_full_curve",own),("fitted_rational",fitted)):
        direct=objective.evaluate(np.asarray(list(map(broad.logq,markers))))
        replays=[]
        for modulus in control.MODULI:
            row=control.replay(selected,markers,modulus,source_count,broad,cells,fft,screen,np,guard)
            difference=abs(row["ordinary"]["log2_selected_run_mass_bound_FLOAT"]-direct["objective"]/LN2)
            assert difference<2e-5
            row.update({"modulus":modulus,"ordinary_two_path_control_bits":difference})
            replays.append(row)
            progress("rho_profile_replay",witness=label,markers=list(map(str,markers)),**row)
        results.append({"witness":label,"markers":list(map(str,markers)),
            "ordinary_log_marker_gradient_FLOAT":[float(x/LN2) for x in direct["gradient"]],
            "tilted_region_means_FLOAT":[float(x/256) for x in direct["means"]],
            "ordinary_covariance_FLOAT":direct["covariance"].tolist(),"replays":replays})
    guard()
    print(json.dumps({"scope":"selected actual r2048 equal profile only, charged its entire indexed run mass",
        "directed":False,"negative_is_not_a_certificate":True,"complete_weight_bound":False,
        "run_class_cover":False,"whole_code_certificate":False,
        "parameters":{"weight":WEIGHT,"fulls":[f,f],"singles":[s,s],"runs":c,
            "z":str(args.z),"alpha":str(args.alpha),"moduli":control.MODULI,
            "inherited_centers":centers,"inherited_center":center,"max_proposals":control.MAX_PROPOSALS,
            "seconds":SECONDS,"memory_budget_bytes":MEMORY},
        "pair_bound":"residual-rank returns; dense33_bucket_sum with old/dense/Jensen ceilings",
        "fixed_source":{"log2_factor_FLOAT":log_source/LN2,"log2_run_mass":math.log(run_mass)/LN2,
            "exact_run_mass_sha256":hashlib.sha256(str(run_mass).encode()).hexdigest(),
            "same_count_in_every_comparison":True,"exact_profile_queries":0,
            "joint_missing_h":missing,"c_h":str(c_h),"eta":str(eta)},
        "fit":{"method":"at most12 analytic coordinate proposals, ordinary objective only",
            "status":status,"evaluations":len(history),"history":history,
            "ordinary_stationarity_is_not_projected_stationarity":True},
        "results":results,"local_rho_report":local_report,"tests":tests,
        "field_tests":caps["tests"],"residual_rank_receipt":caps["residual_rank_receipt"],
        "rho_pin":RHO_PIN,"control_pin":rho.CONTROL_PIN,"rank_pin":rho.RANK_PIN,
        "sizes_pin":control.SIZES_PIN,"moment_pin":control.CALIBRATION_PIN,
        "composition_pin":bank.COMPOSITION_PIN,"correction_pin":composition.CORRECTION_PIN,
        "triple_pin":composition.TRIPLE_PIN,"pair_pin":correction.PAIR_PIN,
        "FFT_pin":cells.FFT_PIN,"core_pin":anchor.CORE_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
