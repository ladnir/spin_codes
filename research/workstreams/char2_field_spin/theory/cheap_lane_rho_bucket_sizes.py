"""Complete r2048 FLOAT screen with residual-rank returns and rho buckets.

Use one coherent raw/powered local cap map at z23/100,alpha13/128. The
frozen bucket builder sums disjoint first-occupancy corrections before
subtracting, retaining the old9/dense33 and degree-Jensen ceilings.
Only its dense33_bucket_sum table enters the complete-size replay.

The original scaled fifteen-center bank gains center768. An explicit
--fitted-center768 option assigns that cell the previously tested rational
markers; all other cells retain the frozen curve. Every target keeps its
OWN profile normalizer, even when sharing a marker. This is a fixed witness,
not an optimizer. Every marginal and ordered profile remains in the frozen
complete indexed-run-mass envelope, including ordinary Cauchy fallbacks.

Six pair doublings retain H64/H63J,4096 route slots and32767 updates plus
the output-only terminal. Source Omega2048 is applied once, counts outside
alpha. FLOAT reconnaissance only, not a certificate or whole-code result.
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
FITTED768=(F(320890,682431),F(251362,1026933))
MEMORY=512*(1<<20)
LN2=math.log(2)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modulus",type=int,choices=(256,512),default=256)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--fitted-center768",action="store_true")
    parser.add_argument("--ordinary-fallback",action=argparse.BooleanOptionalAction,default=True)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    sizes=load("cheap_lane_pair_sizes_screen.py",SIZES_PIN,"rho_sizes_complete")
    rho=load("cheap_lane_rho_bucket_control.py",RHO_PIN,"rho_sizes_buckets")
    rank=rho.load("cheap_lane_residual_rank_return.py",rho.RANK_PIN,"rho_sizes_rank")
    bank=sizes.load("cheap_lane_pair_bank_screen.py",sizes.BANK_PIN,"rho_sizes_bank")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"rho_sizes_broad")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"rho_sizes_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"rho_sizes_triple")
    pair=correction.load_pair()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"rho_sizes_neighbor")
    peel=broad.load("cheap_lane_profile_peeling.py",broad.PEEL_PIN,"rho_sizes_peel")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"rho_sizes_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("complete rho screen deadline; no automatic retry")
        if peak()>MEMORY: raise MemoryError("512MiB complete rho screen budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    centers=tuple(sorted(set(broad.default_centers(WEIGHT))|{768}))
    assert len(centers)==16 and centers[0]>0 and centers[-1]<WEIGHT//2
    tests=broad.exact_tests(np,screen,peel)
    tests["rho"]=rho.exact_tests(pair,correction,np,guard)
    tests["residual_rank"]=rank.exact_tests(guard)
    tests["FFT_exact_toy_coefficients"]=fft.exact_toy(np,guard)
    tests["scale_checks"]=cells.exact_scale_tests()
    caps=rank.augment_caps(composition.build_caps(Z,ALPHA,triple,anchor,guard,deadline),Z,ALPHA,guard)
    verifier,core=caps["verifier"],caps["core"]
    tables,local_report=rho.bucket_pair_logs(caps["powered"],float(ALPHA),pair,old,
        correction,np,guard,caps["raw"],progress)
    selected=tables["dense33_bucket_sum"]
    del tables,caps["raw"],caps["powered"]
    progress("rho_complete_coefficients_ready",actual_weight=WEIGHT,z=str(Z),alpha=str(ALPHA),
        centers=centers,paired_modulus=args.modulus,fitted_center768=args.fitted_center768)
    original_markers=sizes.markers
    def fixed_markers(weight,center):
        if args.fitted_center768 and weight==WEIGHT and center==768: return FITTED768
        return original_markers(weight,center)
    # Private imported module only; restore even if a numerical guard fails.
    sizes.markers=fixed_markers
    try:
        regions,bank_report,stats,fallbacks=sizes.marginals(selected,WEIGHT,centers,args.modulus,
            args.ordinary_fallback,broad,cells,fft,screen,np,guard,progress)
    finally:
        sizes.markers=original_markers
    del selected
    assert len(regions)==1025
    assert regions[768][3:]==(FITTED768 if args.fitted_center768 else original_markers(WEIGHT,768))
    missing,c_h,eta=anchor.source_factor(WEIGHT)
    log_source=float((verifier.exact(ALPHA)*(verifier.exact(eta).log()
        +WEIGHT*verifier.exact(core.BETA).log()-core.LIMIT*verifier.exact(Z).log())).mid())
    progress("rho_complete_marginals_ready",marginals=len(regions),ordinary_fallbacks=len(fallbacks))
    # This control charges all run512 sources to ONE selected profile. It is
    # not a run cover and is never added to the complete envelope below.
    numerator=8192*comb(WEIGHT-1,511)*comb(8192-WEIGHT-1,511)
    run_mass,remainder=divmod(numerator,512); assert remainder==0
    forward,reverse=screen.paired_chains(regions[768],regions[768])
    assert abs(forward-reverse)<1e-7
    selected_control={"actual_full":[768,768],"single":[512,512],"runs":512,
        "markers":list(map(str,regions[768][3:])),
        "uses_ordinary_fallback":any(row["coordinate_full"]==768 for row in fallbacks),
        "log2_selected_full_run_mass_charged_FLOAT":(log_source+math.log(run_mass)+forward)/LN2,
        "scope":"selected-profile control only; not added to complete envelope"}
    progress("rho_selected_anchor_control",**selected_control)
    result=sizes.complete_join(WEIGHT,regions,fallbacks,broad,screen,np,log_source,guard,progress)
    guard()
    print(json.dumps({"scope":"complete actual r2048 FLOAT residual-rank/rho-bucket paired-bank run-mass envelope",
        "directed":False,"negative_is_not_a_certificate":True,"whole_code_certificate":False,
        "other_weights_covered":False,"parameter_optimization":False,
        "parameters":{"actual_weight":WEIGHT,"z":str(Z),"alpha":str(ALPHA),
            "paired_modulus":args.modulus,"centers":centers,"fitted_center768":args.fitted_center768,
            "fitted_marker_witness":list(map(str,FITTED768)) if args.fitted_center768 else None,
            "curve_full":str(sizes.CF),"curve_single":str(sizes.CS),
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "pair_bound":"dense33_bucket_sum: max of valid correction lower bounds; same degree-Jensen ceiling",
        "result":result,"selected_anchor_control":selected_control,"bank":bank_report,"FFT_controls":stats,
        "ordinary_fallback":{"explicitly_enabled":args.ordinary_fallback,"marginals_used":len(fallbacks),
            "first16":fallbacks[:16],"all_coordinate_full_indices":[x["coordinate_full"] for x in fallbacks],
            "method":"ordinary positive Cauchy numerator; no coefficient clipping or profile omission"},
        "joint_missing_for_actual_weight":{"h":missing,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "local_rho_report":local_report,"tests":tests,"field_tests":caps["tests"],
        "residual_rank_receipt":caps["residual_rank_receipt"],
        "sizes_pin":SIZES_PIN,"rho_pin":RHO_PIN,"rank_pin":rho.RANK_PIN,
        "bank_pin":sizes.BANK_PIN,"broad_pin":sizes.BROAD_PIN,"composition_pin":bank.COMPOSITION_PIN,
        "correction_pin":composition.CORRECTION_PIN,"triple_pin":composition.TRIPLE_PIN,
        "pair_pin":correction.PAIR_PIN,"FFT_pin":cells.FFT_PIN,"core_pin":anchor.CORE_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
