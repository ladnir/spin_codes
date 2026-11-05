"""Complete actual-r4096 FLOAT screen with the current ordered pair operator.

Build composition -> residual-rank returns -> ordered129 pair correction
once, then fit a fixed31-center bank and retain every actual matching
profile. Default z1/16,alpha119/1280 matches the historical size4096
diagnostic; explicit rational overrides select ONE witness, not a search.
Each center receives at most twelve frozen coordinate-fit evaluations.

Centers are the doubled r2048 grid, with first center2 instead of0. The
upper endpoint uses fit target2047.5/s1 only to propose positive markers.
All actual targets use their own integer f,s=4096-2f normalizer. There is
no hole reversal at r4096. The frozen generic complete join retains all
4,198,400 ordered profiles and exact run masses summing to C(8192,4096).

The f0 boundary is extracted directly: total singleton degree4096 in64
pair factors forces degree(0,64) in EVERY factor. Thus its pair-majorant
regional coefficients are H_(0,64)^64 and H_(0,64)^63 J_(0,64), with
normalizer1 and no marker. This is an exact projection evaluated in FLOAT,
not a directed certificate. Only f0 and its possible fallback are replaced.

Other unresolved residues use reported ordinary Cauchy bounds if enabled;
no coefficient is clipped or profile omitted. The constructor, routing,
source factor Omega4096, and output-only terminal remain distinct stages.
No new alpha is applied after local pair construction. Parent execution
only:300s default,512MiB,20GiB free,stdout only,zero generated files.
No numerical imports or jobs on import; no frozen helper is edited.
"""

import argparse
from collections import defaultdict
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
ORDERED_PIN="f816addd50969a83ad951567944d905938ece966c9260175e3b4a08b1be53255"
WEIGHT,MEMORY=4096,512*(1<<20)
SMALL_CENTERS=(0,3,11,25,44,69,98,131,169,211,256,304,354,406,458,512,
    566,618,670,720,768,813,855,893,926,955,980,999,1013,1021,1024)
CENTERS=(2,)+tuple(2*x for x in SMALL_CENTERS[1:])
LN2=math.log(2)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def fit_target(center):
    full=F(4095,2) if center==2048 else F(center)
    return full,WEIGHT-2*full


def singleton_boundary(log_coefficients,screen):
    """Marker-free leading singleton coefficient of the chosen pair bound."""
    row=tuple(float(x) for x in log_coefficients[0,64])
    def normalize(values):
        if any(math.isnan(x) or x==math.inf for x in values):
            raise ArithmeticError("invalid boundary log coefficient")
        scale=max(values)
        if not math.isfinite(scale): raise ArithmeticError("zero boundary matrix object")
        entries=tuple(0.0 if x==-math.inf else math.exp(x-scale) for x in values)
        if any(math.isfinite(x) and y==0 for x,y in zip(values,entries)):
            raise ArithmeticError("positive boundary coefficient underflow")
        return screen.normalized_matrix(entries,scale)
    h,j=normalize(row[:4]),normalize(row[4:])
    r=screen.matrix_power(h,64)
    terminal=screen.matrix_vector(screen.matrix_power(h,63),j)
    region=(r,terminal,0.0,F(1),F(1))
    forward,reverse=screen.paired_chains(region,region)
    direct=screen.first_log(screen.matrix_vector(screen.matrix_power(h,16383),j))
    if not math.isfinite(direct) or max(abs(forward-direct),abs(reverse-direct))>2e-5:
        raise ArithmeticError("boundary regional/whole finite-chain disagreement")
    return region,{"actual_full":0,"actual_single":4096,"pair_degree":[0,64],
        "regular":"H_(0,64)^64","terminal":"H_(0,64)^63 J_(0,64)",
        "exact_projection_of_selected_majorant":True,"arithmetic":"FLOAT",
        "marker_powers_applied":False,"log_normalizer":0,"normalizer_exact":1,
        "extra_orientation_factor":False,"indexed_sources_at_profile_0_0":2,
        "finite_chain_control_difference_bits_FLOAT":abs(forward-direct)/LN2}


def exact_boundary_tests(guard):
    """Positive rational top-degree products, including terminal order."""
    def product(a,b):
        out=[]
        for row in a:
            result=[]
            for j in range(len(b[0])):
                entry=defaultdict(F)
                for k,left in enumerate(row):
                    for (f,s),x in left.items():
                        for (ff,ss),y in b[k][j].items(): entry[f+ff,s+ss]+=x*y
                result.append(dict(entry))
            out.append(tuple(result))
        return tuple(out)
    def matrix_product(a,b):
        return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                           for j in range(len(b[0]))) for i in range(len(a)))
    checks=0
    for slots in (1,2):
        h=(( {(0,slots):F(1,3),(1,0):F(2)}, {(0,slots):F(2,5),(0,0):F(7)} ),
           ( {(0,slots):F(3,7)}, {(0,slots):F(5,11),(0,0):F(13)} ))
        j=(( {(0,slots):F(7,13),(0,0):F(17)}, ),( {(0,slots):F(11,19),(1,0):F(23)}, ))
        top_h=tuple(tuple(x[0,slots] for x in row) for row in h)
        top_j=tuple(tuple(x[0,slots] for x in row) for row in j)
        rr,tt=h,j; er,et=top_h,top_j
        for step in range(2):
            guard(); tt=product(rr,tt); rr=product(rr,rr)
            et=matrix_product(er,et); er=matrix_product(er,er)
            degree=slots*(1<<(step+1)); marker=F(3,5)
            for actual,expected in ((rr,er),(tt,et)):
                for i,row in enumerate(actual):
                    for k,p in enumerate(row):
                        assert p.get((0,degree),F(0))==expected[i][k]
                        assert all(f+s<=degree for f,s in p)
                        assert p[0,degree]*marker**degree/(comb(degree,degree)*marker**degree)==expected[i][k]
                        checks+=3
        # The terminal side must be left-multiplied, not a transposed row.
        assert matrix_product(top_h,top_j)[0][0] != matrix_product(
            ((top_j[0][0],top_j[1][0]),),top_h)[0][0]; checks+=1
    assert len(CENTERS)==31 and tuple(sorted(set(CENTERS)))==CENTERS
    assert fit_target(2)==(F(2),F(4092)) and fit_target(2048)==(F(4095,2),F(1))
    assert comb(4096,0)*comb(4096,4096)==1
    assert 8192*comb(4095,4095)*comb(4095,4095)//4096==2
    checks+=4
    return {"exact_rational_boundary_checks":checks,"positive_top_degree_extraction":True,
        "noncommuting_terminal":True,"normalizer_and_indexed_endpoint_count":True,
        "fixed_doubled_grid_with_first_center2":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(1,16))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--modulus",type=int,choices=(256,512),default=256)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--ordinary-fallback",action=argparse.BooleanOptionalAction,default=True)
    args=parser.parse_args()
    if not 0<args.z<1 or not 0<args.alpha<=1: parser.error("require0<z<1,0<alpha<=1")
    if max(abs(n).bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    sizes=load("cheap_lane_pair_sizes_screen.py",SIZES_PIN,"large_sizes")
    rho=load("cheap_lane_rho_bucket_control.py",RHO_PIN,"large_coarse")
    ordered=load("cheap_lane_rho_ordered_buckets.py",ORDERED_PIN,"large_ordered")
    control=rho.load("cheap_lane_pair_profile_calibration.py",rho.CONTROL_PIN,"large_control")
    rank=rho.load("cheap_lane_residual_rank_return.py",rho.RANK_PIN,"large_rank")
    calibration=control.load("cheap_lane_broad_profile_calibration.py",control.CALIBRATION_PIN,"large_moments")
    bank=sizes.load("cheap_lane_pair_bank_screen.py",sizes.BANK_PIN,"large_bank")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"large_broad")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"large_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"large_triple")
    pair=correction.load_pair()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"large_neighbor")
    peel=broad.load("cheap_lane_profile_peeling.py",broad.PEEL_PIN,"large_peel")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"large_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("ordered large-screen deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("512MiB ordered large-screen budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    tests=exact_boundary_tests(guard)
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
    selected=tables["ordered129_bucket_sum"]
    del tables,caps["raw"],caps["powered"]
    endpoint,endpoint_report=singleton_boundary(selected,screen)
    objective=control.Objective(selected,0.0,calibration,np,guard)
    witnesses={}; fit_reports=[]
    for center in CENTERS:
        guard(); fstar,sstar=fit_target(center)
        objective.targets=256*np.asarray([float(fstar),float(sstar)])
        start=sizes.markers(WEIGHT,fstar)
        best,history,status=control.coordinate_fit(objective,np.asarray(list(map(broad.logq,start))),np,guard)
        marker=tuple(max(F(1,1<<30),F(math.exp(float(x))).limit_denominator(1<<20)) for x in best["x"])
        rounded=objective.evaluate(np.asarray(list(map(broad.logq,marker))))
        witnesses[center]=marker
        report={"center":center,"fit_targets":[str(fstar),str(sstar)],"markers":list(map(str,marker)),
            "evaluations":len(history),"status":status,
            "ordinary_objective_gain_bits_FLOAT":history[0]["log2_ordinary_FLOAT"]-rounded["objective"]/LN2,
            "rounded_region_means_FLOAT":[float(x/256) for x in rounded["means"]],
            "rounded_per_region_log_marker_gradient_bits_FLOAT":[float(x/(256*LN2)) for x in rounded["gradient"]]}
        fit_reports.append(report); progress("ordered_large_center_fitted",**report)
    del objective
    original_markers=sizes.markers
    def calibrated_markers(weight,center):
        assert weight==WEIGHT and center in witnesses
        return witnesses[center]
    sizes.markers=calibrated_markers
    try:
        # Permit the inherited routine to finish f0 even if its residue is
        # unreliable; replace it independently, then enforce the caller's
        # fallback policy on all remaining actual marginals.
        regions,bank_report,stats,all_fallbacks=sizes.marginals(selected,WEIGHT,CENTERS,args.modulus,
            True,broad,cells,fft,screen,np,guard,progress)
    finally:
        sizes.markers=original_markers
    del selected
    assert len(regions)==2049
    for f,region in enumerate(regions):
        center=min(CENTERS,key=lambda c:(abs(c-f),c))
        assert region[3:]==witnesses[center]
    removed=[row for row in all_fallbacks if row["coordinate_full"]==0]
    fallbacks=[row for row in all_fallbacks if row["coordinate_full"]!=0]
    regions[0]=endpoint
    if fallbacks and not args.ordinary_fallback:
        raise ArithmeticError("unresolved non-boundary residues; ordinary fallback disabled")
    endpoint_report["replaced_f0_fallbacks"]=removed
    assert regions[0][2:]==(0.0,F(1),F(1))
    missing,c_h,eta=anchor.source_factor(WEIGHT)
    assert missing==16 and c_h==F(1,17)
    log_source=float((verifier.exact(args.alpha)*(verifier.exact(eta).log()
        +WEIGHT*verifier.exact(core.BETA).log()-core.LIMIT*verifier.exact(args.z).log())).mid())
    endpoint_chain=screen.paired_chains(endpoint,endpoint)[0]
    endpoint_report["two_alternating_sources_log2_bound_FLOAT"]=(log_source+math.log(2)+endpoint_chain)/LN2
    progress("ordered_large_marginals_ready",marginals=len(regions),ordinary_fallbacks=len(fallbacks),
        exact_endpoint=endpoint_report)
    result=sizes.complete_join(WEIGHT,regions,fallbacks,broad,screen,np,log_source,guard,progress)
    assert result["ordered_profiles"]==4198400 and result["actual_weight"]==WEIGHT
    guard()
    print(json.dumps({"scope":"complete actual r4096 FLOAT current-ordered-operator run-mass envelope",
        "directed":False,"negative_is_not_a_certificate":True,"whole_code_certificate":False,
        "other_weights_covered":False,"proof_parameter_search":False,"marker_fit_only":True,
        "parameters":{"actual_weight":WEIGHT,"hole_coordinates":False,"z":str(args.z),"alpha":str(args.alpha),
            "paired_modulus":args.modulus,"centers":CENTERS,"max_fit_evaluations":control.MAX_PROPOSALS,
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "pair_bound":"current residual-rank returns and ordered129; N4/coarse/Jensen ceilings retained",
        "fits":fit_reports,"bank":bank_report,"result":result,"exact_projection_boundary":endpoint_report,
        "FFT_controls":stats,"ordinary_fallback":{"explicitly_enabled":args.ordinary_fallback,
            "marginals_used":len(fallbacks),"first16":fallbacks[:16],
            "all_coordinate_full_indices":[row["coordinate_full"] for row in fallbacks],
            "bank_raw_fallback_count_before_exact_boundary":len(all_fallbacks),
            "method":"ordinary positive Cauchy numerator; no clipping or omitted profiles"},
        "fit_scope":"ordinary repeated-profile proposals only; no projected stationarity claim",
        "endpoint_proxy_changes_actual_degrees":False,"source_factor_applied_once":True,
        "joint_missing_for_actual_weight":{"h":missing,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "local_ordered_report":local_report,"tests":tests,"field_tests":caps["tests"],
        "residual_rank_receipt":caps["residual_rank_receipt"],"sizes_pin":SIZES_PIN,"rho_pin":RHO_PIN,
        "ordered_pin":ORDERED_PIN,"rank_pin":rho.RANK_PIN,"control_pin":rho.CONTROL_PIN,
        "moment_pin":control.CALIBRATION_PIN,"composition_pin":bank.COMPOSITION_PIN,
        "correction_pin":composition.CORRECTION_PIN,"triple_pin":composition.TRIPLE_PIN,
        "pair_pin":correction.PAIR_PIN,"FFT_pin":cells.FFT_PIN,"core_pin":anchor.CORE_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
