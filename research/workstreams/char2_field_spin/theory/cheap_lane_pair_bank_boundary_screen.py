"""Complete r1024 FLOAT scan with an explicit old-cone boundary patch.

The fixed original fifteen-center bank supplies refined/corrected paired
regions. Replace f0..6 by independently normalized directed regional
upper bounds from the OLD cone, with exact full degree through6 and
singleton residue512. Both constructions bound the same underlying
regional coefficient at the same z and alpha. No source factor or count
is included in either regional construction.

This is intentionally a MIXED regional witness, not a claim that every
region uses the refined pair caps. The first seven regions use their own
v marker and denominator. Remove only their ordinary-fallback records;
retain and label any remaining fallback. The unchanged frozen join sums
all263168 ordered profiles and preserves every unqueried source via its
exact residual run mass. Directed boundary arithmetic does not make the
subsequent FLOAT matrix/count join a certificate.

No scans, calibration, coefficient files, or frozen source edits.
300s default (225 optional),256MiB,20GiB free; stdout only.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
BANK_PIN="2b5099dc1794a97a2c72c4d09584a62ee248b4982b9889bc791e01e5bd685742"
BOUNDARY_PIN="fdb3094bda1f319827c7f50d63cd68720c460cc66498f54575e81c5e2e518698"
WEIGHT,MEMORY=1024,256*(1<<20)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency hash mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--seconds",type=int,choices=(225,300),default=300)
    parser.add_argument("--count-seconds",type=float,default=30)
    parser.add_argument("--max-queries",type=int,default=20000)
    parser.add_argument("--boundary-precision",type=int,choices=(256,384),default=256)
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
    bank=load("cheap_lane_pair_bank_screen.py",BANK_PIN,"patched_pair_bank")
    boundary=load("cheap_lane_old_boundary_regions.py",BOUNDARY_PIN,"old_boundary_patch")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"patched_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"patched_triple")
    pair=correction.load_pair(); small=triple.load_small()
    calibration=small.load_calibration(); broad=calibration.load_broad()
    centers=tuple(broad.default_centers(WEIGHT))
    assert len(centers)==15 and 0 not in centers and 512 not in centers
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"patched_neighbor")
    peel=broad.load("cheap_lane_profile_peeling.py",broad.PEEL_PIN,"patched_peel_tests")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"patched_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("patched-bank deadline; no complete bound or retry")
        if peak()>MEMORY: raise MemoryError("256MiB patched-bank budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    counts=screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"patched_counts")
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
    selected=tables["concavity_Jensen_minimum"]; del tables
    progress("combined_pair_coefficients_ready",controls=controls,centers=centers,
        planned_old_boundary_full_indices=list(range(7)),ordinary_fallback_explicitly_enabled=True)
    regions,bank_report,stats,fallbacks,assignment=bank.marginals(selected,centers,512,True,
        broad,cells,fft,screen,np,guard,progress)
    paired_control=bank.marginal_control(selected,regions,assignment,fallbacks,512,pair,
        calibration,broad,cells,fft,screen,np,guard)
    del selected
    original_fallbacks=list(fallbacks)
    patch,boundary_report=boundary.old_boundary_regions(args.z,args.alpha,anchor,screen,guard,
        deadline,precision=args.boundary_precision,check_singleton=True)
    assert sorted(patch)==list(range(7))
    assert boundary_report["z"]==str(args.z) and boundary_report["alpha"]==str(args.alpha)
    assert not boundary_report["source_factor_applied"] and not boundary_report["support_count_applied"]
    replacement_rows=[]
    for f,new in patch.items():
        guard(); previous=regions[f]
        # Diagnostic only: a replacement may be tighter or looser. Both
        # are valid regional upper bounds, with independent normalizers.
        ratios=[]
        for (a,la),(b,lb) in zip(new[:2],previous[:2]):
            ratios.extend((math.log(x)+la-new[2]-math.log(y)-lb+previous[2])/math.log(2)
                          for x,y in zip(a,b))
        replacement_rows.append({"full":f,"single":WEIGHT-2*f,
            "previous_method":"ordinary_fallback" if any(x["full"]==f for x in fallbacks) else "paired_residue",
            "new_method":"old_cone_exact_full_le6_singleton_cyclic512",
            "new_minus_previous_entry_log2_range_FLOAT":[min(ratios),max(ratios)],
            "new_markers":[str(new[3]),str(new[4])]})
        regions[f]=new
    fallbacks=[x for x in fallbacks if x["full"] not in patch]
    assert all(x["full"]>=7 for x in fallbacks) and len(regions)==513
    for row in bank_report:
        row["old_boundary_replacements"]=[f for f in patch if row["first_full"]<=f<=row["last_full"]]
        row["ordinary_fallbacks_after_patch"]=sum(x["center"]==row["center"] for x in fallbacks)
    progress("old_boundary_replacements_ready",replacements=replacement_rows,
        previous_fallback_full_indices=[x["full"] for x in original_fallbacks],
        remaining_fallback_full_indices=[x["full"] for x in fallbacks],
        boundary_tests=boundary_report["tests"],exact_f0_control=boundary_report["exact_singleton_control"])
    h,c_h,eta=anchor.source_factor(WEIGHT); exponent=verifier.exact(args.alpha)
    log_source=float((exponent*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()
                               -core.LIMIT*verifier.exact(args.z).log())).mid())
    report,diagnostics=bank.join_with_diagnostics(regions,fallbacks,counts,broad,screen,np,log_source,
        args.max_queries,args.count_seconds,deadline,guard,progress)
    # The join's historical 'resolved_marginals' bucket now includes the
    # directed old-boundary exports; make that distinction explicit.
    for rows in diagnostics["queried_weighted_leaders_by_numerical_method"].values():
        for row in rows: row["old_boundary_full_indices"]=[f for f in row["full"] if f in patch]
    diagnostics["resolved_bucket_scope"]="paired FFT regions plus directed old-boundary exports, never ordinary fallbacks"
    guard()
    print(json.dumps({"scope":"complete r1024 FLOAT scan with MIXED old-boundary and refined/corrected-pair regional bounds",
        "directed_complete_bound":False,"negative_is_not_a_certificate":True,"whole_code_certificate":False,
        "other_weights_covered":False,"fresh_parameter_or_marker_calibration":False,
        "parameters":{"weight":WEIGHT,"z":str(args.z),"alpha":str(args.alpha),"modulus":512,
            "paired_centers":centers,"seconds":args.seconds,"count_seconds":args.count_seconds,
            "max_queries":args.max_queries,"boundary_precision":args.boundary_precision},
        "result":report,"diagnostic_separation":diagnostics,"bank":bank_report,
        "regional_provenance":{"old_boundary_full_indices":list(range(7)),
            "paired_full_indices":[7,512],"paired_range_is_inclusive":True,
            "same_z_alpha_source_and_count":True,"target_own_denominators_retained":True},
        "boundary_report":boundary_report,"boundary_replacement_controls":replacement_rows,
        "ordinary_fallback":{"explicitly_enabled":True,"before_patch":len(original_fallbacks),
            "before_full_indices":[x["full"] for x in original_fallbacks],
            "marginals_used":len(fallbacks),"first16":fallbacks[:16],
            "all_fallback_full_indices":[x["full"] for x in fallbacks]},
        "FFT_controls_for_pair_bank":stats,"frozen_shared_center_control":paired_control,
        "joint_missing":{"h":h,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "controls":controls,"tests":tests,"field_tests":caps["tests"],"block_caps":caps["block_report"],
        "bank_pin":BANK_PIN,"boundary_pin":BOUNDARY_PIN,"composition_pin":bank.COMPOSITION_PIN,
        "count_pin":screen.COUNT_PIN,"seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
