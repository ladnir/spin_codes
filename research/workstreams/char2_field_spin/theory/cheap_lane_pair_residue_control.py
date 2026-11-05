"""Complete r1024 FLOAT residue control at one requested interior modulus.

Pin the successful mixed-bank witness. Keep its z, alpha, original fifteen
marker centers, raw caps, pair corrections, exact indexed source counts,
and final output-only chronology unchanged. Only the paired interior
residue modulus is selected from128/256/512. The old-cone f0..6 boundary
always retains exact full degree and singleton modulus512.

Each regional tuple carries its own denominator. The frozen complete
join retains all263168 ordered profiles, including the full residual run
mass after an exact-count budget expires. An unresolved interior residue
uses its explicitly reported ordinary bound; none is silently dropped.
Directed boundary exports do not make the subsequent FLOAT join directed.

No parameter search, generated files, or frozen edits. One invocation:
300s,256MiB,20GiB free, stdout only. The512 default is the old control.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
MIXED_PIN="15608a5d808df9a180ef7c74fb6c741c69d58c4d998d80b26d2fc60791098927"
Z,ALPHA=F(304397809,625000000),F(119,1280)
WEIGHT,SECONDS,MEMORY=1024,300,256*(1<<20)


def load_mixed():
    path=Path(__file__).with_name("cheap_lane_pair_bank_boundary_screen.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=MIXED_PIN:
        raise RuntimeError("frozen mixed-bank hash mismatch")
    spec=importlib.util.spec_from_file_location("residue_control_mixed",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modulus",type=int,choices=(128,256,512),default=512)
    parser.add_argument("--count-seconds",type=float,default=30)
    parser.add_argument("--max-queries",type=int,default=20000)
    args=parser.parse_args()
    if not 0<=args.count_seconds<=120 or not 0<=args.max_queries<=100000:
        parser.error("count-seconds0..120; max-queries0..100000")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    mixed=load_mixed()
    bank=mixed.load("cheap_lane_pair_bank_screen.py",mixed.BANK_PIN,"residue_control_bank")
    boundary=mixed.load("cheap_lane_old_boundary_regions.py",mixed.BOUNDARY_PIN,"residue_control_boundary")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"residue_control_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"residue_control_triple")
    pair=correction.load_pair(); small=triple.load_small()
    calibration=small.load_calibration(); broad=calibration.load_broad()
    centers=tuple(broad.default_centers(WEIGHT))
    assert len(centers)==15 and 0 not in centers and 512 not in centers
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"residue_control_neighbor")
    peel=broad.load("cheap_lane_profile_peeling.py",broad.PEEL_PIN,"residue_control_peel")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"residue_control_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("residue-control deadline; no complete bound or retry")
        if peak()>MEMORY: raise MemoryError("256MiB residue-control budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    counts=screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"residue_control_counts")
    counts.DEADLINE=deadline
    tests=broad.exact_tests(np,screen,peel)
    tests["count_small_checks"]=counts.checks()
    tests["pair_and_correction_checks"]=correction.exact_tests(pair,np,guard)+pair.exact_tests(np,guard)
    tests["FFT_exact_toy_coefficients"]=fft.exact_toy(np,guard)
    tests["scale_checks"]=cells.exact_scale_tests()
    caps=composition.build_caps(Z,ALPHA,triple,anchor,guard,deadline)
    verifier,core=caps["verifier"],caps["core"]
    tables,single,controls=correction.corrected_pair_logs(caps["powered"],float(ALPHA),
        tuple(range(-16,17,4)),pair,old,np,guard,raw=caps["raw"])
    del caps["raw"],caps["powered"],single
    selected=tables["concavity_Jensen_minimum"]; del tables
    progress("combined_pair_coefficients_ready",paired_interior_modulus=args.modulus,
        boundary_singleton_modulus=512,centers=centers,controls=controls)
    regions,bank_report,stats,fallbacks,assignment=bank.marginals(selected,centers,args.modulus,True,
        broad,cells,fft,screen,np,guard,progress)
    control=bank.marginal_control(selected,regions,assignment,fallbacks,args.modulus,pair,
        calibration,broad,cells,fft,screen,np,guard)
    del selected
    before_fallbacks=list(fallbacks)
    patch,boundary_report=boundary.old_boundary_regions(Z,ALPHA,anchor,screen,guard,deadline,
        precision=256,check_singleton=True)
    assert sorted(patch)==list(range(7))
    assert boundary_report["z"]==str(Z) and boundary_report["alpha"]==str(ALPHA)
    assert boundary_report["singleton_modulus"]==512 and boundary_report["full_limit"]==6
    assert not boundary_report["source_factor_applied"] and not boundary_report["support_count_applied"]
    for f,region in patch.items(): regions[f]=region
    fallbacks=[row for row in fallbacks if row["full"] not in patch]
    assert len(regions)==513 and all(row["full"]>=7 for row in fallbacks)
    for row in bank_report:
        row["old_boundary_replacements"]=[f for f in patch if row["first_full"]<=f<=row["last_full"]]
        row["ordinary_fallbacks_after_patch"]=sum(x["center"]==row["center"] for x in fallbacks)
    progress("directed_boundary_replacements_ready",paired_interior_modulus=args.modulus,
        boundary_singleton_modulus=512,patched_full_indices=list(range(7)),
        remaining_fallback_full_indices=[x["full"] for x in fallbacks],
        exact_f0_control=boundary_report["exact_singleton_control"])
    h,c_h,eta=anchor.source_factor(WEIGHT); exponent=verifier.exact(ALPHA)
    log_source=float((exponent*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()
                               -core.LIMIT*verifier.exact(Z).log())).mid())
    report,diagnostics=bank.join_with_diagnostics(regions,fallbacks,counts,broad,screen,np,log_source,
        args.max_queries,args.count_seconds,deadline,guard,progress)
    for rows in diagnostics["queried_weighted_leaders_by_numerical_method"].values():
        for row in rows: row["old_boundary_full_indices"]=[f for f in row["full"] if f in patch]
    diagnostics["resolved_bucket_scope"]="paired interior FFT plus directed old-boundary exports, not ordinary fallbacks"
    guard()
    print(json.dumps({"scope":"complete r1024 MIXED regional FLOAT residue-control envelope",
        "directed_complete_bound":False,"whole_code_certificate":False,
        "negative_is_not_a_certificate":True,"other_weights_covered":False,
        "parameters":{"weight":WEIGHT,"z":str(Z),"alpha":str(ALPHA),
            "paired_interior_modulus":args.modulus,"boundary_singleton_modulus":512,
            "paired_centers":centers,"seconds":SECONDS,"count_seconds":args.count_seconds,
            "max_queries":args.max_queries},
        "comparison_scope":{"only_interior_residue_modulus_changes":True,
            "same_caps_markers_source_counts_and_terminal":True,"calibration_performed":False,
            "default512_replays_frozen_mixed_algorithm":True},
        "result":report,"diagnostic_separation":diagnostics,"bank":bank_report,
        "regional_provenance":{"old_boundary_full_indices":list(range(7)),
            "paired_interior_full_range_inclusive":[7,512],"target_own_denominators_retained":True},
        "ordinary_fallback":{"explicitly_enabled":True,"before_patch":len(before_fallbacks),
            "before_full_indices":[x["full"] for x in before_fallbacks],"marginals_used":len(fallbacks),
            "first16":fallbacks[:16],"all_fallback_full_indices":[x["full"] for x in fallbacks]},
        "boundary_report":boundary_report,"FFT_controls_for_pair_bank":stats,
        "frozen_shared_center_control":control,"controls":controls,"tests":tests,
        "field_tests":caps["tests"],"block_caps":caps["block_report"],
        "joint_missing":{"h":h,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "mixed_pin":MIXED_PIN,"bank_pin":mixed.BANK_PIN,"boundary_pin":mixed.BOUNDARY_PIN,
        "composition_pin":bank.COMPOSITION_PIN,"count_pin":screen.COUNT_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
