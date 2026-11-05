"""Directed neighboring-source band using the fixed r2048 ordered witness.

Default: one center893 pilot for actual sizes2047..2049, including a
same-table singleton2048 control against the frozen regional engine.
--complete extracts every own-weight coefficient while each cyclic256
table is live, then runs the positive source-scaled band envelope once.
The marked polynomials and local caps are shared, NOT normalized entries.

For Omega_r=eta_r*beta^r, scale each R/V by Omega_r^(alpha/256).
There are255 regular matrices and one terminal vector in every chain.
The only remaining global factor is z^(-alpha*L). Group exact indexed
run masses by d=f0+f1=r-runs. No extra profile multiplicity or new alpha.
The2049 endpoint adds(1024,1024); no r2048 endpoint exclusion is reused.

All markers,z,alpha and local caps remain frozen. No fit or encoder change.
Actual-index bands only, at most17 weights within2032..2064. No holes.
300/1200s,512MiB,20GiB free; parent-only serial jobs, stdout, no data files.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode=True
REFERENCE_PIN="1faaf32f2cf1ced0e78fb9952534e1c79b138cf9f530c7c82b82607e650dd78d"
MULTI_PIN="ab74014cbd1e2e6334ab15acaf2fcbc9dc6d9bc97b19e5e6c495f0d2eb8a5852"
BAND_PIN="127a6a82fcb785bcb339115f0e0ea8d951841340ce32540655b3e86b06a20503"
MEMORY=512*(1<<20)
POINTS=8192


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen source mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def schedule(weights,markers):
    centers=tuple(markers)
    assert len(centers)==31 and centers==tuple(sorted(centers))
    full=max(weights)//2
    assigned={center:[] for center in centers}
    for f in range(full+1):
        center=min(centers,key=lambda c:(abs(c-f),c))
        assigned[center].append(f)
    result=[]
    for center,targets in assigned.items():
        if not targets: continue
        assert targets==list(range(targets[0],targets[-1]+1))
        result.append((center,targets[0],targets[-1],"cyclic"))
    assert [f for _,lo,hi,_ in result for f in range(lo,hi+1)]==list(range(full+1))
    return tuple(result)


def rational_endpoint(x,upper):
    if not x.is_finite(): raise ArithmeticError("nonfinite comparison entry")
    return F(str((x.upper() if upper else x.lower()).fmpq()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complete",action="store_true")
    parser.add_argument("--rlo",type=int,default=2047)
    parser.add_argument("--rhi",type=int,default=2049)
    parser.add_argument("--center",type=int,default=893)
    parser.add_argument("--bits",type=int,choices=(128,192,256),default=192)
    parser.add_argument("--seconds",type=int,choices=(300,1200),default=300)
    args=parser.parse_args()
    if not 2032<=args.rlo<=args.rhi<=2064 or args.rhi-args.rlo>=17:
        parser.error("at most17 consecutive actual sizes within2032..2064")
    weights=tuple(range(args.rlo,args.rhi+1))
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    ref=load("cheap_lane_ordered_regions_directed.py",REFERENCE_PIN,"band_reference")
    if args.center not in ref.WITNESSES: parser.error("center must be one of the31 frozen witnesses")
    anchor=ref.load("cheap_lane_cone_cauchy_anchors.py",ref.ANCHOR_PIN,"band_anchor")
    integer=ref.load("cheap_lane_integer_regions.py",ref.INTEGER_PIN,"band_integer")
    local=ref.load("cheap_lane_ordered_local_directed.py",ref.LOCAL_PIN,"band_local")
    join=ref.load("cheap_lane_directed_run_mass.py",ref.JOIN_PIN,"band_single_join")
    multi=load("cheap_lane_band_regions.py",MULTI_PIN,"band_multi_regions")
    band=load("cheap_lane_band_run_mass.py",BAND_PIN,"band_run_envelope")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("band deadline; no incomplete certificate")
        if peak()>MEMORY: raise MemoryError("512MiB band budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"elapsed_seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    built=local.build_local_pairs(ref.Z,ref.ALPHA,anchor,guard,deadline,384,progress=progress)
    v,core=built["verifier"],built["core"]
    built.pop("coefficients",None)
    tests={"integer":integer.exact_tests(guard),"single_join":join.exact_tests(guard),
           "band":band.exact_tests(join,guard),"multi_extraction":multi.exact_tests(guard)}
    rows=schedule(weights,ref.WITNESSES)
    scales={}; source_report={}; logtwo=v.arb(2).log()
    threshold=v.exact(ref.Z)**v.exact(-ref.ALPHA*core.LIMIT)
    for r in weights:
        guard(); h,c,eta=anchor.source_factor(r)
        omega=v.exact(eta)*v.exact(core.BETA)**r
        scales[r]=omega**v.exact(ref.ALPHA/F(256))
        original=(omega*v.exact(ref.Z)**(-core.LIMIT))**v.exact(ref.ALPHA)
        restored=scales[r]**256*threshold
        assert restored.overlaps(original)
        source_report[str(r)]={"h":h,"c_h":str(c),"eta":str(eta),
            "scale_log2_interval":str(scales[r].log()/logtwo),
            "scale256_times_threshold_overlaps_original":True}
    selected_rows=rows if args.complete else tuple(x for x in rows if x[0]==args.center)
    if not selected_rows: raise ValueError("selected cell has no targets in this band")
    regions={r:{} for r in weights}; tables=[]
    for item in selected_rows:
        partial,report=multi.one_table_many(v,built["rational_coefficients"],item,
            ref.WITNESSES[item[0]],weights,integer,args.bits,guard,progress)
        for r,entries in partial.items():
            assert not set(regions[r]).intersection(entries)
            regions[r].update(entries)
        tables.append(report); del partial; gc.collect(); guard()
        progress("band_table_complete",center=item[0],first_full=item[1],last_full=item[2],
                 actual_weights=weights)
    if args.complete:
        for r in weights:
            assert sorted(regions[r])==list(range(r//2+1))
    # On a selected cell, reproduce the old r2048 normalized table exactly.
    # This is an implementation control, not imported certificate data.
    if not args.complete and 2048 in weights:
        old=ref.load("cheap_lane_pair_regions_integer.py",ref.REGION_PIN,"band_old_region_control")
        driver=SimpleNamespace(WEIGHT=2048,MODULUS=256,marker=lambda c:ref.WITNESSES[c])
        old_bits=old.BITS; old.BITS=args.bits; comparisons=0
        try:
            for item in selected_rows:
                rlo,rhi=item[1],min(item[2],1024)
                if rlo>rhi: continue
                old_item=(item[0],rlo,rhi,item[3])
                control,control_report=old.one_table(v,built["rational_coefficients"],old_item,
                    driver,integer,guard,progress)
                for f,old_entry in control.items():
                    for new_row,old_row in zip(regions[2048][f],old_entry):
                        for a,b in zip(new_row,old_row):
                            for upper in (False,True):
                                assert rational_endpoint(a,upper)==rational_endpoint(b,upper)
                                comparisons+=1
                del control,control_report; gc.collect()
        finally:
            old.BITS=old_bits
        tests["frozen2048_regional_endpoint_equalities"]=comparisons
        progress("singleton2048_regional_control_complete",exact_endpoint_equalities=comparisons)
    envelope=band.build_envelope(regions,scales,guard)
    counts=band.band_counts(weights,points=POINTS)
    common={f for f in set.intersection(*(set(entries) for entries in regions.values()))
            if 2*f<len(counts) and counts[2*f]>0}
    if not common: raise ArithmeticError("no common selected marginal")
    sample=895 if 895 in common else min(common,key=lambda f:(abs(f-args.center),f))
    d=2*sample
    controls=[]; summed=v.arb(0)
    for r in weights:
        c=r-d
        if not 1<=c<=r: continue
        number=join.run_counts(points=POINTS,weight=r)[c]
        chain=join.two_orientations(regions[r][sample],regions[r][sample],diagonal=True)[0]
        individual=v.positive(threshold*scales[r]**256*number*chain)
        summed+=individual
        controls.append({"actual_weight":r,"full":[sample,sample],"single":[c,c],"runs":c,
            "log2_selected_full_run_mass_bound_interval":str(individual.log()/logtwo)})
    if d>=len(counts) or counts[d]==0: raise ArithmeticError("selected band degree is empty")
    selected_chain=join.two_orientations(envelope[sample],envelope[sample],diagonal=True)[0]
    selected=v.positive(threshold*counts[d]*selected_chain)
    assert selected>=summed or selected.overlaps(summed)
    pilot={"full":[sample,sample],"degree_sum":d,
        "individual_actual_run_controls":controls,
        "log2_sum_individual_selected_bounds_interval":str(summed.log()/logtwo),
        "log2_selected_full_band_degree_mass_bound_interval":str(selected.log()/logtwo),
        "envelope_loss_vs_sum_log2_interval":str((selected/summed).log()/logtwo),
        "selected_is_not_a_degree_class_or_complete_band":True}
    result={"scope":"all supports and paths in the actual source-size band at35% distance"
                    if args.complete else "selected source-scaled band-envelope pilot only",
        "directed":True,"all_profiles_in_band_covered":False,"whole_code_certificate":False,
        "parameters":{"actual_weights":weights,"z":str(ref.Z),"alpha":str(ref.ALPHA),
            "regional_moduli":[256,256],"dyadic_bits":args.bits,"Arb_precision":384,
            "budget_seconds":args.seconds,"memory_budget_bytes":MEMORY,"complete":args.complete},
        "tests":tests,"tables":tables,"selected_control":pilot,
        "source_scaling":source_report,"source_factors_per_chain":256,
        "global_factor":"z^(-alpha*LIMIT) only; Omega already distributed across R/V",
        "indexed_band_counts_sha256":hashlib.sha256(":".join(map(str,counts)).encode()).hexdigest(),
        "indexed_band_count_total":str(sum(counts)),"all_profiles_counts_outside_alpha":True,
        "weights_have_own_coefficient_and_normalizer":True,
        "no_old_boundary_substitution":True,"ordinary_fallbacks":0,
        "local_coefficients":built["report"],
        "fixed_schedule":[{"center":c,"lo":lo,"hi":hi,"markers":list(map(str,ref.WITNESSES[c]))}
                          for c,lo,hi,_ in rows]}
    if args.complete:
        joined=band.band_join(envelope,weights,join,guard,progress)
        total=v.positive(threshold*joined["bound"])
        if weights==(2047,2048,2049): assert joined["report"]["ordered_envelope_pairs"]==1050625
        worst=max(range(len(joined["degree_bounds"])),key=lambda d:joined["degree_bounds"][d].upper())
        result.update({"all_profiles_in_band_covered":True,
            "log2_union_bound_interval":str(total.log()/logtwo),
            "band_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
            "band_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
            "join":joined["report"],"worst_degree_sum":{"d":worst,
                "log2_interval":str((threshold*joined["degree_bounds"][worst]).log()/logtwo)}})
    guard()
    result.update({"elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "generated_bytes":0,"reference_pin":REFERENCE_PIN,"multi_pin":MULTI_PIN,"band_pin":BAND_PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
