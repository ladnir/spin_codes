"""Count-aware neighboring-source control, fixed ordered local witness.

Default: one shared center893 table for actual sizes2040..2056. Compare
the old source-only envelope and the count-aware envelope at one ordered
profile, including the nested2047..2049 band. This is NOT band coverage.
--complete uses all31 table cells and the count-aware complete join.
Each actual weight has its own coefficient, normalizer, and source scale.
The exact run count is absorbed with exponent1/256, not alpha/256.
Only z^(-alpha*L) remains outside the join. The encoder is unchanged.

Parent-owned serial numerical jobs;300/1800s,512MiB,20GiB free.
No data files, bytecode, caches, fit, benchmark, or encoder edits.
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

sys.dont_write_bytecode=True
BASE_PIN="e775786311a3840d21a708beb529720bf4a0e32a7a064fc69f038641555c3ce1"
COUNTED_PIN="2c481652ef8a37a14c5bd88ce6af31e0517a3f51706b7a07e3aa6447906ab74b"
MEMORY=512*(1<<20)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen source mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complete",action="store_true")
    parser.add_argument("--rlo",type=int,default=2040)
    parser.add_argument("--rhi",type=int,default=2056)
    parser.add_argument("--center",type=int,default=893)
    parser.add_argument("--bits",type=int,choices=(192,256),default=192)
    parser.add_argument("--seconds",type=int,choices=(300,1800),default=300)
    args=parser.parse_args()
    if not 2032<=args.rlo<=args.rhi<=2064 or args.rhi-args.rlo>=17:
        parser.error("at most17 consecutive actual sizes within2032..2064")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    weights=tuple(range(args.rlo,args.rhi+1))
    started=time.monotonic(); deadline=started+args.seconds
    base=load("cheap_lane_ordered_band_directed.py",BASE_PIN,"counted_base")
    ref=base.load("cheap_lane_ordered_regions_directed.py",base.REFERENCE_PIN,"counted_reference")
    if args.center not in ref.WITNESSES: parser.error("unknown frozen center")
    anchor=ref.load("cheap_lane_cone_cauchy_anchors.py",ref.ANCHOR_PIN,"counted_anchor")
    local=ref.load("cheap_lane_ordered_local_directed.py",ref.LOCAL_PIN,"counted_local")
    integer=ref.load("cheap_lane_integer_regions.py",ref.INTEGER_PIN,"counted_integer")
    join=ref.load("cheap_lane_directed_run_mass.py",ref.JOIN_PIN,"counted_single_join")
    multi=base.load("cheap_lane_band_regions.py",base.MULTI_PIN,"counted_multi")
    oldband=base.load("cheap_lane_band_run_mass.py",base.BAND_PIN,"counted_oldband")
    counted=load("cheap_lane_counted_band_run_mass.py",COUNTED_PIN,"counted_newband")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("count-aware band deadline; no incomplete certificate")
        if peak()>MEMORY: raise MemoryError("count-aware band512MiB limit")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"elapsed_seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    built=local.build_local_pairs(ref.Z,ref.ALPHA,anchor,guard,deadline,384,progress=progress)
    v,core=built["verifier"],built["core"]; built.pop("coefficients",None)
    tests=counted.exact_tests(join,guard)
    scales={}; sources={}; logtwo=v.arb(2).log()
    threshold=v.exact(ref.Z)**v.exact(-ref.ALPHA*core.LIMIT)
    for r in weights:
        h,c,eta=anchor.source_factor(r)
        omega=v.exact(eta)*v.exact(core.BETA)**r
        scales[r]=omega**v.exact(ref.ALPHA/F(256))
        assert (scales[r]**256*threshold).overlaps(
            (omega*v.exact(ref.Z)**(-core.LIMIT))**v.exact(ref.ALPHA))
        sources[str(r)]={"h":h,"c_h":str(c),"eta":str(eta)}
    rows=base.schedule(weights,ref.WITNESSES)
    selected_rows=rows if args.complete else tuple(x for x in rows if x[0]==args.center)
    if not selected_rows: raise ValueError("selected center has no targets")
    regions={r:{} for r in weights}; tables=[]
    for item in selected_rows:
        partial,report=multi.one_table_many(v,built["rational_coefficients"],item,
            ref.WITNESSES[item[0]],weights,integer,args.bits,guard,progress)
        for r,entries in partial.items():
            assert not set(regions[r]).intersection(entries)
            regions[r].update(entries)
        tables.append(report); del partial; gc.collect(); guard()
        progress("counted_band_table_complete",center=item[0],first_full=item[1],last_full=item[2])
    bands=[weights]
    if set((2047,2048,2049)).issubset(weights) and weights!=(2047,2048,2049):
        bands.insert(0,(2047,2048,2049))
    pilots=[]
    for band_weights in bands:
        own={r:regions[r] for r in band_weights}; own_scales={r:scales[r] for r in band_weights}
        counts=oldband.band_counts(band_weights)
        common={f for f in set.intersection(*(set(rows) for rows in own.values()))
                if 2*f<len(counts) and counts[2*f]>0}
        if not common: raise ArithmeticError("no eligible common selected profile")
        sample=895 if 895 in common else min(common,key=lambda f:(abs(f-args.center),f))
        d=2*sample; per_source=[]; individual_sum=v.arb(0)
        for r in band_weights:
            c=r-d
            if not 1<=c<=r: continue
            number=join.run_counts(points=8192,weight=r)[c]
            chain=join.two_orientations(own[r][sample],own[r][sample],diagonal=True)[0]
            value=v.positive(threshold*scales[r]**256*number*chain)
            individual_sum+=value
            per_source.append({"actual_weight":r,"runs":c,"log2_bound_interval":str(value.log()/logtwo)})
        original=oldband.build_envelope(own,own_scales,guard)
        old_chain=join.two_orientations(original[sample],original[sample],diagonal=True)[0]
        old=v.positive(threshold*counts[d]*old_chain)
        new=counted.selected_envelope(own,own_scales,d,join,guard)
        new_region=new["envelope"][sample]
        new_chain=join.two_orientations(new_region,new_region,diagonal=True)[0]
        updated=v.positive(threshold*new["countfactor"]*new_chain)
        assert old>=individual_sum or old.overlaps(individual_sum)
        assert updated>=individual_sum or updated.overlaps(individual_sum)
        pilots.append({"actual_weights":band_weights,"full":[sample,sample],"degree_sum":d,
            "individual_controls":per_source,"sum_individual_log2_interval":str(individual_sum.log()/logtwo),
            "source_only_envelope_log2_interval":str(old.log()/logtwo),
            "count_aware_envelope_log2_interval":str(updated.log()/logtwo),
            "count_aware_loss_vs_individual_sum_log2_interval":str((updated/individual_sum).log()/logtwo),
            "report":new["report"],"not_a_complete_degree_or_band":True})
    result={"scope":"complete actual-source band at35%" if args.complete else "selected count-aware band controls only",
        "whole_code_certificate":False,"all_profiles_in_band_covered":False,"directed":True,
        "actual_weights":weights,"z":str(ref.Z),"alpha":str(ref.ALPHA),"Arb_precision":384,
        "dyadic_bits":args.bits,"tests":tests,"tables":tables,"source_factors":sources,
        "pilots":pilots,"local_coefficients":built["report"],"encoder_changed":False,
        "global_factor":"z^(-alpha*LIMIT) only; source and exact counts already inside join",
        "budget_seconds":args.seconds,"memory_budget_bytes":MEMORY}
    if args.complete:
        for r in weights: assert sorted(regions[r])==list(range(r//2+1))
        joined=counted.counted_join(regions,scales,join,guard,progress=progress)
        total=v.positive(threshold*joined["bound"])
        result.update({"all_profiles_in_band_covered":True,"join":joined["report"],
            "log2_union_bound_interval":str(total.log()/logtwo),
            "band_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
            "band_below_2_pow_minus_60":bool(total<v.arb(2)**(-60))})
    guard()
    result.update({"elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "generated_bytes":0,"base_pin":BASE_PIN,"counted_pin":COUNTED_PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
