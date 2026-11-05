"""Directed r2048 ordered-ratio witness, with bounded integer regional replay.

Use the frozen round28 rational marker bank at z23/100,alpha23/256.
No fit/search occurs here. The new local constructor exports coherent
directed pair coefficient upper bounds for the actual refresh distribution.
Reuse the audited integer regional engine: cyclic256 by256, six doublings,
H64/H63J, own integer target normalizers. The final run join is explicitly
called with weight2048,points8192; source applied once, counts outside alpha.

Default is one center893 pilot, with dyadic192 and a256-bit control.
--complete covers all1025 marginals/1050624 ordered profiles, with no old
r1024 boundary substitution or numerical fallback. No incomplete cover
can claim a source-size certificate. All exact input markers are fixed.
300/1200s,512MiB,20GiB free, stdout only; no coefficient files or caches.
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
LOCAL_PIN="1ed1d59d9beed36b9a1bd1a2bc65af043efe342e95e2afd5047f90c696c04c47"
REGION_PIN="a2a2a6d784b5ae23f26b93e2c4e3caec87f7667598f8bd76d70ab3b49614e9a0"
INTEGER_PIN="9dc583ae54caed0a78f71e409cc37aa398713b419f65538370f10c3f876f40bf"
JOIN_PIN="54f2ff06a1476ccf84fa4db1b757f1d5f8e4e814a1f0f733dc73ccf18b2ddabb"
ANCHOR_PIN="319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
WEIGHT,POINTS,MODULUS=2048,8192,256
Z,ALPHA=F(23,100),F(23,256)
MEMORY=512*(1<<20)
WITNESSES={
    0: (F("255/646883"), F("623632/483565")),
    3: (F("2419/1024041"), F("293785/228638")),
    11: (F("8764/1015915"), F("76375/60142")),
    25: (F("16268/835553"), F("990191/796033")),
    44: (F("19097/562552"), F("1188291/982718")),
    69: (F("20840/396259"), F("488307/419335")),
    98: (F("13149/178492"), F("2571/2308")),
    131: (F("81591/841483"), F("237530/224483")),
    169: (F("27031/219901"), F("995648/999729")),
    211: (F("86273/572830"), F("942700/1014063")),
    256: (F("14896/83141"), F("497696/577845")),
    304: (F("145616/698575"), F("268421/339217")),
    354: (F("145441/611765"), F("436473/605084")),
    406: (F("98552/369131"), F("566040/868703")),
    458: (F("292881/992593"), F("333124/569733")),
    512: (F("212837/658769"), F("338729/653828")),
    566: (F("246817/705174"), F("255105/561733")),
    618: (F("122742/327335"), F("268248/679163")),
    670: (F("352629/883663"), F("268101/793189")),
    720: (F("106361/252384"), F("138583/485930")),
    768: (F("461297/1043126"), F("96565/408826")),
    813: (F("458463/994147"), F("185747/968869")),
    855: (F("14701/30731"), F("58540/386667")),
    893: (F("470924/954065"), F("75529/651815")),
    926: (F("45263/89353"), F("48961/571044")),
    955: (F("457383/883348"), F("47125/788193")),
    980: (F("483733/917311"), F("29909/791045")),
    999: (F("450692/843147"), F("674/31575")),
    1013: (F("422884/783397"), F("3992/427051")),
    1021: (F("119281/219747"), F("649/255264")),
    1024: (F("96022/176593"), F("432/1020359")),
}


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen source mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def schedule():
    centers=tuple(WITNESSES)
    assert len(centers)==31 and centers==tuple(sorted(centers))
    assigned={center:[] for center in centers}
    for f in range(WEIGHT//2+1):
        center=min(centers,key=lambda c:(abs(c-f),c))
        assigned[center].append(f)
    result=[]
    for center,targets in assigned.items():
        assert targets==list(range(targets[0],targets[-1]+1))
        assert all(x>0 for x in WITNESSES[center])
        result.append((center,targets[0],targets[-1],"cyclic"))
    assert [f for _,lo,hi,_ in result for f in range(lo,hi+1)]==list(range(1025))
    return tuple(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complete",action="store_true")
    parser.add_argument("--center",type=int,choices=tuple(WITNESSES),default=893)
    parser.add_argument("--bits",type=int,choices=(128,192,256),default=192)
    parser.add_argument("--control-bits",type=int,choices=(128,192,256),default=256)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    parser.add_argument("--seconds",type=int,choices=(300,1200),default=300)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    anchor=load("cheap_lane_cone_cauchy_anchors.py",ANCHOR_PIN,"ordered_regions_anchor")
    region=load("cheap_lane_pair_regions_integer.py",REGION_PIN,"ordered_regions_reuse")
    integer=load("cheap_lane_integer_regions.py",INTEGER_PIN,"ordered_regions_integer")
    join=load("cheap_lane_directed_run_mass.py",JOIN_PIN,"ordered_regions_join")
    local=load("cheap_lane_ordered_local_directed.py",LOCAL_PIN,"ordered_regions_local")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("ordered regional deadline; no incomplete certificate")
        if peak()>MEMORY: raise MemoryError("512MiB ordered regional budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"elapsed_seconds":time.monotonic()-started,
                          "peak_process_bytes":peak(),**values}),flush=True)
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]
    built.pop("coefficients",None)
    tests={"integer":integer.exact_tests(guard),"join":join.exact_tests(guard)}
    counts=join.run_counts(points=POINTS,weight=WEIGHT)
    tests["actual2048_run_count_sum"]=str(sum(counts))
    rows=schedule(); tests["all1025_integer_marginals_in_fixed_schedule"]=True
    driver=SimpleNamespace(WEIGHT=WEIGHT,MODULUS=MODULUS,marker=lambda c:WITNESSES[c])
    h,c_h,eta=anchor.source_factor(WEIGHT)
    source=(v.exact(eta)*v.exact(core.BETA)**WEIGHT*v.exact(Z)**(-core.LIMIT))**v.exact(ALPHA)
    logtwo=v.arb(2).log()
    runs=(args.bits,) if args.complete or args.control_bits==args.bits else (args.bits,args.control_bits)
    reports=[]; chosen=[]
    original_bits=region.BITS
    regions=[None]*1025
    try:
        for bits in runs:
            region.BITS=bits
            for item in rows if args.complete else tuple(x for x in rows if x[0]==args.center):
                guard()
                partial,report=region.one_table(v,built["rational_coefficients"],item,driver,
                                                integer,guard,progress)
                report["dyadic_mantissa_bits"]=bits
                reports.append(report)
                if args.complete:
                    for f,entry in partial.items():
                        assert regions[f] is None
                        regions[f]=entry
                sample=895 if 895 in partial else item[0]
                c=WEIGHT-2*sample
                if 1<=c<=WEIGHT:
                    forward,reverse=join.two_orientations(partial[sample],partial[sample],diagonal=True)
                    selected=v.positive(source*counts[c]*forward)
                    chosen.append({"bits":bits,"full":[sample,sample],"singles":[c,c],"runs":c,
                        "markers":list(map(str,WITNESSES[item[0]])),
                        "log2_selected_whole_run_mass_bound_interval":str(selected.log()/logtwo),
                        "selected_is_not_a_run_or_source_size_cover":True})
                progress("ordered_integer_cell_complete",center=item[0],bits=bits,
                         first_full=item[1],last_full=item[2])
                del partial; gc.collect()
    finally:
        region.BITS=original_bits
    result={"scope":"all actual2048-group source supports and paths at35% distance" if args.complete
                    else "selected directed regional-cell and precision control only",
        "directed":True,"all_profiles_at_weight_covered":False,"whole_code_certificate":False,
        "parameters":{"weight":WEIGHT,"points":POINTS,"z":str(Z),"alpha":str(ALPHA),
            "regional_moduli":[256,256],"Arb_precision":args.precision,
            "dyadic_mantissas":runs,"complete":args.complete,
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "tables":reports,"selected_controls":chosen,"tests":tests,
        "local_coefficients":built["report"],
        "fixed_schedule":[{"center":c,"lo":lo,"hi":hi,"method":method,
                           "markers":list(map(str,WITNESSES[c]))} for c,lo,hi,method in rows],
        "source_log2_interval":str(source.log()/logtwo),
        "joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "historical_boundary_substitution":False,"ordinary_fallbacks":0,
        "source_once_counts_outside_alpha":True}
    if args.complete:
        assert all(entry is not None for entry in regions)
        joined=join.directed_join(regions,guard,progress,weight=WEIGHT,points=POINTS)
        assert joined["report"]["ordered_profiles"]==1050624
        total=v.positive(source*joined["bound"])
        worst=max(range(1,WEIGHT+1),key=lambda c:joined["run_bounds"][c-1].upper())
        result.update({"all_profiles_at_weight_covered":True,
            "log2_union_bound_interval":str(total.log()/logtwo),
            "actual2048_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
            "actual2048_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
            "join":joined["report"],
            "worst_run":{"runs":worst,
                "log2_interval":str((source*joined["run_bounds"][worst-1]).log()/logtwo)}})
    guard()
    result.update({"elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "generated_bytes":0,"local_pin":LOCAL_PIN,"regional_pin":REGION_PIN,
        "integer_pin":INTEGER_PIN,"join_pin":JOIN_PIN,"anchor_pin":ANCHOR_PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()

