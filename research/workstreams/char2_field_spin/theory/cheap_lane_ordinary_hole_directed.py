"""Directed ordinary-Cauchy cover of actual6144-group supports at35%.

Freeze the current hole-screen rational markers, z1/64 and alpha119/1280.
Reverse the complete local H/J degree table once; evaluate ordinary H64
and H63J, retaining own hole-target normalizers. No FFT, fitting, fallback
or boundary replacement. The directed join uses hole weight2048 only for
its identical indexed support counts and matching coordinates; the source
factor remains Omega6144. This is not a whole-code certificate.

Parent only, -B,300s,512MiB,20GiB free,stdout only,no generated files.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time
sys.dont_write_bytecode=True
ORDINARY_PIN="cb2cc15b9209a95cdfd57b51d98baed99aa205a82099595075c2e72563b826ed"
HOLE_PIN="cb9e6201b3fda3d7598b66339be99ba4584c8c87524cd7400d96089f87320854"
SIZES_PIN="7a7024b78912abb78eaaef9e90f766eefdd7e7f807d8dfe6143f38c86c233afd"
WEIGHT,HOLES,MEMORY=6144,2048,512*(1<<20)
Z,ALPHA=F(1,64),F(119,1280)
BANK=(
    (0,0,1,F("119/626497"),F("180827/206245")),
    (3,2,7,F("461/405003"),F("856819/980830")),
    (11,8,18,F("2610/627829"),F("166773/193159")),
    (25,19,34,F("4343/462848"),F("76283/90192")),
    (44,35,56,F("12450/760951"),F("100427/122138")),
    (69,57,83,F("25103/990435"),F("713115/900532")),
    (98,84,114,F("27919/786600"),F("431687/569873")),
    (131,115,150,F("29811/638452"),F("738756/1026797")),
    (169,151,190,F("34753/587603"),F("60363/89159")),
    (211,191,233,F("26084/360429"),F("39443/62438")),
    (256,234,280,F("59788/695717"),F("474727/811574")),
    (304,281,329,F("94661/948996"),F("294557/548539")),
    (354,330,380,F("58110/512261"),F("499340/1021089")),
    (406,381,432,F("110131/867498"),F("145649/330118")),
    (458,433,485,F("10579/75693"),F("399346/1010081")),
    (512,486,539,F("46069/302382"),F("289918/828949")),
    (566,540,592,F("166043/1010864"),F("268185/876302")),
    (618,593,644,F("183320/1046893"),F("218090/820873")),
    (670,645,695,F("182029/981865"),F("190399/839051")),
    (720,696,744,F("160038/821659"),F("153233/801887")),
    (768,745,790,F("53797/264560"),F("149525/946639")),
    (813,791,834,F("105761/501203"),F("117657/919498")),
    (855,835,874,F("207211/951078"),F("31499/312316")),
    (893,875,909,F("23208/103685"),F("79822/1035953")),
    (926,910,940,F("185584/811011"),F("37659/661667")),
    (955,941,967,F("241811/1037423"),F("36049/909790")),
    (980,968,989,F("221438/935691"),F("25217/1007951")),
    (999,990,1006,F("240087/1003264"),F("5444/385937")),
    (1013,1007,1017,F("62871/260636"),F("2235/362189")),
    (1021,1018,1022,F("203507/839890"),F("1561/930660")),
    (1024,1023,1024,F("253611/1045226"),F("243/870170")),
)

def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free")
    started=time.monotonic();deadline=started+args.seconds
    ordinary=load("cheap_lane_ordered_ordinary_directed.py",ORDINARY_PIN,"hole_ordinary")
    hole=load("cheap_lane_ordered_hole_screen.py",HOLE_PIN,"hole_screen_pinned")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"hole_anchor")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"hole_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"hole_normalizers")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"hole_join")
    sizes=large.load("cheap_lane_pair_sizes_screen.py",SIZES_PIN,"hole_sizes")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"hole_broad")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline:raise TimeoutError("hole ordinary deadline; no partial cover")
        if peak()>MEMORY:raise MemoryError("hole ordinary512MiB exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    tests={"ordinary":ordinary.exact_tests(guard),
        "complement":sizes.complement_tests(broad,guard),"join":join.exact_tests(guard)}
    assert [f for _,lo,hi,_,_ in BANK for f in range(lo,hi+1)]==list(range(1025))
    assert comb(8192,WEIGHT)==comb(8192,HOLES)
    tests["complete_bank_and_count_complement"]=True
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,384,progress=progress)
    v,core=built["verifier"],built["core"];logtwo=v.arb(2).log()
    coefficients={(64-f-s,s):row for (f,s),row in built["coefficients"].items()}
    assert len(coefficients)==len(built["coefficients"])==2145
    assert all(row is built["coefficients"][64-f-s,s] for (f,s),row in coefficients.items())
    regions=[None]*1025
    for center,lo,hi,u,w in BANK:
        guard()
        h,j=ordinary.evaluate_pairs(v,coefficients,(u,w),guard)
        regular=v.matpow(h,64);terminal=v.matmul(v.matpow(h,63),j)
        numerator=tuple(tuple(v.positive(x) for row in matrix for x in row)
            for matrix in (regular,terminal))
        for f in range(lo,hi+1):
            guard();s=HOLES-2*f;actual_full=WEIGHT-4096+f
            denominator=extraction.target_denominator(HOLES,f,(u,w))
            assert comb(4096,f)*comb(4096-f,s)==comb(4096,actual_full)*comb(4096-actual_full,s)
            dd=v.exact(denominator)
            regions[f]=tuple(tuple(v.positive(x/dd) for x in row) for row in numerator)
        progress("ordinary_hole_cell_complete",center=center,lo=lo,hi=hi)
    assert all(x is not None for x in regions)
    h,c_h,eta=anchor.source_factor(WEIGHT);assert h==24 and c_h==F(1,25)
    source=(v.exact(eta)*v.exact(core.BETA)**WEIGHT*v.exact(Z)**(-core.LIMIT))**v.exact(ALPHA)
    joined=join.directed_join(regions,guard,progress,weight=HOLES,points=8192)
    assert joined["report"]["ordered_profiles"]==1050624
    total=v.positive(source*joined["bound"])
    worst=max(range(1,HOLES+1),key=lambda c:joined["run_bounds"][c-1].upper())
    guard()
    print(json.dumps({"scope":"all actual6144-group source supports and all histories at35%",
        "directed_arithmetic":True,"all_profiles_at_weight_covered":True,"whole_code_certificate":False,
        "parameters":{"actual_weight":WEIGHT,"coordinate_weight":HOLES,"z":str(Z),"alpha":str(ALPHA),
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,
            "memory_budget_bytes":MEMORY},
        "regional_method":"complete degree reversal then ordinary H64/H63J with own hole normalizer",
        "degree_reversals":1,"source_factor_weight":WEIGHT,"source_once_counts_outside_alpha":True,
        "boundary_replacement":False,"residue_projection":False,"ordinary_fallbacks":0,
        "log2_union_bound_interval":str(total.log()/logtwo),
        "actual6144_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
        "worst_run":{"runs":worst,"log2_interval":str((source*joined["run_bounds"][worst-1]).log()/logtwo)},
        "source_log2_interval":str(source.log()/logtwo),
        "joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "join":joined["report"],"tests":tests,"local_coefficients":built["report"],
        "bank":[{"center":c,"lo":lo,"hi":hi,"u":str(u),"v":str(w)} for c,lo,hi,u,w in BANK],
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "nonnegative_intersections":v.intersections,"generated_bytes":0,
        "pins":{"ordinary":ORDINARY_PIN,"hole_screen":HOLE_PIN,"sizes":SIZES_PIN,"local":large.LOCAL_PIN,
            "join":large.JOIN_PIN,"normalizers":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)

if __name__=="__main__": main()

