"""Directed replay of one fixed current-operator witness at actual weight4096.

The local constructor supplies rational upper coefficients for H=M*M and
J=M*T at z1/16,alpha119/1280. Each fixed marker table uses the positive
integer cyclic256 engine. Six doublings produce H^64 and H^63 J. Each
actual target(f,4096-2f) has its own multinomial and marker denominator.
No local power, source factor, or support count is reapplied there.

The full-degree-zero endpoint is independent of the cyclic calculation.
The target(0,4096) forces every pair factor to have degree(0,64), since
each local pair degree satisfies f+s<=64. Its regional objects are
H_(0,64)^64 and H_(0,64)^63 J_(0,64), evaluated from exact rational local
coefficients with directed Arb arithmetic. Its normalizer is1; no marker
or additional orientation factor appears. Only marginal f0 is replaced.

Default execution builds center2048 and this endpoint. It reports both
ordered chains(2048,0)/(0,2048), each charged with the exact full run2048
mass. These selected channels do not bound the entire run or source size.
--complete builds all31 frozen marker cells and joins every2049 marginal
and4,198,400 ordered profiles. The exact run masses cover C(8192,4096).
Omega4096 and the output threshold are applied once, outside that join.

Numerical receipts are indexed separately; source alone is not evidence.
Import performs no job.
Parent execution only, using -B:300/1200s,512MiB,20GiB free, stdout only.
There is no calibration, fallback, generated coefficient file, or claim
about other source sizes or the whole code.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
LOCAL_PIN="1ed1d59d9beed36b9a1bd1a2bc65af043efe342e95e2afd5047f90c696c04c47"
REGION_PIN="ab74014cbd1e2e6334ab15acaf2fcbc9dc6d9bc97b19e5e6c495f0d2eb8a5852"
INTEGER_PIN="9dc583ae54caed0a78f71e409cc37aa398713b419f65538370f10c3f876f40bf"
JOIN_PIN="54f2ff06a1476ccf84fa4db1b757f1d5f8e4e814a1f0f733dc73ccf18b2ddabb"
ANCHOR_PIN="319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
SCREEN_PIN="1a0655fde36dbb4f090a143998833bf8ccaf8c25f029df42b1b8837a1db670ed"
WEIGHT,POINTS,MODULUS,MEMORY=4096,8192,256,512*(1<<20)
Z,ALPHA=F(1,16),F(119,1280)
WITNESSES={
    2:(F("1122471/826715"),F("875755627/369801")),
    6:(F("758343/558581"),F("252205960/320183")),
    22:(F("1291737/951776"),F("129345153/607261")),
    50:(F("533307/393113"),F("91709801/993070")),
    88:(F("761933/561807"),F("34419869/669028")),
    138:(F("1109402/818029"),F("21556591/674266")),
    196:(F("1289879/950714"),F("13165028/602819")),
    262:(F("1244493/917411"),F("14080387/892076")),
    338:(F("1379794/1015271"),F("5804643/494740")),
    422:(F("1130641/829657"),F("6693913/747494")),
    512:(F("1321807/966349"),F("3853289/551113")),
    608:(F("984580/716421"),F("2092151/377731")),
    708:(F("508593/367999"),F("4635407/1043131")),
    812:(F("1297014/932449"),F("326287/90889")),
    916:(F("1322979/944614"),F("1164073/397509")),
    1024:(F("1127399/799042"),F("1482464/622447")),
    1132:(F("823769/579478"),F("331974/171377")),
    1236:(F("1151639/804284"),F("1525871/965500")),
    1340:(F("1125734/780597"),F("1295747/1014554")),
    1440:(F("1440390/992077"),F("857333/836248")),
    1536:(F("1458177/998003"),F("758412/933113")),
    1626:(F("919362/625603"),F("323062/508581")),
    1710:(F("403679/273250"),F("263913/543679")),
    1786:(F("496339/334403"),F("292893/810671")),
    1852:(F("1550843/1040668"),F("112427/430349")),
    1910:(F("871429/582689"),F("21432/119959")),
    1960:(F("363569/242348"),F("103615/932231")),
    1998:(F("1265131/841250"),F("26368/425425")),
    2026:(F("127417/84567"),F("4246/157873")),
    2042:(F("1318522/874133"),F("7309/1004494")),
    2048:(F("646355/428343"),F("53/87651")),
}


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def schedule():
    """Nearest fixed center, with ties assigned to the lower center."""
    centers=tuple(WITNESSES)
    assert len(centers)==31 and centers==tuple(sorted(set(centers)))
    assigned={center:[] for center in centers}
    for f in range(WEIGHT//2+1):
        assigned[min(centers,key=lambda c:(abs(c-f),c))].append(f)
    result=[]
    for center,targets in assigned.items():
        assert targets==list(range(targets[0],targets[-1]+1))
        assert all(value>0 for value in WITNESSES[center])
        result.append((center,targets[0],targets[-1],"cyclic"))
    assert [f for _,lo,hi,_ in result for f in range(lo,hi+1)]==list(range(2049))
    return tuple(result)


def singleton_boundary(v,degrees,join,guard):
    """Directed exact top-coefficient projection, with no cyclic aliases."""
    guard()
    if any(min(f,s)<0 or f+s>64 or len(row)!=6 for (f,s),row in degrees.items()):
        raise ValueError("expected the local pair total-degree downset")
    row=degrees[0,64]
    if any(not isinstance(x,(F,int)) or isinstance(x,bool) or x<0 for x in row):
        raise ValueError("exact nonnegative rational boundary channels required")
    h=[[v.exact(row[0]),v.exact(row[1])],[v.exact(row[2]),v.exact(row[3])]]
    j=[[v.exact(row[4])],[v.exact(row[5])]]
    r=v.matpow(h,64); terminal=v.matmul(v.matpow(h,63),j)
    region=tuple(tuple(v.positive(x) for line in matrix for x in line) for matrix in (r,terminal))
    assert tuple(map(len,region))==(4,2)
    assert comb(4096,0)*comb(4096,4096)==1
    forward,reverse=join.two_orientations(region,region,diagonal=True)
    direct=v.positive(v.matmul(v.matpow(h,16383),j)[0][0])
    if not forward.overlaps(direct) or not reverse.overlaps(direct):
        raise ArithmeticError("endpoint regional/direct finite-chain intervals disagree")
    guard()
    return region,{"actual_full":0,"actual_single":4096,"pair_degree":[0,64],
        "regular":"H_(0,64)^64","terminal":"H_(0,64)^63 J_(0,64)",
        "coefficient_projection_exact":True,"arithmetic":"directed Arb from exact rational local upper caps",
        "normalizer_exact":1,"marker_powers_applied":False,"residue_aliases":False,
        "extra_orientation_factor":False,"direct_chain_exponent":16383,
        "finite_chain_overlap":True,"indexed_sources_at_profile_0_0":2}


def exact_tests(screen,region,join,guard):
    """Frozen rational endpoint/extraction toys plus this fixed schedule."""
    tests={"boundary":screen.exact_boundary_tests(guard),
        "extraction":region.exact_tests(guard),"ordered_join":join.exact_tests(guard)}
    rows=schedule()
    assert tuple(WITNESSES)==screen.CENTERS
    assert rows[-1][0]==2048 and rows[-1][1]<=2048<=rows[-1][2]
    assert WEIGHT==64*64 and 256*64-1==16383
    tests["own_schedule_and_boundary_dimensions"]=True
    return tests


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complete",action="store_true")
    parser.add_argument("--bits",type=int,choices=(128,192,256),default=192)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    parser.add_argument("--seconds",type=int,choices=(300,1200),default=300)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    anchor=load("cheap_lane_cone_cauchy_anchors.py",ANCHOR_PIN,"large_directed_anchor")
    local=load("cheap_lane_ordered_local_directed.py",LOCAL_PIN,"large_directed_local")
    region=load("cheap_lane_band_regions.py",REGION_PIN,"large_directed_regions")
    integer=load("cheap_lane_integer_regions.py",INTEGER_PIN,"large_directed_integer")
    join=load("cheap_lane_directed_run_mass.py",JOIN_PIN,"large_directed_join")
    screen=load("cheap_lane_ordered_large_screen.py",SCREEN_PIN,"large_directed_boundary_tests")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("large directed deadline; no incomplete cover")
        if peak()>MEMORY: raise MemoryError("512MiB large directed budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    tests=exact_tests(screen,region,join,guard)
    tests["integer_arithmetic"]=integer.exact_tests(guard)
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]
    built.pop("coefficients",None)
    counts=join.run_counts(points=POINTS,weight=WEIGHT)
    assert sum(counts)==comb(POINTS,WEIGHT) and counts[WEIGHT]==2
    tests["all_actual4096_run_counts"]=True
    h,c_h,eta=anchor.source_factor(WEIGHT)
    assert h==16 and c_h==F(1,17)
    source=(v.exact(eta)*v.exact(core.BETA)**WEIGHT*v.exact(Z)**(-core.LIMIT))**v.exact(ALPHA)
    logtwo=v.arb(2).log()
    boundary,boundary_report=singleton_boundary(v,built["rational_coefficients"],join,guard)
    rows=schedule(); regions=[None]*2049; reports=[]
    chosen=rows if args.complete else tuple(item for item in rows if item[0]==2048)
    for item in chosen:
        guard()
        partial,report=region.one_table_many(v,built["rational_coefficients"],item,
            WITNESSES[item[0]],(WEIGHT,),integer,args.bits,guard,progress)
        own=partial[WEIGHT]
        assert tuple(own)==tuple(range(item[1],item[2]+1))
        for f,entry in own.items():
            assert regions[f] is None
            regions[f]=entry
        reports.append(report)
        progress("large_directed_cell_complete",center=item[0],first_full=item[1],last_full=item[2])
        del partial,own; gc.collect()
    # Only f0 changes. The local pair source/count/alpha remain common to
    # cyclic and exact endpoint projections of the same regional object.
    regions[0]=boundary
    forward,reverse=join.two_orientations(regions[2048],regions[0])
    selected=[]
    for full,value in (((2048,0),forward),((0,2048),reverse)):
        charged=v.positive(source*counts[2048]*value)
        selected.append({"full":full,"singles":[WEIGHT-2*f for f in full],"runs":2048,
            "log2_conditional_chain_interval":str(value.log()/logtwo),
            "log2_selected_chain_with_full_run_mass_interval":str(charged.log()/logtwo),
            "full_run_count_decimal_digits":len(str(counts[2048])),
            "not_a_run_or_source_size_cover":True})
    endpoint_chain=join.two_orientations(boundary,boundary,diagonal=True)[0]
    boundary_report["log2_profile_0_0_bound_interval"]=str((source*2*endpoint_chain).log()/logtwo)
    boundary_report["profile_0_0_diagnostic_not_added_to_join"]=True
    result={"scope":"all actual4096-group source supports and all histories at35% distance" if args.complete
                    else "selected center2048 and exact full-degree-zero endpoint only",
        "directed_arithmetic":True,"all_profiles_at_weight_covered":False,
        "whole_code_certificate":False,"ordinary_fallbacks":0,
        "parameters":{"weight":WEIGHT,"points":POINTS,"z":str(Z),"alpha":str(ALPHA),
            "regional_moduli":[MODULUS,MODULUS],"requested_Arb_precision":args.precision,
            "working_Arb_precision":built["report"]["working_precision"],
            "integer_mantissa_bits":args.bits,"complete":args.complete,
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "tables":reports,"selected_leading_channels":selected,"boundary":boundary_report,
        "tests":tests,"local_coefficients":built["report"],
        "fixed_schedule":[{"center":center,"lo":lo,"hi":hi,"method":kind,
            "markers":list(map(str,WITNESSES[center]))} for center,lo,hi,kind in rows],
        "source_log2_interval":str(source.log()/logtwo),
        "joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "source_once_counts_outside_alpha":True,"no_extra_endpoint_support_count":True,
        "source_count_checksum":hashlib.sha256(str(sum(counts)).encode()).hexdigest()}
    if args.complete:
        if any(entry is None for entry in regions): raise ArithmeticError("incomplete actual-weight marginal cover")
        joined=join.directed_join(regions,guard,progress,weight=WEIGHT,points=POINTS)
        assert joined["report"]["ordered_profiles"]==4198400
        total=v.positive(source*joined["bound"])
        worst=max(range(1,WEIGHT+1),key=lambda c:joined["run_bounds"][c-1].upper())
        result.update({"all_profiles_at_weight_covered":True,
            "log2_union_bound_interval":str(total.log()/logtwo),
            "actual4096_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
            "actual4096_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
            "join":joined["report"],"worst_run":{"runs":worst,
                "log2_interval":str((source*joined["run_bounds"][worst-1]).log()/logtwo)}})
    guard()
    result.update({"elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "nonnegative_intersections":v.intersections,"generated_bytes":0,
        "pins":{"local":LOCAL_PIN,"regions":REGION_PIN,"integer":INTEGER_PIN,
            "join":JOIN_PIN,"anchor":ANCHOR_PIN,"boundary_tests":SCREEN_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
