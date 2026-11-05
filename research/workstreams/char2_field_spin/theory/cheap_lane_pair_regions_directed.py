"""Directed hybrid regional replay for all actual source supports of weight1024.

Use one coherent directed local binary-pair coefficient table. The low
three marker cells retain full degree exactly and singleton residue256;
four middle cells retain both residues256. One exact rectangle F<=512,
S<=158 replaces the eight high-full marker cells. Every retained target
is normalized by its OWN multinomial and marker monomial. Exact downset
truncation removes only aliases and is no weaker than the corresponding
residue bound for the same local coefficient majorant. No interpolation
or profile monotonicity is used. The eight tables cover f7..512.

The old-cone f0..6 boundary retains singleton residue512. Its outward
FLOAT exports are lifted as EXACT binary rationals, multiplied by
exp(integer scale - exact binary denominator lower bound). No printed
interval or decimal result is parsed into the proof.

The pair matrix H represents two updates; J=M*T has one update and the
output-only terminal. SIX doublings give H64 and H63J. The final directed
join uses exact indexed run counts times a maximum over every ordered
profile in that run class. The common source factor is applied ONCE.
This can certify one weight at35% distance, not the entire code.

Default cooperative budget300s;1200s requires an explicit option after a
measured pilot.512MiB/20GiB guards, stdout only, no coefficient files.
--test-center constructs only one selected table, including the high
rectangle at472; it makes no complete-weight claim and skips the join.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
import math
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
LOCAL_PIN="7284d1cd125d2cf1d1a3c1a38aa3dfd8b2ef7ceae4386ea2975406ea3e64f81f"
JOIN_PIN="54f2ff06a1476ccf84fa4db1b757f1d5f8e4e814a1f0f733dc73ccf18b2ddabb"
BOUNDARY_PIN="fdb3094bda1f319827c7f50d63cd68720c460cc66498f54575e81c5e2e518698"
HYBRID_PIN="985c3fe04f055811976d8060a4b5a93063a6cabd81dd6f4c285769f60e096727"
FINE_PIN="e4fcb2bbb4a73183d31f52ab2b86d5ee9e4906045dd045158fd75966d9aba1a2"
ANCHOR_PIN="319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
WITNESS_PIN="b2beb358bdd3d69b14718fec30b9b8f343ede5ea5134f2855f7b8321aa4dd801"
Z,ALPHA=F(304397809,625000000),F(119,1280)
CURVE_FULL,CURVE_SINGLE=F(1562175,990029),F(1140210,880327)
WEIGHT,MODULUS,MEMORY=1024,256,512*(1<<20)
SCHEDULE=((32,7,64,"mixed"),(96,65,128,"mixed"),(160,129,192,"mixed"),
          (224,193,256,"cyclic"),(288,257,320,"cyclic"),
          (352,321,384,"cyclic"),(416,385,432,"cyclic"),(472,433,512,"exact"))


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen source mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def marker(center):
    empty=4096-WEIGHT+center
    return CURVE_FULL*center/empty,CURVE_SINGLE*(WEIGHT-2*center)/empty


def lift_boundary(region,v):
    """Lift the outward exported tuple, not a rounded displayed interval."""
    regular,terminal,denominator,_,_=region
    if not math.isfinite(denominator): raise ArithmeticError("nonfinite boundary log denominator")
    denominator_exact=v.exact(F(*denominator.as_integer_ratio()))
    lifted=[]
    for entries,scale in (regular,terminal):
        if not math.isfinite(scale) or not scale.is_integer() or abs(scale)>2**53:
            raise ArithmeticError("boundary scale must be an exactly represented integer")
        factor=(v.exact(F(int(scale)))-denominator_exact).exp()
        row=[]
        for value in entries:
            if not math.isfinite(value) or value<=0:
                raise ArithmeticError("boundary export must be strictly positive")
            row.append(v.positive(v.exact(F(*value.as_integer_ratio()))*factor))
        lifted.append(tuple(row))
    if len(lifted[0])!=4 or len(lifted[1])!=2: raise ValueError("boundary R4/V2 required")
    return tuple(lifted)


def exact_controls(v,adapter,base,hybrid,fine,join,guard):
    coverage=[f for _,lo,hi,_ in SCHEDULE for f in range(lo,hi+1)]
    assert coverage==list(range(7,513))
    assert 1024-2*433==158 and SCHEDULE[-1][2]==512
    # Match every unchanged low/middle assignment to the original15-bank.
    original=(32,96,160,224,288,352,416,448,464,472,480,488,496,504,508)
    for center,lo,hi,kind in SCHEDULE:
        assert min(marker(center))>0
        if kind!="exact":
            assert all(min(original,key=lambda c:(abs(c-f),c))==center for f in range(lo,hi+1))
    toy=((tuple((.5,.25,.125,.0625)),2.0),((.375,.1875),-1.0),.125,F(1),F(1))
    lifted=lift_boundary(toy,v)
    for values,(entries,scale) in zip(lifted,toy[:2]):
        for actual,value in zip(values,entries):
            expected=v.exact(F(value))*(v.exact(F(scale)-F(toy[2]))).exp()
            assert actual.overlaps(expected)
    return {"cyclic_packed_toys":base.exact_tests(v,adapter),
        "mixed_packed_toys":hybrid.mixed_tests(v,adapter),
        "exact_rectangle_toys":fine.exact_tests(v,adapter),
        "exact_run_join_tests":join.exact_tests(guard),
        "target_coverage_and_original_marker_assignment":True,
        "boundary_binary_rational_lift_checks":6}


def one_table(v,degrees,item,base,hybrid,fine,guard,progress):
    """Return normalized flat(R4,V2) entries for this assigned cell only."""
    center,lo,hi,kind=item; u,w=marker(center)
    up=[v.exact(u)**k for k in range(65)]; wp=[v.exact(w)**k for k in range(65)]
    scaled={(f,s):tuple(up[f]*wp[s]*x for x in row) for (f,s),row in degrees.items()}
    plain=[sum((row[k] for row in scaled.values()),v.arb(0)) for k in range(6)]
    plain_h=[[plain[0],plain[1]],[plain[2],plain[3]]]; plain_j=[[plain[4]],[plain[5]]]
    evaluated=(v.matpow(plain_h,64),v.matmul(v.matpow(plain_h,63),plain_j))
    if kind=="mixed":
        full,single,stride=hi,MODULUS-1,2*hi+1
        rr,tt=hybrid.mixed_pack(v,scaled,full,MODULUS)
        multiply=lambda a,b: hybrid.mixed_product(v,a,b,full,MODULUS)
    elif kind=="cyclic":
        full,single,stride=MODULUS-1,MODULUS-1,2*MODULUS-1
        rr,tt=base.fold(v,scaled,MODULUS,MODULUS)
        multiply=lambda a,b: base.multiply(v,a,b,MODULUS,MODULUS)
    else:
        full,single,stride=512,158,1025
        rr,tt=fine.rectangle_pack(v,scaled,full,single)
        multiply=lambda a,b: fine.rectangle_product(v,a,b,full,single)
    progress("table_packed",center=center,method=kind,first_full=lo,last_full=hi,
        retained_coefficients_per_entry=single*stride+full+1,
        maximum_product_coefficients=2*full+stride*(2*single)+1)
    del scaled,up,wp,plain,plain_h,plain_j
    for step in range(6):
        guard(); tt=multiply(rr,tt); rr=multiply(rr,rr)
        progress("pair_doubling_complete",center=center,doubling=step+1,
            paired_factors=1<<(step+1))
    dc_checks=0
    for table,ordinary in zip((rr,tt),evaluated):
        for i,row in enumerate(table):
            for j,poly in enumerate(row):
                guard()
                dc=sum((poly[f+stride*s] for s in range(single+1) for f in range(full+1)),v.arb(0))
                if kind=="cyclic": assert dc.overlaps(ordinary[i][j])
                else: assert dc<=ordinary[i][j] or dc.overlaps(ordinary[i][j])
                dc_checks+=1
    result={}
    for f in range(lo,hi+1):
        guard(); s=WEIGHT-2*f
        ff=f%MODULUS if kind=="cyclic" else f
        ss=s if kind=="exact" else s%MODULUS
        assert ff<=full and ss<=single
        index=ff+stride*ss
        denominator=v.arb(comb(4096,f)*comb(4096-f,s))*v.exact(u)**f*v.exact(w)**s
        assert denominator>0
        result[f]=(tuple(v.positive(poly[index])/denominator for row in rr for poly in row),
                   tuple(v.positive(poly[index])/denominator for row in tt for poly in row))
    selected=lo+(hi-lo)//2
    report={"center":center,"method":kind,"first_full":lo,"last_full":hi,
        "full_exact_limit":full if kind!="cyclic" else None,
        "full_residue":MODULUS if kind=="cyclic" else None,
        "singleton_exact_limit":single if kind=="exact" else None,
        "singleton_residue":MODULUS if kind!="exact" else None,
        "packing_stride":stride,"u":str(u),"v":str(w),"dc_checks":dc_checks,
        "pair_doublings":6,"regular_packets":128,"terminal_updates":127,
        "sample_full":selected,"sample_entry_log2_intervals":[str(x.log()/v.arb(2).log())
            for values in result[selected] for x in values]}
    # Do not let a loop alias retain one large polynomial across cells.
    del rr,tt,table,ordinary,poly,row,evaluated,multiply
    gc.collect(); guard()
    return result,report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision",type=int,choices=(256,384),default=256)
    parser.add_argument("--seconds",type=int,choices=(300,1200),default=300)
    parser.add_argument("--test-center",type=int,choices=tuple(x[0] for x in SCHEDULE))
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    anchor=load("cheap_lane_cone_cauchy_anchors.py",ANCHOR_PIN,"pair_regions_anchor")
    local=load("cheap_lane_pair_local_directed.py",LOCAL_PIN,"pair_regions_local")
    join=load("cheap_lane_directed_run_mass.py",JOIN_PIN,"pair_regions_join")
    boundary=load("cheap_lane_old_boundary_regions.py",BOUNDARY_PIN,"pair_regions_boundary")
    hybrid=load("cheap_lane_cauchy_hybrid_directed.py",HYBRID_PIN,"pair_regions_hybrid")
    fine=load("cheap_lane_cauchy_fine.py",FINE_PIN,"pair_regions_fine")
    # Witness pin is checked but never imported: no FLOAT proof path.
    witness=Path(__file__).with_name("cheap_lane_pair_residue_control.py")
    if hashlib.sha256(witness.read_bytes()).hexdigest()!=WITNESS_PIN:
        raise RuntimeError("fixed FLOAT witness source changed")
    base=hybrid.load_directed().load_base(); adapter=anchor.load_adapter()
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("directed paired-region deadline; no retry or incomplete claim")
        if peak()>MEMORY: raise MemoryError("512MiB directed paired-region budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"elapsed_seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]
    del built["rational_coefficients"]
    tests=exact_controls(v,adapter,base,hybrid,fine,join,guard)
    progress("local_pair_coefficients_ready",rational_coefficient_sha256=built["report"]["rational_coefficient_sha256"],
        precision=args.precision,test_center=args.test_center)
    regions=[None]*513; reports=[]
    schedule=SCHEDULE if args.test_center is None else tuple(x for x in SCHEDULE if x[0]==args.test_center)
    for item in schedule:
        partial,report=one_table(v,built["coefficients"],item,base,hybrid,fine,guard,progress)
        for f,region in partial.items():
            assert regions[f] is None; regions[f]=region
        reports.append(report); del partial
        progress("regional_table_complete",**report)
    result={"scope":"directed paired regional table pilot only" if args.test_center is not None
        else "all actual weight1024 source supports, matching profiles and state paths at35% distance",
        "directed":True,"all_profiles_at_weight_covered":False,"whole_code_certificate":False,
        "parameters":{"weight":WEIGHT,"z":str(Z),"alpha":str(ALPHA),"precision_bits":args.precision,
            "paired_modulus":256,"boundary_singleton_modulus":512,"test_center":args.test_center,
            "budget_seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "tables":reports,"local_coefficients":built["report"],"tests":tests,
        "hybrid_refinement":"positive exact full/downset truncation; no extra alpha; each target owns its normalizer"}
    if args.test_center is None:
        patch,boundary_report=boundary.old_boundary_regions(Z,ALPHA,anchor,None,guard,deadline,
            precision=args.precision,check_singleton=True)
        assert sorted(patch)==list(range(7))
        assert boundary_report["z"]==str(Z) and boundary_report["alpha"]==str(ALPHA)
        assert not boundary_report["source_factor_applied"] and not boundary_report["support_count_applied"]
        for f,region in patch.items():
            assert regions[f] is None; regions[f]=lift_boundary(region,v)
        assert all(x is not None for x in regions)
        progress("directed_boundary_lift_complete",full_indices=list(range(7)),
            exact_binary_input_values=True,source_factor_applied=False)
        h,c_h,eta=anchor.source_factor(WEIGHT)
        source=(v.exact(eta)*v.exact(core.BETA)**WEIGHT*v.exact(Z)**(-core.LIMIT))**v.exact(ALPHA)
        joined=join.directed_join(regions,guard,progress)
        total=v.positive(source*joined["bound"]); assert total>0
        logtwo=v.arb(2).log()
        worst=max(range(1,1025),key=lambda c:joined["run_bounds"][c-1].upper())
        result.update({"all_profiles_at_weight_covered":True,"log2_union_bound_interval":str(total.log()/logtwo),
            "weight1024_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
            "weight1024_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
            "boundary":boundary_report,"boundary_lift":"exact binary upper entries times exp(integer scale minus exact binary lower log denominator)",
            "join":joined["report"],"source_log2_interval":str(source.log()/logtwo),
            "joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
            "worst_run":{"runs":worst,"log2_interval":str((source*joined["run_bounds"][worst-1]).log()/logtwo)}})
    guard()
    result.update({"elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "nonnegative_intersections":v.intersections,"local_pin":LOCAL_PIN,"join_pin":JOIN_PIN,
        "boundary_pin":BOUNDARY_PIN,"hybrid_pin":HYBRID_PIN,"fine_pin":FINE_PIN,
        "base_pin":fine.BASE_PIN,"anchor_pin":ANCHOR_PIN,"fixed_witness_pin":WITNESS_PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
