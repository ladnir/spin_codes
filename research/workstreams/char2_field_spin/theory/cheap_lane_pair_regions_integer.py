"""Directed r1024 regional cover using exact integer products and upward dyadics.

Pin the frozen local-cap table, eight-table hybrid schedule, source factor,
and indexed-run join. Replace only the expensive regional Arb polynomial
arithmetic. Nonnegative integer polynomial products are exact. Positive
folding/downset truncation is followed by upward dyadic quantization at
128bits, with a separate common exponent for each matrix or vector. This
is not a claim of128 relative bits for each tiny coefficient. Every
exported coefficient is therefore a genuine upper bound; no assumed FFT
error model, discarded small coefficient, or new probability law enters.

Six doublings give H64 and H63J with J=M*T. Each target removes its own
exact multinomial/marker normalizer. The unchanged f0..6 old-cone boundary
is lifted from its exact binary rational upper exports. A directed run-mass
join covers all263168 ordered profiles, source applied once outside.

The cyclic DC control is now an UPPER comparison, not an equality check.
For exact/truncated axes, rounding and truncation have opposite effects;
no comparison with the untruncated DC is asserted. Exact rational toys
check the integer arithmetic and upward-rounding invariant instead.

--test-center selects only one table and makes no complete-weight claim.
300s default/1200s explicit;512MiB/20GiB guards; no generated files.
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
DRIVER_PIN="d3a6e858437e19bede76e215a9d9b31f26a49b9b601ecd60705f1a2b99f0f213"
INTEGER_PIN="9dc583ae54caed0a78f71e409cc37aa398713b419f65538370f10c3f876f40bf"
MEMORY=512*(1<<20)
BITS=128


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen source mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def one_table(v,degrees,item,driver,integer,guard,progress):
    center,lo,hi,kind=item; u,w=driver.marker(center)
    up=[u**k for k in range(65)]; wp=[w**k for k in range(65)]
    scaled={(f,s):tuple(up[f]*wp[s]*x for x in row) for (f,s),row in degrees.items()}
    if kind=="mixed":
        full,single,cf,cs=hi,255,False,True
    elif kind=="cyclic":
        full,single,cf,cs=255,255,True,True
    else:
        full,single,cf,cs=512,158,False,False
    arithmetic=integer.IntegerArithmetic(full,single,cyclic_full=cf,cyclic_single=cs,
        bits=BITS,guard=guard)
    # The reference DC uses the original coherent exact rational pair
    # coefficients, BEFORE the integer engine rounds or truncates them.
    plain=[v.exact(sum((row[k] for row in scaled.values()),F(0))) for k in range(6)]
    plain_h=[[plain[0],plain[1]],[plain[2],plain[3]]]; plain_j=[[plain[4]],[plain[5]]]
    evaluated=(v.matpow(plain_h,64),v.matmul(v.matpow(plain_h,63),plain_j))
    rr,tt=arithmetic.pack(scaled)
    progress("integer_table_packed",center=center,method=kind,first_full=lo,last_full=hi,
        matrix_exponent=rr.exponent,terminal_exponent=tt.exponent,arithmetic=arithmetic.report())
    del scaled,up,wp,plain,plain_h,plain_j
    for step in range(6):
        guard(); tt=arithmetic.multiply(rr,tt); rr=arithmetic.multiply(rr,rr)
        progress("integer_pair_doubling_complete",center=center,doubling=step+1,
            paired_factors=1<<(step+1),matrix_exponent=rr.exponent,terminal_exponent=tt.exponent)
    dc_checks=0; dc_ratios=[]
    for matrix,ordinary in zip((rr,tt),evaluated):
        upper=arithmetic.mass(matrix)
        for i,row in enumerate(upper):
            for j,value in enumerate(row):
                guard(); dc=v.exact(value)
                if kind=="cyclic":
                    assert dc>=ordinary[i][j] or dc.overlaps(ordinary[i][j])
                    dc_checks+=1
                if dc>0 and ordinary[i][j]>0:
                    dc_ratios.append(str((dc/ordinary[i][j]).log()/v.arb(2).log()))
    result={}
    for f in range(lo,hi+1):
        guard(); s=driver.WEIGHT-2*f
        ff=f%driver.MODULUS if cf else f
        ss=s%driver.MODULUS if cs else s
        assert ff<=full and ss<=single
        denominator=F(comb(4096,f)*comb(4096-f,s))*u**f*w**s
        assert denominator>0
        result[f]=tuple(tuple(v.exact(x/denominator) for row in arithmetic.extract(matrix,ff,ss) for x in row)
                        for matrix in (rr,tt))
    selected=lo+(hi-lo)//2
    report={"center":center,"method":kind,"first_full":lo,"last_full":hi,
        "full_exact_limit":full if not cf else None,"full_residue":full+1 if cf else None,
        "singleton_exact_limit":single if not cs else None,"singleton_residue":single+1 if cs else None,
        "u":str(u),"v":str(w),"arithmetic":arithmetic.report(),"pair_doublings":6,
        "regular_packets":128,"terminal_updates":127,"cyclic_DC_upper_checks":dc_checks,
        "DC_vs_untruncated_log2_ratios":dc_ratios,
        "DC_ratio_scope":"diagnostic only for truncated modes; no <= or equality assertion",
        "sample_full":selected,"sample_entry_log2_intervals":[str(x.log()/v.arb(2).log())
            for values in result[selected] for x in values]}
    del rr,tt,matrix,ordinary,upper,row,evaluated,arithmetic
    gc.collect(); guard()
    return result,report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision",type=int,choices=(256,384),default=256)
    parser.add_argument("--seconds",type=int,choices=(300,1200),default=300)
    parser.add_argument("--test-center",type=int,choices=(32,96,160,224,288,352,416,472))
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    driver=load("cheap_lane_pair_regions_directed.py",DRIVER_PIN,"integer_regions_frozen_driver")
    integer=load("cheap_lane_integer_regions.py",INTEGER_PIN,"integer_regions_arithmetic")
    anchor=driver.load("cheap_lane_cone_cauchy_anchors.py",driver.ANCHOR_PIN,"integer_regions_anchor")
    local=driver.load("cheap_lane_pair_local_directed.py",driver.LOCAL_PIN,"integer_regions_local")
    join=driver.load("cheap_lane_directed_run_mass.py",driver.JOIN_PIN,"integer_regions_join")
    boundary=driver.load("cheap_lane_old_boundary_regions.py",driver.BOUNDARY_PIN,"integer_regions_boundary")
    hybrid=driver.load("cheap_lane_cauchy_hybrid_directed.py",driver.HYBRID_PIN,"integer_regions_hybrid")
    fine=driver.load("cheap_lane_cauchy_fine.py",driver.FINE_PIN,"integer_regions_fine")
    base=hybrid.load_directed().load_base(); adapter=anchor.load_adapter()
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("integer regional deadline; no retry or incomplete claim")
        if peak()>MEMORY: raise MemoryError("512MiB integer regional budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"elapsed_seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    built=local.build_local_pairs(driver.Z,driver.ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]
    del built["coefficients"]
    tests=driver.exact_controls(v,adapter,base,hybrid,fine,join,guard)
    tests["integer_arithmetic"]=integer.exact_tests(guard)
    progress("local_pair_coefficients_ready",rational_coefficient_sha256=built["report"]["rational_coefficient_sha256"],
        precision=args.precision,dyadic_bits=BITS,test_center=args.test_center)
    regions=[None]*513; reports=[]
    schedule=driver.SCHEDULE if args.test_center is None else tuple(x for x in driver.SCHEDULE if x[0]==args.test_center)
    for item in schedule:
        partial,report=one_table(v,built["rational_coefficients"],item,driver,integer,guard,progress)
        for f,region in partial.items():
            assert regions[f] is None; regions[f]=region
        reports.append(report); del partial
        progress("integer_regional_table_complete",**report)
    result={"scope":"directed integer regional table pilot only" if args.test_center is not None
        else "all actual weight1024 source supports, matching profiles and state paths at35% distance",
        "directed":True,"all_profiles_at_weight_covered":False,"whole_code_certificate":False,
        "arithmetic":"exact nonnegative fmpz_poly products with coefficientwise upward128-bit dyadic quantization",
        "parameters":{"weight":driver.WEIGHT,"z":str(driver.Z),"alpha":str(driver.ALPHA),
            "local_and_final_Arb_precision_bits":args.precision,"regional_upward_dyadic_bits":BITS,
            "paired_modulus":256,"boundary_singleton_modulus":512,"test_center":args.test_center,
            "budget_seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "tables":reports,"local_coefficients":built["report"],"tests":tests,
        "unchanged_frozen_witness":{"schedule":driver.SCHEDULE,"local_caps":True,
            "source_and_counts":True,"terminal":True,"no_new_probability_distribution":True}}
    if args.test_center is None:
        patch,boundary_report=boundary.old_boundary_regions(driver.Z,driver.ALPHA,anchor,None,guard,deadline,
            precision=args.precision,check_singleton=True)
        assert sorted(patch)==list(range(7))
        assert boundary_report["z"]==str(driver.Z) and boundary_report["alpha"]==str(driver.ALPHA)
        assert not boundary_report["source_factor_applied"] and not boundary_report["support_count_applied"]
        for f,region in patch.items():
            assert regions[f] is None; regions[f]=driver.lift_boundary(region,v)
        assert all(x is not None for x in regions)
        progress("directed_boundary_lift_complete",full_indices=list(range(7)),exact_binary_input_values=True)
        h,c_h,eta=anchor.source_factor(driver.WEIGHT)
        source=(v.exact(eta)*v.exact(core.BETA)**driver.WEIGHT*v.exact(driver.Z)**(-core.LIMIT))**v.exact(driver.ALPHA)
        joined=join.directed_join(regions,guard,progress)
        total=v.positive(source*joined["bound"]); assert total>0
        logtwo=v.arb(2).log()
        worst=max(range(1,1025),key=lambda c:joined["run_bounds"][c-1].upper())
        result.update({"all_profiles_at_weight_covered":True,"log2_union_bound_interval":str(total.log()/logtwo),
            "weight1024_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
            "weight1024_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
            "boundary":boundary_report,"join":joined["report"],
            "source_log2_interval":str(source.log()/logtwo),
            "joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
            "worst_run":{"runs":worst,"log2_interval":str((source*joined["run_bounds"][worst-1]).log()/logtwo)}})
    guard()
    result.update({"elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "nonnegative_intersections":v.intersections,"frozen_driver_pin":DRIVER_PIN,
        "integer_arithmetic_pin":INTEGER_PIN,"local_pin":driver.LOCAL_PIN,"join_pin":driver.JOIN_PIN,
        "boundary_pin":driver.BOUNDARY_PIN,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
