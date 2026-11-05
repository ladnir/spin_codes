"""Directed old-cone boundary regions for actual r1024, full degree0..6.

Retain full degree exactly through6 and singleton degree modulo512.
This positive quotient removes full-degree aliases but keeps singleton
aliases. It does not claim the selected coefficients are exact. The
full marker is1; the singleton marker is the declared endpoint witness.
Seven positive mixed-axis doublings produce P128 and P127*T.

Use the OLD frozen cone, powered before route averaging. Its regional
entries majorize the later refined/corrected EXACT paired regional object,
not necessarily that object's different cyclic-residue upper bound.
Each target divides by its own C4096f*C(4096-f,1024-2f)*v^(1024-2f).
No source factor or support count enters the reusable boundary function.
Its seven returned tuples have the broad.marginals interface(P,T,logden,u,v).
Arb upper endpoints are used to construct scaled FLOAT entries; any
subsequent complete FLOAT join remains a diagnostic, not a certificate.

The optional cheap f0 control computes the exact univariate coefficient
and checks that the cyclic boundary encloses it from above. Frozen
mixed-axis and singleton toys check truncation, carries and terminal order.
180s/256MiB/20GiB, stdout only; no generated data or frozen edits.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
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
HYBRID_PIN="985c3fe04f055811976d8060a4b5a93063a6cabd81dd6f4c285769f60e096727"
SINGLETON_PIN="88f5156ba031f348b23d2e824fd37886822ac1b63874634063272014d2d5577a"
ANCHOR_PIN="319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
R,FULL,MODULUS,STRIDE=1024,6,512,13
SCALE=F(1140210,880327)/3
SECONDS,MEMORY=180,256*(1<<20)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen boundary dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def local_records(guard):
    """Actual paired halves, with integer fair-orientation numerators."""
    half=[]
    for f in range(FULL+1):
        for s in range(17-f):
            for left in range(s+1):
                half.append((f,s,f+left,f+s-left,
                    comb(16,f)*comb(16-f,s)*comb(s,left)*(1<<(16-s))))
    records=defaultdict(int)
    for i,(f,s,a,b,mass) in enumerate(half):
        guard()
        for ff,ss,c,d,other in half:
            if f+ff<=FULL:
                records[f+ff,s+ss,tuple(sorted((a,b,c,d)))]+=mass*other
    totals=defaultdict(int)
    for (f,s,_),mass in records.items(): totals[f,s]+=mass
    for f in range(FULL+1):
        for s in range(33-f):
            assert totals[f,s]==comb(32,f)*comb(32-f,s)*(1<<32)
    return dict(records)


def upper_float(value,v):
    """Round an Arb upper endpoint outward to a finite binary FLOAT."""
    if not value.is_finite(): raise ArithmeticError("nonfinite boundary enclosure")
    endpoint=value.upper().fmpq(); result=float(endpoint)
    if not math.isfinite(result): raise ArithmeticError("boundary FLOAT overflow")
    while v.fmpq(*result.as_integer_ratio())<endpoint:
        result=math.nextafter(result,math.inf)
    return result


def normalized_upper(items,v):
    """Positive scaled tuple; no screen re-normalization can round it down."""
    if not all(x>0 for x in items):
        raise ArithmeticError("boundary coefficient not provably positive")
    logs=[x.log() for x in items]
    scale=math.ceil(max(upper_float(x,v) for x in logs))
    exact_scale=v.exact(F(scale))
    entries=tuple(upper_float((x-exact_scale).exp(),v) for x in logs)
    if not all(0<x<=1+1e-12 for x in entries):
        raise ArithmeticError("boundary normalized entry underflow or scale error")
    return entries,float(scale)


def old_boundary_regions(z,alpha,anchor,screen,guard,deadline,
                         precision=256,single_scale=SCALE,check_singleton=True):
    """Return({f: broad-compatible region}, directed provenance report).

    The caller owns resource limits. No source/count factor is applied.
    screen is accepted for interface compatibility; no FLOAT screen
    arithmetic is used to construct these regional upper enclosures.
    """
    if not 0<z<1 or not 0<alpha<=1 or single_scale<=0:
        raise ValueError("invalid boundary witnesses")
    if precision not in (256,384): raise ValueError("precision must be256 or384")
    hybrid=load("cheap_lane_cauchy_hybrid_directed.py",HYBRID_PIN,"boundary_hybrid")
    singleton=load("cheap_lane_singleton_endpoint_directed.py",SINGLETON_PIN,"boundary_singleton")
    adapter=anchor.load_adapter(); core=adapter.load_verifier()
    core.Z,core.ALPHA=z,alpha
    from flint import ctx
    ctx.prec=precision
    class BoundaryVerifier(core.Verifier):
        def guard(self): guard()
    v=BoundaryVerifier(2,deadline)
    tests={"mixed_exact_toy":hybrid.mixed_tests(v,adapter),
           "singleton_exact_toy":singleton.exact_tests(v,adapter)}
    records=local_records(guard)
    shapes=sorted({occ for _,_,occ in records})
    caps=v.local_caps(shapes); exponent=v.exact(alpha)
    powered={occ:tuple(v.positive(x)**exponent for x in row) for occ,row in caps.items()}
    degrees={(f,s):[v.arb(0) for _ in range(6)]
             for f in range(FULL+1) for s in range(33-f)}
    for (f,s,occ),numerator in records.items():
        guard(); weight=v.exact(F(numerator,1<<32)*single_scale**s)
        for k,value in enumerate(powered[occ]): degrees[f,s][k]+=weight*value
    del caps,powered,records
    plain=[sum((row[k] for row in degrees.values()),v.arb(0)) for k in range(6)]
    plain_p=[[plain[0],plain[1]],[plain[2],plain[3]]]
    plain_t=[[plain[4]],[plain[5]]]
    ordinary=(v.matpow(plain_p,128),v.matmul(v.matpow(plain_p,127),plain_t))
    rr,tt=hybrid.mixed_pack(v,degrees,FULL,MODULUS)
    for _ in range(7):
        tt=hybrid.mixed_product(v,rr,tt,FULL,MODULUS)
        rr=hybrid.mixed_product(v,rr,rr,FULL,MODULUS)
    dc_checks=0
    for table,evaluated in zip((rr,tt),ordinary):
        for i,row in enumerate(table):
            for j,item in enumerate(row):
                dc=sum((item[f+STRIDE*s] for f in range(FULL+1) for s in range(MODULUS)),v.arb(0))
                assert dc<=evaluated[i][j] or dc.overlaps(evaluated[i][j])
                dc_checks+=1
    regions,intervals={},{}
    for f in range(FULL+1):
        guard(); s=R-2*f; index=f+STRIDE*(s%MODULUS)
        numerator_r=[v.positive(item[index]) for row in rr for item in row]
        numerator_t=[v.positive(item[index]) for row in tt for item in row]
        denominator=v.exact(F(comb(4096,f)*comb(4096-f,s))*single_scale**s)
        assert denominator>0
        # The downstream join subtracts this denominator. Use a lower
        # enclosure here, and upper enclosures for every numerator.
        den_lower=denominator.log().lower().fmpq()
        den_float=float(den_lower)
        while v.fmpq(*den_float.as_integer_ratio())>den_lower:
            den_float=math.nextafter(den_float,-math.inf)
        regions[f]=(normalized_upper(numerator_r,v),normalized_upper(numerator_t,v),
                    den_float,F(1),single_scale)
        intervals[f]={"full":f,"single":s,
            "regular_log_intervals":[str((x/denominator).log()) for x in numerator_r],
            "terminal_log_intervals":[str((x/denominator).log()) for x in numerator_t],
            "denominator_log_interval":str(denominator.log())}
    exact_control=None
    if check_singleton:
        polys=[v.poly([degrees[0,s][k] for s in range(33)]) for k in range(6)]
        ar,at=[[polys[0],polys[1]],[polys[2],polys[3]]],[[polys[4]],[polys[5]]]
        for _ in range(7):
            at=singleton.multiply(v,ar,at,R)
            ar=singleton.multiply(v,ar,ar,R)
        control_checks=0
        denominator=v.exact(F(comb(4096,R))*single_scale**R)
        exact_r,exact_t=[],[]
        for a,b,destination in ((ar,rr,exact_r),(at,tt,exact_t)):
            for i,row in enumerate(a):
                rr_out=[]
                for j,item in enumerate(row):
                    exact=v.positive(item[R]); alias=v.positive(b[i][j][0])
                    assert exact<=alias or exact.overlaps(alias)
                    rr_out.append(exact/denominator); control_checks+=1
                destination.append(rr_out)
        chain=v.positive(v.matmul(v.matpow(exact_r,255),exact_t)[0][0])
        exact_control={"coefficient_checks":control_checks,
            "exact_f0_conditional_log2_interval":str(chain.log()/v.arb(2).log()),
            "scope":"exact f0 old-cone conditional chain only; source/count omitted"}
    guard()
    report={"scope":"old-cone regional upper bounds, actual r1024 full degrees0..6",
        "directed_regional_arithmetic":True,"subsequent_FLOAT_join_is_not_a_certificate":True,
        "z":str(z),"alpha":str(alpha),"full_limit":FULL,"singleton_modulus":MODULUS,
        "full_scale":"1","singleton_scale":str(single_scale),"stride":STRIDE,
        "retained_polynomial_slots":6650,"maximum_product_slots":13299,
        "regional_packets":128,"terminal_updates":127,"tests":tests,
        "local_degree_normalizers":210,"truncated_DC_checks":dc_checks,
        "exact_singleton_control":exact_control,"intervals":intervals,
        "nonnegative_intersections":v.intersections,"source_factor_applied":False,
        "support_count_applied":False,"hybrid_pin":HYBRID_PIN,"singleton_pin":SINGLETON_PIN,
        "core_pin":anchor.CORE_PIN,"precision_bits":precision,"generated_bytes":0}
    return regions,report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--v",type=F,default=SCALE)
    parser.add_argument("--precision",type=int,choices=(256,384),default=256)
    args=parser.parse_args()
    if not 0<args.z<1 or not 0<args.alpha<=1 or args.v<=0: parser.error("invalid witnesses")
    if max(n.bit_length() for x in (args.z,args.alpha,args.v) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    anchor=load("cheap_lane_cone_cauchy_anchors.py",ANCHOR_PIN,"boundary_anchor")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("180-second boundary deadline; no retry")
        if peak()>MEMORY: raise MemoryError("256MiB boundary budget exceeded")
    regions,report=old_boundary_regions(args.z,args.alpha,anchor,None,guard,deadline,
                                        args.precision,args.v)
    assert sorted(regions)==list(range(7))
    report.update({"FLOAT_export_regions":{f:[x[0],x[1],x[2],str(x[3]),str(x[4])] for f,x in regions.items()},
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    print(json.dumps(report,indent=2),flush=True)


if __name__=="__main__":
    main()
