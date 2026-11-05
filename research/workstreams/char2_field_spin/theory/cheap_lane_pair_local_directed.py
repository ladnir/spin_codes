"""Directed local coefficients for the refined binary-pair cone.

build_local_pairs returns six coefficient channels at each degree(f,s),
f+s<=64. The first four channels form H=M*M; the last two form J=M*T.
The caller supplies regional marker powers, normalization, source factors,
and source counts. None is applied here. A region has64 regular pairs,
or63 regular pairs followed by J; this module does not form that chain.

First freeze every genuine ell1/2/3 raw cap to its exact rational upper
endpoint. All powers and corrections use THIS coherent cap table. For
each occupancy and fixed rho, form the five powered minima before route
averaging. Exact dyadic route weights and positive packed convolutions
then give a lower bound on the integrated path minimum. Subtract only
its rational lower endpoint, times a rational lower bound on2-2^alpha,
from a rational upper old-pair coefficient. Intersect with an independent
upper degree-Jensen coefficient. The resulting rational coefficients are
upper bounds, not intervals claimed to enclose the optimal pair bound.

Packing f+65*s is carry-free through the two-packet product: f<=64.
There are12 old-pair,12 raw-pair, and6 products per rho (default9).
No FLOAT value enters an exported coefficient. The pinned cap builder
has a FLOAT representability guard, used only to abort unsuitable inputs.
No numerical files, encoder edits, source counts, or probability claim.
Standalone budgets:180--300s,256MiB,20GiB free; stdout only.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
import json
from math import comb, isqrt
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
COMPOSITION_PIN="57c3abe184ce6ee643779c48106b79f1d39c85bfb4175e28404ec86b88bc7169"
CORRECTION_PIN="804eeb5b9e7df90cba0ffa57651b0e411784e5bb775e6e3d62b532e994728cf1"
ANCHOR_PIN="319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
STRIDE,PACKET_SIZE,PAIR_SIZE=65,32,64
PACKET_LENGTH=STRIDE*PACKET_SIZE+1
CHANNELS=("Q","A","E","L","T0","U")
PATHS=((0,0,1,2),(0,1,1,3),(2,0,3,2),(2,1,3,3),(0,4,1,5),(2,4,3,5))
MEMORY=256*(1<<20)


def load(name,digest,module_name):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
        raise RuntimeError("frozen helper hash mismatch: "+name)
    spec=importlib.util.spec_from_file_location(module_name,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def endpoint(value,upper):
    """Exact rational endpoint; no decimal or binary64 conversion."""
    if not value.is_finite(): raise ArithmeticError("nonfinite directed coefficient")
    return F(str((value.upper() if upper else value.lower()).fmpq()))


def enclose_point(value,point):
    assert endpoint(value,False)<=point<=endpoint(value,True)


def route_geometry(shapes,guard):
    """Grouped exact records (shape index,packed degree,integer/2^32)."""
    indices={occ:i for i,occ in enumerate(shapes)}
    half=[]
    for f in range(17):
        for s in range(17-f):
            for left in range(s+1):
                a,b=f+left,f+s-left
                if a<=b:
                    numerator=comb(16,f)*comb(16-f,s)*comb(s,left)
                    numerator*=(1 if a==b else 2)*(1<<(16-s))
                    half.append((f,s,a,b,numerator))
    records=defaultdict(int); raw_records=0
    for i,(f,s,a,b,w) in enumerate(half):
        guard()
        for j in range(i,len(half)):
            ff,ss,c,d,ww=half[j]
            records[indices[tuple(sorted((a,b,c,d)))],f+ff+STRIDE*(s+ss)]+=w*ww*(1 if i==j else 2)
            raw_records+=1
    by_shape=[[] for _ in shapes]; masses=defaultdict(int)
    for (oi,index),numerator in sorted(records.items()):
        by_shape[oi].append((index,numerator)); masses[index]+=numerator
    checks=0
    for f in range(33):
        for s in range(33-f):
            assert masses[f+STRIDE*s]==comb(32,f)*comb(32-f,s)*(1<<32)
            checks+=1
    assert len(masses)==561
    return by_shape,{"raw_route_records":raw_records,"grouped_route_records":len(records),
                     "exact_packet_normalizers":checks,"route_denominator":1<<32}


def route_polys(values,geometry,v,guard):
    """Average supplied per-occupancy values; never take a minimum here."""
    n=len(values[0]); coeff=[[v.arb(0) for _ in range(PACKET_LENGTH)] for _ in range(n)]
    for row,records in zip(values,geometry):
        guard()
        for index,numerator in records:
            for k,value in enumerate(row):
                coeff[k][index]+=value*numerator
    denominator=1<<32
    return tuple(v.poly([x/denominator for x in channel]) for channel in coeff)


def pair_products(packet,guard):
    out=[]
    for a,b,c,d in PATHS:
        guard(); out.append(packet[a]*packet[b]+packet[c]*packet[d])
    return tuple(out)


def cross_products(minima,guard):
    out=[]
    for a,b in ((0,2),(0,3),(1,2),(1,3),(0,4),(1,4)):
        guard(); out.append(minima[a]*minima[b])
    return tuple(out)


def occupancy_minima(row,rho):
    q,a,e,l,t,u=row
    return ((rho*q).min(a),(rho*e).min(l),(q/rho).min(e),
            (a/rho).min(l),(t/rho).min(u))


def sqrt_bounds(x):
    scale=1<<100; k=isqrt(x.numerator*scale*scale//x.denominator)
    lo=F(k,scale)
    return lo,lo if lo*lo==x else F(k+1,scale)


def exact_tests(v,guard):
    """Rational route toys, packed carries, terminal order and alpha1."""
    checks=0
    for a,b,c,d,rho in product((F(0),F(1,3),F(3,2)),repeat=5):
        if not rho: continue
        lower=min(rho*a,c)*min(b/rho,d)
        assert 0<=lower<=min(a*b,c*d)<=(a*b+c*d)/2
        checks+=1
    atoms=(((0,0),F(2),tuple(map(F,(1,2,3,5,7,11)))),
           ((1,0),F(3),tuple(map(F,(2,3,5,7,11,13)))),
           ((0,1),F(1),tuple(map(F,(3,5,7,11,13,17)))))
    coeff=[[v.arb(0) for _ in range(STRIDE+2)] for _ in CHANNELS]
    raw_coeff=[[v.arb(0) for _ in range(STRIDE+2)] for _ in CHANNELS]
    expected=defaultdict(lambda:[F(0)]*6); raw=defaultdict(lambda:[F(0)]*6)
    target_lo=defaultdict(lambda:[F(0)]*6); route_mass=defaultdict(F)
    for (f,s),w,row in atoms:
        for j,x in enumerate(row):
            coeff[j][f+STRIDE*s]+=v.exact(w*x)
            raw_coeff[j][f+STRIDE*s]+=v.exact(w*x*x)
    packed=pair_products(tuple(v.poly(x) for x in coeff),guard)
    raw_packed=pair_products(tuple(v.poly(x) for x in raw_coeff),guard)
    for (da,wa,a),(db,wb,b) in product(atoms,repeat=2):
        index=da[0]+db[0]+STRIDE*(da[1]+db[1]); w=wa*wb
        route_mass[index]+=w
        for j,(ai,bi,ci,di) in enumerate(PATHS):
            x,y=a[ai]*b[bi],a[ci]*b[di]
            expected[index][j]+=w*(x+y); raw[index][j]+=w*(x*x+y*y)
            target_lo[index][j]+=w*sqrt_bounds(x*x+y*y)[0]
    best=defaultdict(lambda:[F(0)]*6)
    for rho in (F(1,16),F(1,2),F(1),F(4),F(16)):
        m=[[v.arb(0) for _ in range(STRIDE+2)] for _ in range(5)]
        ex=defaultdict(lambda:[F(0)]*6)
        for (f,s),w,row in atoms:
            vals=(min(rho*row[0],row[1]),min(rho*row[2],row[3]),
                  min(row[0]/rho,row[2]),min(row[1]/rho,row[3]),min(row[4]/rho,row[5]))
            for j,x in enumerate(vals): m[j][f+STRIDE*s]+=v.exact(w*x)
        actual=cross_products(tuple(v.poly(x) for x in m),guard)
        for (da,wa,a),(db,wb,b) in product(atoms,repeat=2):
            index=da[0]+db[0]+STRIDE*(da[1]+db[1])
            for j,(ai,bi,ci,di) in enumerate(PATHS):
                ex[index][j]+=wa*wb*min(rho*a[ai],a[ci])*min(b[bi]/rho,b[di])
        for index,row in ex.items():
            for j,x in enumerate(row):
                enclose_point(actual[j][index],x); best[index][j]=max(best[index][j],x)
                checks+=1
    _,sqrt2hi=sqrt_bounds(F(2))
    for index,row in expected.items():
        for j,x in enumerate(row):
            enclose_point(packed[j][index],x); enclose_point(raw_packed[j][index],raw[index][j])
            assert 2*best[index][j]<=x
            correction_upper=x-(2-sqrt2hi)*best[index][j]
            assert target_lo[index][j]<=correction_upper
            # Exact rational root enclosures also cover equality cases.
            assert target_lo[index][j]**2<=route_mass[index]*raw[index][j]
            # alpha1: c=0, both old and raw/Jensen recover M*M and M*T.
            assert x-F(0)*best[index][j]==x
            checks+=5
    # Full-axis extreme32+32 must not carry into singleton degree1.
    p=[v.arb(0)]*(32+STRIDE*32+1)
    p[32]=v.exact(F(2,7)); p[STRIDE*32]=v.exact(F(3,11))
    square=v.poly(p)**2
    for index,value in ((64,F(4,49)),(32+STRIDE*32,F(12,77)),(STRIDE*64,F(9,121))):
        enclose_point(square[index],value); checks+=1
    enclose_point(square[STRIDE],F(0)); checks+=1
    return checks


def build_local_pairs(z,alpha,anchor,guard,deadline,precision=256,
                      rho_exponents=tuple(range(-16,17,4)),progress=None):
    """Export rational upper coefficients and Arb representations thereof.

    The caller supplies Fractions z/alpha and a resource guard. The result
    contains coefficients[(f,s)] and rational_coefficients[(f,s)], each
    a six-tuple in channel order Q,A,E,L,T0,U. These are PAIR channels.
    No marker, regional denominator, source count or source factor appears.
    No numerical job is started merely by importing this module.
    """
    z,alpha=F(z),F(alpha); exponents=tuple(rho_exponents)
    if not 0<z<1 or not 0<alpha<=1: raise ValueError("require0<z<1,0<alpha<=1")
    if precision not in (192,256,384): raise ValueError("precision192/256/384 only")
    if not exponents or len(exponents)>17 or len(set(exponents))!=len(exponents):
        raise ValueError("require1..17 distinct fixed rho exponents")
    if any(type(e) is not int or abs(e)>32 for e in exponents):
        raise ValueError("integer rho exponents in[-32,32]")
    from flint import ctx
    ctx.prec=precision
    emit=(lambda phase,**kw: None) if progress is None else progress
    composition=load("cheap_lane_pair_concavity_triple.py",COMPOSITION_PIN,"directed_pair_composition")
    load("cheap_lane_pair_concavity.py",CORRECTION_PIN,"directed_pair_correction_pin")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"directed_pair_triple")
    caps=composition.build_caps(z,alpha,triple,anchor,guard,deadline)
    v=caps["verifier"]; tests=exact_tests(v,guard)
    shapes=sorted(caps["raw"]); rational_raw=[]; raw=[]; powered=[]
    raw_hash=hashlib.sha256(); power=v.exact(alpha)
    for occ in shapes:
        guard(); frozen=tuple(endpoint(v.positive(x),True) for x in caps["raw"][occ])
        assert all(x>=0 for x in frozen)
        rational_raw.append(frozen)
        raw.append(tuple(v.exact(x) for x in frozen))
        powered.append(tuple(v.arb(0) if x==0 else v.exact(x)**power for x in frozen))
        raw_hash.update((str(occ)+":"+":".join(map(str,frozen))+"\n").encode())
    del caps["raw"],caps["powered"],rational_raw
    geometry,route_report=route_geometry(shapes,guard)
    emit("coherent_raw_table_ready",occupancies=len(shapes),tests=tests,**route_report)
    powered_packet=route_polys(powered,geometry,v,guard)
    old_pair=pair_products(powered_packet,guard); del powered_packet
    raw_packet=route_polys(raw,geometry,v,guard)
    raw_pair=pair_products(raw_packet,guard); del raw_packet,raw
    degrees=[(f,s) for f in range(65) for s in range(65-f)]
    best={degree:[F(0)]*6 for degree in degrees}
    for e in exponents:
        guard(); rho=v.exact(F(2)**e)
        minima=[occupancy_minima(row,rho) for row in powered]
        averaged=route_polys(minima,geometry,v,guard); del minima
        candidate=cross_products(averaged,guard); del averaged
        for f,s in degrees:
            index=f+STRIDE*s
            for j in range(6):
                low=max(F(0),endpoint(candidate[j][index],False))
                if low>best[f,s][j]: best[f,s][j]=low
        del candidate
        emit("rho_correction_complete",rho_exponent=e)
    del powered,geometry
    c_lower=max(F(0),endpoint(2-v.exact(F(2))**power,False))
    assert 0<=c_lower<=1
    rational_coefficients={}; coefficients={}; result_hash=hashlib.sha256()
    corrected_entries=jensen_entries=0
    for f,s in degrees:
        guard(); index=f+STRIDE*s; mass=comb(64,f)*comb(64-f,s)
        jensen_factor=v.exact(F(mass))**v.exact(1-alpha)
        result=[]
        for j in range(6):
            upper=endpoint(v.positive(old_pair[j][index]),True)
            lower=best[f,s][j]
            if 2*lower>upper: raise ArithmeticError("correction exceeds half of coherent old pair")
            corrected=upper-c_lower*lower
            raw_coefficient=v.positive(raw_pair[j][index])
            jensen=endpoint(v.positive(jensen_factor*raw_coefficient**power),True)
            final=min(corrected,jensen)
            assert 0<=final<=upper
            corrected_entries+=corrected<upper; jensen_entries+=jensen<corrected
            result.append(final)
        rational_coefficients[f,s]=tuple(result)
        coefficients[f,s]=tuple(v.exact(x) for x in result)
        result_hash.update((str((f,s))+":"+":".join(map(str,result))+"\n").encode())
    guard()
    report={"scope":"directed local binary-pair coefficient upper bounds only",
        "directed":True,"whole_code_certificate":False,"source_factor_applied":False,
        "support_count_applied":False,"regional_marker_applied":False,
        "pair_channels":list(CHANNELS),"regular_pair":"M*M","terminal_pair":"M*T",
        "pair_slots":64,"packet_slots":32,"packing_stride":STRIDE,
        "precision":precision,"z":str(z),"alpha":str(alpha),"rho_exponents":list(exponents),
        "coherent_raw_table_sha256":raw_hash.hexdigest(),"rational_coefficient_sha256":result_hash.hexdigest(),
        "coefficient_degrees":len(degrees),"coefficients_per_degree":6,
        "old_pair_products":12,"raw_pair_products":12,"correction_products":6*len(exponents),
        "correction_constant_lower":str(c_lower),"corrected_entries":corrected_entries,
        "Jensen_further_improved_entries":jensen_entries,"local_exact_tests":tests,
        "field_tests":caps["tests"],"block_caps":caps["block_report"],
        "route_geometry":route_report,"nonnegative_intersections":v.intersections,
        "composition_pin":COMPOSITION_PIN,"correction_pin":CORRECTION_PIN,
        "triple_pin":composition.TRIPLE_PIN,"generated_bytes":0}
    return {"rational_coefficients":rational_coefficients,"coefficients":coefficients,
            "verifier":v,"core":caps["core"],"report":report}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--precision",type=int,choices=(192,256,384),default=256)
    parser.add_argument("--seconds",type=int,choices=(180,240,300),default=300)
    parser.add_argument("--rho-radius",type=int,choices=(4,8,12,16,20,24,28,32),default=16)
    args=parser.parse_args()
    if not 0<args.z<1 or not 0<args.alpha<=1: parser.error("require0<z<1,0<alpha<=1")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    anchor=load("cheap_lane_cone_cauchy_anchors.py",ANCHOR_PIN,"directed_pair_anchor")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("local pair deadline; no retry")
        if peak()>MEMORY: raise MemoryError("256MiB local pair budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)
    built=build_local_pairs(args.z,args.alpha,anchor,guard,deadline,args.precision,
        tuple(range(-args.rho_radius,args.rho_radius+1,4)),progress)
    guard()
    print(json.dumps({**built["report"],"seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        indent=2),flush=True)


if __name__=="__main__":
    main()
