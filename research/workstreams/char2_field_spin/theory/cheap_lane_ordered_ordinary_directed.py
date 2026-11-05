"""Directed ordinary-Cauchy A/B at actual source weight4096.

Keep the frozen current local ordered-pair bounds, z1/16,alpha119/1280,
and31 rational marker witnesses. Evaluate the positive local H(u,v),J(u,v)
instead of projecting regional coefficients. Ordinary H^64 and H^63 J
bound each target after its OWN exact multinomial/marker denominator.
This retains matching degrees and their routing probabilities; it is not
the discarded route-free maximum-occupancy relaxation.

Only f0 uses the frozen exact top-degree boundary, for a controlled A/B.
All2049 actual marginals and4,198,400 ordered profiles enter the directed
run-mass join. Actual Omega4096 and the output correction appear once.
No local alpha is reapplied. Negative/positive outcomes are both reported;
this source alone is not evidence of a completed numerical certificate.

Parent execution only:180/300s,512MiB,20GiB free,-B,stdout only. No fit,
regional FFT or polynomial products, data files, caches, or encoder edits.
No numerical work runs on import. Other source sizes are not covered.
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
LARGE_PIN="9b2cef00e449dc10c0cd58eca062caf1986ca3f7eed5003087abe32e53e6b6d6"
MEMORY=512*(1<<20)


def load_large():
    path=Path(__file__).with_name("cheap_lane_ordered_large_directed.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=LARGE_PIN:
        raise RuntimeError("frozen large-weight directed driver mismatch")
    spec=importlib.util.spec_from_file_location("ordinary_frozen_large",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def evaluate_pairs(v,coefficients,marker,guard):
    """Positive ordinary H/J evaluation; coefficients already include alpha."""
    u,w=map(F,marker)
    if min(u,w)<=0: raise ValueError("strictly positive rational markers required")
    up=[v.exact(F(1))]; wp=[v.exact(F(1))]
    uu,ww=v.exact(u),v.exact(w)
    for _ in range(64):
        up.append(up[-1]*uu); wp.append(wp[-1]*ww)
    values=[v.exact(F(0)) for _ in range(6)]
    for (f,s),row in coefficients.items():
        guard()
        if min(f,s)<0 or f+s>64 or len(row)!=6:
            raise ValueError("expected six pair channels with f+s<=64")
        scale=up[f]*wp[s]
        for k in range(6): values[k]+=row[k]*scale
    values=tuple(v.positive(x) for x in values)
    return [[values[0],values[1]],[values[2],values[3]]],[[values[4]],[values[5]]]


def exact_tests(guard):
    """Independent rational ordinary/exact-coefficient/terminal controls."""
    def mm(a,b):
        return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                           for j in range(len(b[0]))) for i in range(len(a)))
    def pm(a,b):
        out={}
        for (f,s),x in a.items():
            for (g,t),y in b.items():
                key=f+g,s+t; term=mm(x,y)
                if key not in out: out[key]=tuple(tuple(F(0) for _ in row) for row in term)
                out[key]=tuple(tuple(a+b for a,b in zip(ar,br)) for ar,br in zip(out[key],term))
        return out
    h={(0,0):((F(1,3),F(2,5)),(F(3,7),F(5,11))),
       (1,0):((F(2,3),F(1,7)),(F(5,13),F(3,11))),
       (0,2):((F(7,17),F(2,19)),(F(3,23),F(11,29)))}
    j={(0,0):((F(2,7),),(F(3,11),)),
       (1,0):((F(5,13),),(F(7,17),)),
       (0,2):((F(11,19),),(F(13,23),))}
    rows={key:tuple(x for row in h[key] for x in row)+tuple(x for row in j[key] for x in row)
          for key in h}
    class Exact:
        @staticmethod
        def exact(x): return F(x)
        @staticmethod
        def positive(x):
            assert x>=0
            return x
    def evaluate(poly,marker):
        first=next(iter(poly.values())); u,w=marker
        return tuple(tuple(sum((matrix[i][k]*u**f*w**s for (f,s),matrix in poly.items()),F(0))
                           for k in range(len(first[0]))) for i in range(len(first)))
    checks=0
    for marker in ((F(2,3),F(3,5)),(F(7,4),F(5,2)),(F(1),F(1))):
        guard(); ph,pj=evaluate_pairs(Exact(),rows,marker,guard)
        ph,pj=tuple(map(tuple,ph)),tuple(map(tuple,pj))
        assert ph==evaluate(h,marker) and pj==evaluate(j,marker); checks+=2
        rr,tt=h,j; er,et=ph,pj; top_h,top_j=h[0,2],j[0,2]
        tr,ttop=top_h,top_j
        for step in range(2):
            guard(); tt=pm(rr,tt); rr=pm(rr,rr)
            et=mm(er,et); er=mm(er,er)
            ttop=mm(tr,ttop); tr=mm(tr,tr)
            slots=2*(1<<(step+1))
            for poly,ordinary,top in ((rr,er,tr),(tt,et,ttop)):
                assert evaluate(poly,marker)==ordinary; checks+=1
                assert poly[0,slots]==top; checks+=1
                for (f,s),value in poly.items():
                    assert f+s<=slots
                    mass=comb(slots,f)*comb(slots-f,s)
                    denominator=F(mass)*marker[0]**f*marker[1]**s
                    for i,row in enumerate(value):
                        for k,x in enumerate(row):
                            assert ordinary[i][k]/denominator>=x/mass; checks+=1
    assert mm(h[0,0],h[1,0])!=mm(h[1,0],h[0,0]); checks+=1
    return {"exact_fraction_checks":checks,"ordinary_bounds_every_target_coefficient":True,
        "positive_marker_cancellation":True,"noncommuting_terminal_and_top_degree":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    large=load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"ordinary_anchor")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"ordinary_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"ordinary_normalizers")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"ordinary_join")
    screen=large.load("cheap_lane_ordered_large_screen.py",large.SCREEN_PIN,"ordinary_boundary_tests")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("ordinary directed deadline; no incomplete cover")
        if peak()>MEMORY: raise MemoryError("ordinary directed512MiB budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    tests=exact_tests(guard)
    tests["frozen_boundary_extraction_join"]=large.exact_tests(screen,extraction,join,guard)
    built=local.build_local_pairs(large.Z,large.ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    assert large.WEIGHT==4096 and large.POINTS==8192
    coefficients=built["coefficients"]; degrees=built["rational_coefficients"]
    boundary,boundary_report=large.singleton_boundary(v,degrees,join,guard)
    regions=[None]*2049; reports=[]; boundary_checks=0
    for center,lo,hi,_ in large.schedule():
        guard(); marker=large.WITNESSES[center]
        h,j=evaluate_pairs(v,coefficients,marker,guard)
        regular=v.matpow(h,64); terminal=v.matmul(v.matpow(h,63),j)
        numerator=tuple(tuple(v.positive(x) for row in matrix for x in row) for matrix in (regular,terminal))
        assert tuple(map(len,numerator))==(4,2)
        for f in range(lo,hi+1):
            guard(); denominator=extraction.target_denominator(large.WEIGHT,f,marker)
            dd=v.exact(denominator)
            own=tuple(tuple(v.positive(x/dd) for x in row) for row in numerator)
            assert regions[f] is None
            if f==0:
                for upper_row,exact_row in zip(own,boundary):
                    for upper,exact in zip(upper_row,exact_row):
                        assert upper>=exact or upper.overlaps(exact)
                        boundary_checks+=1
            regions[f]=own
        report={"center":center,"first_full":lo,"last_full":hi,
            "u":str(marker[0]),"v":str(marker[1]),"marginals":hi-lo+1,
            "method":"ordinary positive H64/H63J with own exact target denominator"}
        reports.append(report); progress("ordinary_directed_cell_complete",**report)
        del h,j,regular,terminal,numerator,denominator,dd,own
        gc.collect()
    regions[0]=boundary
    assert boundary_checks==6 and all(x is not None for x in regions)
    tests["ordinary_dominates_exact_boundary_entries"]=boundary_checks
    counts=join.run_counts(points=large.POINTS,weight=large.WEIGHT)
    assert sum(counts)==comb(8192,4096) and counts[4096]==2
    h,c_h,eta=anchor.source_factor(large.WEIGHT)
    assert h==16 and c_h==F(1,17)
    source=(v.exact(eta)*v.exact(core.BETA)**4096*v.exact(large.Z)**(-core.LIMIT))**v.exact(large.ALPHA)
    forward,reverse=join.two_orientations(regions[2048],regions[0])
    selected=[{"full":full,"runs":2048,
        "log2_selected_with_full_run_mass_interval":str((source*counts[2048]*value).log()/logtwo),
        "not_a_run_or_size_cover":True}
        for full,value in (((2048,0),forward),((0,2048),reverse))]
    boundary_report["profile_0_0_diagnostic_not_added_to_join"]=True
    boundary_report["log2_profile_0_0_bound_interval"]=str(
        (source*2*join.two_orientations(boundary,boundary,diagonal=True)[0]).log()/logtwo)
    joined=join.directed_join(regions,guard,progress,weight=4096,points=8192)
    assert joined["report"]["ordered_profiles"]==4198400
    total=v.positive(source*joined["bound"])
    worst=max(range(1,4097),key=lambda c:joined["run_bounds"][c-1].upper())
    guard()
    print(json.dumps({"scope":"complete actual4096-group ordinary-Cauchy directed run envelope",
        "directed_arithmetic":True,"all_profiles_at_weight_covered":True,
        "whole_code_certificate":False,"other_weights_covered":False,
        "parameters":{"weight":4096,"z":str(large.Z),"alpha":str(large.ALPHA),
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,
            "memory_budget_bytes":MEMORY,"fixed_marker_cells":len(reports)},
        "regional_method":"ordinary H(u,v)^64 and H(u,v)^63 J(u,v), own target denominators",
        "matching_degrees_retained":True,"route_free_occupancy_relaxation":False,
        "exact_f0_replacement":True,"residue_projection":False,"ordinary_fallbacks":0,
        "source_once_counts_outside_alpha":True,"new_local_alpha":False,
        "log2_union_bound_interval":str(total.log()/logtwo),
        "actual4096_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
        "actual4096_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
        "worst_run":{"runs":worst,"log2_interval":str((source*joined["run_bounds"][worst-1]).log()/logtwo)},
        "selected_controls":selected,"boundary":boundary_report,"tables":reports,"join":joined["report"],
        "source_log2_interval":str(source.log()/logtwo),"joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "tests":tests,"local_coefficients":built["report"],"generated_bytes":0,
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "nonnegative_intersections":v.intersections,"large_directed_pin":LARGE_PIN,
        "pins":{"local":large.LOCAL_PIN,"normalizers":large.REGION_PIN,"join":large.JOIN_PIN,"anchor":large.ANCHOR_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
