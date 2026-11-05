"""Directed hole-band envelope with normalized marker cells.

Keep the frozen current local caps, z1/64,alpha119/1280 and31 rational
markers. Assign a matching target a at hole size k=8192-r by x=2048*a/k.
The continuous cuts are0, each original cell's upper integer endpoint
plus1/2, and1024. Ties go to the lower cell for the integer coverage check;
the proof uses closed cells, whose shared boundaries do not add counts.

Each normalized cell is a polytope in(k,a,b), since its cuts are fixed.
The unchanged directed log-denominator and vertex-minimum calculation
therefore applies. The wrapper reverses the complete local table once,
retains both ordered terminal numerators, and uses the actual-r source
factor for every integer size, including odd sizes. Counts stay outside
alpha. No endpoint replacement, extra profile factor, new markers, or
new local power is introduced. Numerical receipts are indexed separately.

Parent execution only:300s/512MiB/20GiB free,-B,stdout,no generated files.
No computation runs on import. This is one band, not whole-code closure.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
BAND_PIN="ca4147fc8eb19aed544aca76eabc1ace357974b2bcd026474cac1465afb753f4"
GEOMETRY_PIN="00a4fde7841efed7c9c58b199e958c9a799bb190e4724ee9332d197e8ad9dee4"
GEOMETRY_FILE="cheap_lane_normalized_marker_geometry.py"
MEMORY=512*(1<<20)


def load_band():
    path=Path(__file__).with_name("cheap_lane_ordinary_hole_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=BAND_PIN:
        raise RuntimeError("frozen ordinary hole-band source mismatch")
    spec=importlib.util.spec_from_file_location("normalized_frozen_band",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def normalized_cells(bank):
    """All fixed witnesses, continuous closed cuts in x=2048*a/k."""
    result=[]; lower=F(0)
    for index,(center,lo,hi,u,w) in enumerate(bank):
        upper=F(1024) if index==len(bank)-1 else F(2*hi+1,2)
        if not 0<=lower<upper<=1024 or min(u,w)<=0:
            raise ValueError("invalid fixed normalized marker bank")
        result.append((center,lower,upper,u,w))
        lower=upper
    if len(result)!=31 or result[0][1]!=0 or result[-1][2]!=1024:
        raise ArithmeticError("normalized cells must cover0..1024 with all31 markers")
    return tuple(result)


def coverage(cells,klo,khi,guard):
    """Check every integer k without enumerating all matching targets."""
    if not 1<=klo<=khi<=4096: raise ValueError("hole range must be1..4096")
    total=0; active=[0]*len(cells); checksum=hashlib.sha256()
    for k in range(klo,khi+1):
        guard(); next_full=0; count=0
        for index,(center,lower,upper,_,_) in enumerate(cells):
            # Lower boundary is strict except in the first cell. At a cut,
            # the integer target belongs to the lower-indexed witness.
            lo=0 if index==0 else (k*lower/2048).__floor__()+1
            hi=min(k//2,(k*upper/2048).__floor__())
            if lo<=hi:
                assert lo==next_full
                assert lower<=F(2048*lo,k)<=F(2048*hi,k)<=upper
                if index: assert lower<F(2048*lo,k)
                size=hi-lo+1; count+=size; active[index]+=size; next_full=hi+1
        assert count==k//2+1 and next_full==k//2+1
        total+=count; checksum.update((str(k)+":"+str(count)+"\n").encode("ascii"))
    return {"hole_sizes_checked":khi-klo+1,"integer_marginals_checked":total,
        "coverage_sha256":checksum.hexdigest(),"tie_policy":"lower cell",
        "continuous_domain":["0","1024"],
        "marginals_by_center":{str(row[0]):n for row,n in zip(cells,active)}}


def normalized_tests(cells,bank,geometry,guard):
    checks=0
    # At reference k2048, the continuous cuts reproduce the original
    # integer marker assignment exactly, not a refitted witness bank.
    for index,(cell,old) in enumerate(zip(cells,bank)):
        center,lower,upper,u,w=cell
        lo=0 if index==0 else lower.__floor__()+1
        hi=upper.__floor__()
        assert (center,lo,hi,u,w)==old; checks+=1
    reports=[]
    for klo,khi in ((1,9),(1792,1793),(2047,2049),(2303,2304),(4095,4096)):
        reports.append(coverage(cells,klo,khi,guard)); checks+=1
    for k,a,b in ((5,1,2),(6,0,3),(7,2,2)):
        x,y=F(2048*a,k),F(2048*b,k)
        vertices=geometry.vertices(k,k,x,x,y,y,guard)
        assert vertices==((F(k),F(a),F(b)),); checks+=1
    return {"normalized_assignment_checks":checks,"reference2048_assignment_unchanged":True,
        "coverage_controls":reports,"normalized_geometry":geometry.exact_tests(guard)}


def numerator_tables(v,coefficients,cells,ordinary,guard,progress):
    """No target/source scaling: ordinary B=H64 and V=H63J only."""
    numerators=[]; logs=[]
    for center,lower,upper,u,w in cells:
        guard(); h,j=ordinary.evaluate_pairs(v,coefficients,(u,w),guard)
        rr=v.matpow(h,64); tt=v.matmul(v.matpow(h,63),j)
        numerators.append((tuple(v.positive(x) for row in rr for x in row),
                           tuple(v.positive(x) for row in tt for x in row)))
        logs.append((v.exact(u).log(),v.exact(w).log()))
        if progress: progress("normalized_band_marker",center=center,normalized_lower=str(lower),
                              normalized_upper=str(upper))
    return tuple(numerators),tuple(logs)


def cell_envelope(v,band,geometry,den,numerators,marker_logs,cells,klo,khi,join,guard,progress):
    """Closed normalized geometry once per unordered pair, both orders."""
    maximum=None; worst=None; worst_upper=None; global_min=None; global_min_at=None
    vertex_hash=hashlib.sha256(); vertex_count=nonempty=empty=ordered=0
    controls=[]; logtwo=v.arb(2).log()
    for i,left in enumerate(cells):
        for j in range(i,len(cells)):
            guard(); right=cells[j]
            vertices=geometry.vertices(klo,khi,left[1],left[2],right[1],right[2],guard)
            vertex_hash.update((str(i)+","+str(j)+":"+repr(vertices)+"\n").encode("ascii"))
            if not vertices:
                empty+=1
                continue
            nonempty+=1; vertex_count+=len(vertices); ordered+=1 if i==j else 2
            low,point,bounds=band.evaluate_cell(v,den,vertices,marker_logs[i],marker_logs[j],
                numerators[i],numerators[j],join,diagonal=i==j)
            if global_min is None or low<global_min:
                global_min,global_min_at=low,(left[0],right[0],point)
            pair_upper=None; winner=None
            for direction,logbound in bounds:
                maximum=logbound if maximum is None else maximum.max(logbound)
                if pair_upper is None or logbound.upper()>pair_upper.upper():
                    pair_upper,winner=logbound,direction
                if worst_upper is None or logbound.upper()>worst_upper:
                    worst_upper=logbound.upper()
                    k,a,b=point
                    worst={"centers":[left[0],right[0]],"order":direction,
                        "normalized_first_cut":list(map(str,left[1:3])),
                        "normalized_second_cut":list(map(str,right[1:3])),
                        "minimum_vertex_hole_full_full":list(map(str,point)),
                        "normalized_vertex":[str(F(2048*a,k)),str(F(2048*b,k))],
                        "actual_size_at_minimum_vertex":str(8192-k),
                        "lower_log_denominator_sum":str(v.exact(low)),
                        "log2_conditional_cell_upper_interval":str(logbound/logtwo)}
            controls.append((pair_upper.upper(),{"centers":[left[0],right[0]],"vertices":len(vertices),
                "larger_order":winner,"log2_cell_upper_interval":str(pair_upper/logtwo)}))
        if progress: progress("normalized_band_geometry_row",row=i,center=left[0],
            nonempty_cell_pairs=nonempty,empty_cell_pairs=empty,vertices=vertex_count)
    if maximum is None: raise ArithmeticError("no nonempty normalized cells")
    report={"nonempty_unordered_cells":nonempty,"empty_unordered_cells":empty,
        "ordered_cells":ordered,"rational_vertices":vertex_count,"vertices_sha256":vertex_hash.hexdigest(),
        "global_minimum_log_denominator_sum":str(v.exact(global_min)),
        "global_minimum_at":str(global_min_at),"loggamma_cache_entries":len(den.gamma),
        "log_denominator_evaluations":den.calls,"fractional_vertices_preserved":True,
        "largest_cell_controls":[row for _,row in sorted(controls,key=lambda x:x[0],reverse=True)[:8]]}
    return maximum,worst,report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first",type=int,default=5888)
    parser.add_argument("--last",type=int,default=6400)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    args=parser.parse_args()
    if not 4097<=args.first<=args.last<=7551:
        raise ValueError("require inclusive actual band4097..7551")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    band=load_band()
    geometry=band.load(GEOMETRY_FILE,GEOMETRY_PIN,"normalized_band_geometry")
    old_geometry=band.load(band.GEOMETRY_FILE,band.GEOMETRY_PIN,"normalized_frozen_rectangle_geometry")
    hole=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"normalized_frozen_hole")
    ordinary=hole.load("cheap_lane_ordered_ordinary_directed.py",hole.ORDINARY_PIN,"normalized_ordinary")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"normalized_anchor")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"normalized_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"normalized_denominator")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"normalized_join")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("normalized band deadline; no incomplete cover")
        if peak()>MEMORY: raise MemoryError("normalized band512MiB budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    klo,khi=8192-args.last,8192-args.first
    cells=normalized_cells(hole.BANK)
    cover=coverage(cells,klo,khi,guard)
    tests={"frozen_band":band.exact_tests(hole.BANK,old_geometry,join,guard),
        "ordinary":ordinary.exact_tests(guard),"normalized":normalized_tests(cells,hole.BANK,geometry,guard)}
    built=local.build_local_pairs(hole.Z,hole.ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    coefficients={(64-f-s,s):row for (f,s),row in built["coefficients"].items()}
    assert len(coefficients)==len(built["coefficients"])==2145
    assert all(row is built["coefficients"][64-f-s,s] for (f,s),row in coefficients.items())
    den=band.LogDenominators(v,core.BETA,hole.ALPHA,guard)
    tests["directed"]=band.directed_tests(v,den,extraction,hole.BANK,hole.ALPHA,core.BETA,join,guard)
    numerators,marker_logs=numerator_tables(v,coefficients,cells,ordinary,guard,progress)
    maximum,worst,geometry_report=cell_envelope(v,band,geometry,den,numerators,marker_logs,
        cells,klo,khi,join,guard,progress)
    weighted,source_report=band.outer_source_sum(v,anchor,hole.ALPHA,args.first,args.last,guard)
    if not weighted>0: raise ArithmeticError("nonpositive weighted source sum")
    log_count=weighted.log(); log_output=-v.exact(hole.ALPHA)*core.LIMIT*v.exact(hole.Z).log()
    log_total=log_count+log_output+maximum
    guard()
    print(json.dumps({"scope":"all indexed supports in the inclusive actual-size band and all histories at35%",
        "directed_arithmetic":True,"whole_code_certificate":False,"all_integer_sizes_in_band":True,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"hole_first":klo,"hole_last":khi,
            "z":str(hole.Z),"alpha":str(hole.ALPHA),"working_precision":built["report"]["working_precision"],
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "cell_coordinate":"x=2048*hole_full/hole_weight","continuous_cut_rule":"old hi+1/2; endpoints0,1024",
        "log2_union_bound_interval":str(log_total/logtwo),
        "band_below_2_pow_minus_40":bool(log_total<-40*logtwo),
        "band_below_2_pow_minus_60":bool(log_total<-60*logtwo),
        "log2_eta_weighted_support_count_interval":str(log_count/logtwo),
        "log2_output_correction_interval":str(log_output/logtwo),
        "worst_cell":worst,"geometry":geometry_report,"coverage":cover,"source_counts":source_report,
        "source_beta_only_in_denominator_tilt":True,"output_factor_once":True,"degree_reversals":1,
        "new_local_alpha":False,"boundary_replacement":False,"extra_profile_or_orientation_factor":False,
        "bank":[{"center":c,"lower":str(lo),"upper":str(hi),"u":str(u),"v":str(w)} for c,lo,hi,u,w in cells],
        "tests":tests,"local_coefficients":built["report"],"elapsed_seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"generated_bytes":0,"nonnegative_intersections":v.intersections,
        "pins":{"band":BAND_PIN,"normalized_geometry":GEOMETRY_PIN,"old_geometry":band.GEOMETRY_PIN,
            "hole":band.HOLE_PIN,"ordinary":hole.ORDINARY_PIN,"local":large.LOCAL_PIN,
            "anchor":large.ANCHOR_PIN,"join":large.JOIN_PIN,"normalizer":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
