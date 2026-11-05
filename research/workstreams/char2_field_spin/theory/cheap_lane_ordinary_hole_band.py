"""Directed ordinary-Cauchy envelope over an inclusive actual-size band.

Use the frozen z1/64,alpha119/1280 and31 hole-marker witnesses. Reverse
the complete local H/J table once. Each marker gives unnormalized B=H64
and V=H63J. For each unordered marker-cell pair, evaluate both finite
ordered chains and minimize its source-tilted log denominator over exact
rational vertices. Concavity of log multinomial plus affine marker/source
terms makes these vertex lower bounds valid throughout the cell.

The final bound is z^(-alpha L) times the exact indexed support count
weighted by eta_r^alpha, times the maximum source-tilted conditional cell
bound. beta^r appears only through the256 denominator tilts. Counts remain
outside alpha. The source parameter r is the ACTUAL size; k=8192-r is only
the hole coordinate. Odd sizes are included. No endpoint replacement,
profile count, marker union factor, or new local alpha is introduced.

Parent execution only:300s/512MiB/20GiB free,-B,stdout,no generated files.
No numerical work runs on import. Numerical receipts are indexed
separately; this source alone proves no numerical band or whole code.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import comb, factorial
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
HOLE_PIN="a430c6c72bdc4c5b19ccf3d99da49fbc0e1385ee6361bc020469f8bc1806b217"
GEOMETRY_PIN="53685cda5b56f70f70820bb195bbd23d8453635a539360257b6d5f5a51dfa01b"
GEOMETRY_FILE="cheap_lane_marker_band_geometry.py"
POINTS,SLOTS,MEMORY=8192,4096,512*(1<<20)


def load(name,pin,label):
    if not isinstance(pin,str) or len(pin)!=64:
        raise RuntimeError("dependency has not been frozen: "+name)
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def schedule(bank,klo,khi):
    """Retain original integer cells; extend only the last, or clip high f."""
    if not 1<=klo<=khi<=4096: raise ValueError("hole range must be1..4096")
    maximum=khi//2; result=[]
    for index,(center,lo,hi,u,w) in enumerate(bank):
        hi=maximum if index==len(bank)-1 else min(hi,maximum)
        if lo<=hi: result.append((center,lo,hi,u,w))
    covered=[f for _,lo,hi,_,_ in result for f in range(lo,hi+1)]
    if covered!=list(range(maximum+1)):
        raise ArithmeticError("marker cells do not partition every integer marginal")
    return tuple(result)


def exact_tests(bank,geometry,join,guard):
    """Rational scaling, indexed counts, cell coverage and terminal controls."""
    checks=0
    for klo,khi in ((1792,2304),(1,1),(641,642),(4095,4095),(2048,2048)):
        rows=schedule(bank,klo,khi)
        cover={f for _,lo,hi,_,_ in rows for f in range(lo,hi+1)}
        for k in range(klo,khi+1):
            guard(); assert set(range(k//2+1))<=cover; checks+=1
    assert schedule(bank,1792,2304)[-1][1:3]==(1023,1152); checks+=1
    for n in (4,6,8):
        weights=[0]*(n+1)
        for mask in range(1<<n): weights[mask.bit_count()]+=1
        for r in range(1,n):
            assert weights[r]==comb(n,r)==weights[n-r]; checks+=1
        assert sum(weights[1:-1])==sum(comb(n,r) for r in range(1,n)); checks+=1
    scale=F(7,5); d0,d1=F(11,3),F(13,4); chain=F(17,19)
    assert scale**256*chain/(d0*d1)**128==chain/((d0/scale)*(d1/scale))**128
    checks+=1
    for k,a,b in ((5,1,2),(6,0,3),(7,2,2)):
        vertices=geometry.vertices(k,k,a,a,b,b,guard)
        assert vertices==((F(k),F(a),F(b)),); checks+=1
    return {"rational_scaling_count_coverage_checks":checks,
        "source_tilt_has_exactly256_factors":True,"odd_sizes_tested":True,
        "singleton_cell_vertices":True,"geometry":geometry.exact_tests(guard),
        "finite_noncommuting_terminal":join.exact_tests(guard)}


class LogDenominators:
    """Cached directed logGamma values; all coordinates remain exact rationals."""
    def __init__(self,v,beta,alpha,guard):
        self.v,self.guard=v,guard
        self.alpha=v.exact(alpha); self.logbeta=v.exact(beta).log()
        self.gamma={}; self.logfactorial=self.lgamma(F(4097)); self.calls=0

    def lgamma(self,x):
        x=F(x)
        if x<=0: raise ValueError("logGamma domain must be positive")
        if x not in self.gamma:
            self.guard(); value=self.v.exact(x).lgamma()
            if not value.is_finite(): raise ArithmeticError("nonfinite directed logGamma")
            self.gamma[x]=value
        return self.gamma[x]

    def value(self,k,a,logu,logw):
        k,a=F(k),F(a); single=k-2*a; empty=4096-k+a
        if min(a,single,empty)<0 or not 1<=k<=4096:
            raise ValueError("invalid exact regional counts")
        self.calls+=1
        return (self.logfactorial-self.lgamma(a+1)-self.lgamma(single+1)
                -self.lgamma(empty+1)+self.v.exact(a)*logu+self.v.exact(single)*logw
                -self.alpha*self.v.exact(F(8192-k,256))*self.logbeta)


def lower_fraction(value):
    if not value.is_finite(): raise ArithmeticError("nonfinite log denominator")
    endpoint=value.lower().fmpq()
    return F(int(endpoint.numerator),int(endpoint.denominator))


def evaluate_cell(v,den,vertices,left_logs,right_logs,left_region,right_region,join,diagonal=False):
    """Return one symmetric denominator minimum and both ordered log bounds."""
    if not vertices: raise ValueError("cannot evaluate an empty cell")
    low=None; low_vertex=None
    for k,a,b in vertices:
        den.guard()
        candidate=lower_fraction(den.value(k,a,*left_logs)+den.value(k,b,*right_logs))
        if low is None or candidate<low: low,low_vertex=candidate,(k,a,b)
    forward,reverse=join.two_orientations(left_region,right_region,diagonal=diagonal)
    bounds=[]
    for direction,numerator in (("forward",forward),("reverse",reverse)):
        if direction=="reverse" and diagonal: continue
        if not numerator.is_finite() or not numerator>0:
            raise ArithmeticError("ordered ordinary chain must be finite and strictly positive")
        bounds.append((direction,numerator.log()-128*v.exact(low)))
    return low,low_vertex,tuple(bounds)


def directed_tests(v,den,extraction,bank,alpha,beta,join,guard):
    checks=0
    for n in (0,1,2,3,7,16,32,4096):
        guard()
        assert den.lgamma(F(n+1)).overlaps(v.exact(F(factorial(n))).log()); checks+=1
    marker=(bank[15][3],bank[15][4]); logu,logw=(v.exact(x).log() for x in marker)
    for k,a in ((5,1),(6,3),(2047,511),(2048,512),(2304,1152)):
        exact=extraction.target_denominator(k,a,marker)
        direct=v.exact(exact).log()-v.exact(alpha)*v.exact(F(8192-k,256))*v.exact(beta).log()
        evaluated=den.value(k,a,logu,logw)
        assert evaluated.overlaps(direct); checks+=1
        low=v.exact(lower_fraction(evaluated))
        assert low<=direct or low.overlaps(direct); checks+=1
    # A singleton geometry cell must reproduce the ordinary normalized
    # finite chain, apart from the deliberate downward endpoint rounding.
    markers=(marker,(bank[16][3],bank[16][4]))
    regions=tuple((tuple(v.exact(F(n,17)) for n in (2+t,3,5,7+t)),
                   tuple(v.exact(F(n,19)) for n in (11+t,13))) for t in (0,1))
    logs=tuple(tuple(v.exact(x).log() for x in m) for m in markers)
    point=(F(5),F(1),F(2))
    low,vertex,bounds=evaluate_cell(v,den,(point,),logs[0],logs[1],regions[0],regions[1],join)
    assert vertex==point; checks+=1
    normalized=[]
    for region,a,m in zip(regions,(1,2),markers):
        divisor=v.exact(extraction.target_denominator(5,a,m))
        normalized.append(tuple(tuple(x/divisor for x in row) for row in region))
    direct_chains=join.two_orientations(*normalized)
    exact_log_den=den.value(5,1,*logs[0])+den.value(5,2,*logs[1])
    for (_,bound),chain in zip(bounds,direct_chains):
        direct=chain.log()+v.exact(alpha)*8187*v.exact(beta).log()
        assert bound>=direct or bound.overlaps(direct); checks+=1
        # Removing only endpoint slack must enclose the identical formula.
        adjusted=bound-128*(exact_log_den-v.exact(low))
        assert adjusted.overlaps(direct); checks+=1
    return {"loggamma_integer_and_singleton_cell_checks":checks,
        "singleton_cell_finite_chain_and_both_terminal_orders":True,
        "fractional_vertices_not_rounded":True,"lower_endpoint_is_exact_rational":True}


def outer_source_sum(v,anchor,alpha,rlo,rhi,guard):
    """Counts are exact and outside alpha; beta is intentionally absent."""
    value=comb(8192,rlo); exact_sum=0; independent_sum=0
    weighted=v.exact(F(0)); count_hash=hashlib.sha256(); factors={}
    min_c=None; min_c_r=None
    for r in range(rlo,rhi+1):
        guard(); direct=comb(8192,r)
        assert value==direct
        exact_sum+=value; independent_sum+=direct
        count_hash.update((str(r)+":"+str(value)+"\n").encode("ascii"))
        h,c_h,eta=anchor.source_factor(r)
        if c_h<=0 or eta<=0: raise ArithmeticError("invalid actual-size source factor")
        weighted+=value*v.exact(eta)**v.exact(alpha)
        if min_c is None or c_h<min_c: min_c,min_c_r=c_h,r
        entry=factors.setdefault(h,{"first_actual":r,"last_actual":r,"sizes":0})
        entry["last_actual"]=r; entry["sizes"]+=1
        if r<rhi:
            value,remainder=divmod(value*(8192-r),r+1)
            assert remainder==0
    assert exact_sum==independent_sum and sum(x["sizes"] for x in factors.values())==rhi-rlo+1
    return weighted,{"actual_first":rlo,"actual_last":rhi,"sizes":rhi-rlo+1,
        "odd_sizes":sum(r%2 for r in range(rlo,rhi+1)),"indexed_support_count":str(exact_sum),
        "indexed_count_sha256":count_hash.hexdigest(),"independent_binomial_checksum":True,
        "log2_indexed_support_count_interval":str(v.exact(F(exact_sum)).log()/v.arb(2).log()),
        "source_h_ranges":factors,"minimum_c_h":str(min_c),"minimum_c_h_actual":min_c_r,
        "eta_recomputed_for_each_actual_size":True,"counts_outside_alpha":True,
        "beta_factor_applied_here":False}


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
    hole=load("cheap_lane_ordinary_hole_directed.py",HOLE_PIN,"ordinary_band_frozen_hole")
    geometry=load(GEOMETRY_FILE,GEOMETRY_PIN,"ordinary_band_vertices")
    ordinary=hole.load("cheap_lane_ordered_ordinary_directed.py",hole.ORDINARY_PIN,"ordinary_band_evaluation")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"ordinary_band_anchor")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"ordinary_band_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"ordinary_band_normalizer")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"ordinary_band_join")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("ordinary band deadline; no incomplete band claim")
        if peak()>MEMORY: raise MemoryError("ordinary band512MiB budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    klo,khi=8192-args.last,8192-args.first
    cells=schedule(hole.BANK,klo,khi)
    tests=exact_tests(hole.BANK,geometry,join,guard)
    tests["ordinary"]=ordinary.exact_tests(guard)
    built=local.build_local_pairs(hole.Z,hole.ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    coefficients={(64-f-s,s):row for (f,s),row in built["coefficients"].items()}
    assert len(coefficients)==len(built["coefficients"])==2145
    assert all(row is built["coefficients"][64-f-s,s] for (f,s),row in coefficients.items())
    den=LogDenominators(v,core.BETA,hole.ALPHA,guard)
    tests["directed"]=directed_tests(v,den,extraction,hole.BANK,hole.ALPHA,core.BETA,join,guard)
    numerators=[]; marker_logs=[]
    for center,lo,hi,u,w in cells:
        guard(); h,j=ordinary.evaluate_pairs(v,coefficients,(u,w),guard)
        rr=v.matpow(h,64); tt=v.matmul(v.matpow(h,63),j)
        numerators.append((tuple(v.positive(x) for row in rr for x in row),
                           tuple(v.positive(x) for row in tt for x in row)))
        marker_logs.append((v.exact(u).log(),v.exact(w).log()))
        progress("ordinary_band_marker",center=center,first_full=lo,last_full=hi)
    maximum=None; worst=None; worst_upper=None; global_min=None; global_min_at=None
    vertex_hash=hashlib.sha256(); vertex_count=nonempty=empty=ordered=0
    cell_reports=[]
    for i,left in enumerate(cells):
        for j in range(i,len(cells)):
            guard(); right=cells[j]
            vertices=geometry.vertices(klo,khi,left[1],left[2],right[1],right[2],guard)
            vertex_hash.update((str(i)+","+str(j)+":"+repr(vertices)+"\n").encode("ascii"))
            if not vertices:
                empty+=1
                continue
            nonempty+=1; vertex_count+=len(vertices); ordered+=1 if i==j else 2
            low,low_vertex,bounds=evaluate_cell(v,den,vertices,marker_logs[i],marker_logs[j],
                numerators[i],numerators[j],join,diagonal=i==j)
            if global_min is None or low<global_min:
                global_min,global_min_at=low,(left[0],right[0],low_vertex)
            pair_bound=None; winner=None
            for direction,logbound in bounds:
                maximum=logbound if maximum is None else maximum.max(logbound)
                if pair_bound is None or logbound.upper()>pair_bound.upper():
                    pair_bound,winner=logbound,direction
                if worst_upper is None or logbound.upper()>worst_upper:
                    worst_upper=logbound.upper()
                    worst={"centers":[left[0],right[0]],"order":direction,
                        "first_interval":[left[1],left[2]],"second_interval":[right[1],right[2]],
                        "minimum_vertex_hole_full_full":list(map(str,low_vertex)),
                        "actual_size_at_minimum_vertex":str(8192-low_vertex[0]),
                        "lower_log_denominator_sum":str(v.exact(low)),
                        "log2_conditional_cell_upper_interval":str(logbound/logtwo)}
            cell_reports.append((pair_bound.upper(),{"centers":[left[0],right[0]],"vertices":len(vertices),
                "larger_order":winner,"log2_cell_upper_interval":str(pair_bound/logtwo)}))
        progress("ordinary_band_geometry_row",row=i,center=left[0],nonempty_cell_pairs=nonempty,
                 empty_cell_pairs=empty,vertices=vertex_count)
    if maximum is None: raise ArithmeticError("no nonempty source cells")
    weighted,source_report=outer_source_sum(v,anchor,hole.ALPHA,args.first,args.last,guard)
    if not weighted>0: raise ArithmeticError("nonpositive weighted source sum")
    log_count=weighted.log()
    log_output=-v.exact(hole.ALPHA)*core.LIMIT*v.exact(hole.Z).log()
    log_outer=log_count+log_output
    log_total=log_outer+maximum
    guard()
    print(json.dumps({"scope":"all indexed supports in the inclusive actual-size band and all histories at35%",
        "directed_arithmetic":True,"whole_code_certificate":False,"all_integer_sizes_in_band":True,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"hole_first":klo,"hole_last":khi,
            "z":str(hole.Z),"alpha":str(hole.ALPHA),"working_precision":built["report"]["working_precision"],
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "log2_union_bound_interval":str(log_total/logtwo),
        "band_below_2_pow_minus_40":bool(log_total<-40*logtwo),
        "band_below_2_pow_minus_60":bool(log_total<-60*logtwo),
        "worst_cell":worst,"log2_outer_factor_interval":str(log_outer/logtwo),
        "log2_eta_weighted_support_count_interval":str(log_count/logtwo),
        "log2_output_correction_interval":str(log_output/logtwo),
        "geometry":{"nonempty_unordered_cells":nonempty,"empty_unordered_cells":empty,
            "ordered_cells":ordered,"rational_vertices":vertex_count,"vertices_sha256":vertex_hash.hexdigest(),
            "global_minimum_log_denominator_sum":str(v.exact(global_min)),
            "global_minimum_at":str(global_min_at),"loggamma_cache_entries":len(den.gamma),
            "log_denominator_evaluations":den.calls,"fractional_vertices_preserved":True},
        "source_counts":source_report,"source_beta_only_in_denominator_tilt":True,
        "output_factor_once":True,"degree_reversals":1,"new_local_alpha":False,
        "boundary_replacement":False,"extra_profile_or_orientation_factor":False,
        "integer_marginal_coverage":khi//2+1,
        "bank":[{"center":c,"lo":lo,"hi":hi,"u":str(u),"v":str(w)} for c,lo,hi,u,w in cells],
        "largest_cell_controls":[entry for _,entry in sorted(cell_reports,key=lambda x:x[0],reverse=True)[:8]],
        "tests":tests,"local_coefficients":built["report"],
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "nonnegative_intersections":v.intersections,"pins":{"hole":HOLE_PIN,"geometry":GEOMETRY_PIN,
            "ordinary":hole.ORDINARY_PIN,"local":large.LOCAL_PIN,"anchor":large.ANCHOR_PIN,
            "join":large.JOIN_PIN,"normalizer":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
