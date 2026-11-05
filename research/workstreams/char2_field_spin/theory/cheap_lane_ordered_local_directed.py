r"""Directed local ordered-ratio pair coefficients, with residual-rank returns.

Freeze every powered local cap UPWARD to an exact multiple of2^-256.
The old pair and every ordered correction use this same rational table.
For its two path products P,Q, G(P,Q)=max(P,Q)+(2^alpha-1)min(P,Q)
is coordinatewise increasing and bounds the original pooled target.
Consequently P+Q-kappa_lower*L remains an upper bound whenever
0<=L<=min(P,Q) and kappa_lower<=2-2^alpha. This does not assume that
the bucket formula itself is monotone under independent atom rounding.
An independent degree-Jensen bound from upward-frozen RAW caps may be
intersected afterward: both bounds concern the same actual pair target.

All ratio comparisons use integer powered atoms. Shared nearest-log
cells certify exact off-diagonal minima; same-cell terms use rho=m/n.
Write the diagonal factors as min(m*a,n*c),min(n*b,m*d), then divide
their product by m*n. One exact LCM clears all grid denominators.
Exact dyadic route weights and carry-safe stride65 fmpz_poly products
give the integrated correction with no interval or FLOAT subtraction.
Only kappa_lower and the independent Jensen upper use Arb endpoints.

The ordered correction dominates every witness at a grid point. Within
one cell centered rho, its normalized witness is min(r,s,rho,r*s/rho).
The cell boundaries imply r*s>=rho_prev*rho and r*s<=rho*rho_next;
therefore no smaller or larger grid point improves this witness. Across
different cells the minimum is exact. The129 grid contains the old
dyadic33 and9 grids, so their additional products need not be rebuilt.

Export API matches cheap_lane_pair_local_directed.build_local_pairs:
six rational upper coefficients at each(f,s),f+s<=64, for H=M*M and
J=M*T. No marker, regional normalization, source count/factor, or extra
alpha is applied. This is only a local bound, not a whole-code result.
No numerical work on import. Parent runs all jobs. Standalone budgets:
180/300s,512MiB,20GiB free, stdout only, zero generated files.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
import json
from math import comb, gcd
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
LOCAL_PIN="7284d1cd125d2cf1d1a3c1a38aa3dfd8b2ef7ceae4386ea2975406ea3e64f81f"
RANK_PIN="3efa28f205a01cafe47a8058ca4d3a69e065d5e22e75ea00a8fe79ed4a5dcf7c"
ORDERED_PIN="f816addd50969a83ad951567944d905938ece966c9260175e3b4a08b1be53255"
ANCHOR_PIN="319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
STRIDE,PACKET_LENGTH,POWER_BITS=65,2081,256
MEMORY=512*(1<<20)
PATHS=((0,0,1,2),(0,1,1,3),(2,0,3,2),(2,1,3,3),(0,4,1,5),(2,4,3,5))
CROSS=((0,0),(0,1),(1,0),(1,1),(0,2),(1,2))


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def ceil_dyadic(value,bits=POWER_BITS):
    """Integer numerator of an upward multiple of2^-bits; preserve zero."""
    if value<0: raise ValueError("nonnegative rational required")
    numerator=value.numerator<<bits
    q,r=divmod(numerator,value.denominator)
    return q+bool(r)


def denominator_bits(value):
    d=value.denominator
    if d&(d-1): raise ArithmeticError("Arb endpoint was not an exact dyadic")
    return d.bit_length()-1


def nearest_integer(a,c,boundaries):
    """Classify c/a by exact squared comparisons; ties go lower."""
    if a<0 or c<0: raise ValueError("nonnegative integer atoms required")
    if not a or not c: return None
    aa,cc=a*a,c*c; lo,hi=0,len(boundaries)
    while lo<hi:
        middle=(lo+hi)//2; numerator,denominator=boundaries[middle]
        if cc*denominator<=aa*numerator: hi=middle
        else: lo=middle+1
    return lo


def grid_data(grid):
    assert grid and tuple(sorted(set(grid)))==tuple(grid) and grid[0]>0
    boundaries=tuple(((a*b).numerator,(a*b).denominator) for a,b in zip(grid,grid[1:]))
    denominator=1
    for rho in grid:
        x=rho.numerator*rho.denominator
        denominator=denominator//gcd(denominator,x)*x
    return boundaries,denominator


def aggregate_rows(rows,geometry,poly,guard):
    """Exact route numerators; common denominator2^32 stays outside."""
    coeff=[[0]*PACKET_LENGTH for _ in range(6)]
    for row,records in zip(rows,geometry):
        guard()
        nonzero=tuple((j,value) for j,value in enumerate(row) if value)
        for index,weight in records:
            for j,value in nonzero: coeff[j][index]+=weight*value
    return tuple(poly(x) for x in coeff)


def pair_products(packet,poly,guard,counters):
    out=[]
    for a,b,c,d in PATHS:
        guard(); out.append(packet[a]*packet[b]+packet[c]*packet[d])
        counters["pair_polynomial_products"]+=2
    return tuple(out)


def bucket_products(rows,geometry,grid,poly,guard,emit):
    """Exact positive ordered correction numerator over the returned LCM."""
    boundaries,common=grid_data(grid); count=len(grid)
    left=[[[] for _ in grid] for _ in range(2)]
    right=[[[] for _ in grid] for _ in range(3)]
    zeros=[0]*5; clipped=[0]*5; classified=[0]*5
    for shape,row in enumerate(rows):
        guard()
        for channel in range(2):
            a,c=row[2*channel:2*channel+2]
            index=nearest_integer(a,c,boundaries)
            if index is None: zeros[channel]+=1; continue
            rho=grid[index]; m,n=rho.numerator,rho.denominator
            left[channel][index].append((shape,(a,c,min(m*a,n*c))))
            classified[channel]+=1
            clipped[channel]+=c*grid[0].denominator<a*grid[0].numerator or c*grid[-1].denominator>a*grid[-1].numerator
        for channel,(bi,di) in enumerate(((0,2),(1,3),(4,5))):
            b,d=row[bi],row[di]; index=nearest_integer(d,b,boundaries)
            if index is None: zeros[2+channel]+=1; continue
            rho=grid[index]; m,n=rho.numerator,rho.denominator
            right[channel][index].append((shape,(b,d,min(n*b,m*d))))
            classified[2+channel]+=1
            clipped[2+channel]+=b*grid[0].denominator<d*grid[0].numerator or b*grid[-1].denominator>d*grid[-1].numerator
    assert all(classified[i]+zeros[i]==len(rows) for i in range(5))
    zero=poly([])
    def route(entries):
        if not entries: return (zero,zero,zero)
        coefficients=[[0]*PACKET_LENGTH for _ in range(3)]
        for shape,values in entries:
            guard()
            for index,weight in geometry[shape]:
                coefficients[0][index]+=weight*values[0]
                coefficients[1][index]+=weight*values[1]
                coefficients[2][index]+=weight*values[2]
        return tuple(poly(x) for x in coefficients)
    lp=[]; rp=[]
    for i in range(count):
        guard(); l=[route(left[j][i]) for j in range(2)]
        r=[route(right[j][i]) for j in range(3)]
        lp.append(tuple(l[j][k] for k in range(3) for j in range(2)))
        rp.append(tuple(r[j][k] for k in range(3) for j in range(3)))
    del left,right
    suffix=[None]*(count+1); suffix[count]=(zero,zero,zero)
    for i in range(count-1,-1,-1):
        guard(); suffix[i]=tuple(rp[i][3+j]+suffix[i+1][j] for j in range(3))
    prefix=[zero,zero,zero]; correction=[zero]*6; products=occupied=0
    def multiply(a,b):
        nonlocal products
        guard()
        if not len(a) or not len(b): return zero
        products+=1; return a*b
    for i,rho in enumerate(grid):
        guard(); l,r=lp[i],rp[i]; m,n=rho.numerator,rho.denominator
        diagonal_scale=common//(m*n); assert diagonal_scale*(m*n)==common
        if any(len(x) for x in l[:2]):
            occupied+=1
            for channel,(a,b) in enumerate(CROSS):
                off=multiply(l[a],prefix[b])+multiply(l[2+a],suffix[i+1][b])
                diagonal=multiply(l[4+a],r[6+b])
                correction[channel]=correction[channel]+common*off+diagonal_scale*diagonal
        prefix=[prefix[j]+r[j] for j in range(3)]
        if i%16==0 or i+1==count:
            emit("directed_ordered_bucket_complete",through_bucket=i+1,
                 grid_size=count,integer_products=products)
    return tuple(correction),common,{"grid_size":count,"occupied_left_buckets":occupied,
        "integer_correction_products":products,"diagonal_denominator_lcm":common,
        "classified_by_row_column":classified,"zeros_by_row_column":zeros,
        "clipped_by_row_column":clipped,"exact_integer_ordering":True}


def exact_tests(base,v,poly,guard):
    """Independent rational route reference and powered-freeze controls."""
    checks=0
    for x in (F(0),F(1,3),F(1,1<<300),F(2),F(7,5)):
        for bits in (8,256):
            n=ceil_dyadic(x,bits); assert F(n,1<<bits)>=x
            assert (n==0)==(x==0)
            if n: assert F(n-1,1<<bits)<x
            checks+=3
    grid=(F(1),F(4),F(16)); boundaries,common=grid_data(grid)
    for ratio,expected in ((F(1,8),0),(F(2),0),(F(8),1),(F(128),2)):
        assert nearest_integer(ratio.denominator,ratio.numerator,boundaries)==expected; checks+=1
    atoms=(((0,0),2,(1,2,3,6,4,8)),((1,0),3,(4,1,0,5,2,0)),
           ((0,1),1,(2,10,6,2,0,7)),((1,1),5,(1,128,128,1,3,7)))
    rows=[x[2] for x in atoms]
    geometry=[[(degree[0]+STRIDE*degree[1],weight)] for degree,weight,row in atoms]
    coefficients=aggregate_rows(rows,geometry,poly,guard)
    counters={"pair_polynomial_products":0}; old=pair_products(coefficients,poly,guard,counters)
    correction,denom,report=bucket_products(rows,geometry,grid,poly,guard,lambda *a,**k:None)
    assert denom==common
    expected=defaultdict(lambda:[F(0)]*6); old_expected=defaultdict(lambda:[F(0)]*6)
    target_hi=defaultdict(lambda:[F(0)]*6)
    offdiag=zeros=ties=clipped=0
    for (da,wa,a),(db,wb,b) in product(atoms,repeat=2):
        index=da[0]+db[0]+STRIDE*(da[1]+db[1]); weight=wa*wb
        for channel,(ai,bi,ci,di) in enumerate(PATHS):
            aa,cc,bb,dd=a[ai],a[ci],b[bi],b[di]
            p,q=aa*bb,cc*dd; ideal=min(p,q)
            i,j=nearest_integer(aa,cc,boundaries),nearest_integer(dd,bb,boundaries)
            if i is None or j is None: lower=F(0); zeros+=1; assert ideal==0
            elif j<i: lower=F(p); offdiag+=1
            elif j>i: lower=F(q); offdiag+=1
            else:
                rho=grid[i]; lower=min(rho*aa,cc)*min(F(bb)/rho,dd)
            assert 0<=lower<=ideal
            for rho in grid:
                assert min(rho*aa,cc)*min(F(bb)/rho,dd)<=lower; checks+=1
            if aa and cc:
                ratio=F(cc,aa); ties+=any(ratio*ratio==x*y for x,y in zip(grid,grid[1:]))
                clipped+=ratio<grid[0] or ratio>grid[-1]
            expected[index][channel]+=weight*lower
            old_expected[index][channel]+=weight*(p+q)
            target_hi[index][channel]+=weight*base.sqrt_bounds(F(p*p+q*q))[1]
            # A conservative rational kappa lower bound for alpha1/2.
            assert (F(p+q)-F(1,2)*lower)**2>=p*p+q*q; checks+=1
    for index,row in expected.items():
        for j,value in enumerate(row):
            assert F(int(correction[j][index]),common)==value
            total=int(old[j][index]); assert total==old_expected[index][j] and 2*value<=total
            assert F(total)-F(1,2)*value>=target_hi[index][j]
            # alpha1 has kappa0: old packet sums and terminal order remain.
            assert F(total)-F(0)*value==old_expected[index][j]; checks+=4
    assert offdiag and zeros and ties and clipped
    # Upward powered freezing may alter ratios, but G is monotone. Check
    # a rational lower-kappa surrogate that is monotone as well.
    for a,b,c,d in product((F(0),F(1,3),F(5,4)),repeat=4):
        frozen=[F(ceil_dyadic(x,8),1<<8) for x in (a,b,c,d)]
        p,q=a*b,c*d; pp,qq=frozen[0]*frozen[1],frozen[2]*frozen[3]
        majorant=pp+qq-F(1,2)*min(pp,qq)
        assert majorant>=p+q-F(1,2)*min(p,q) and majorant*majorant>=p*p+q*q
        checks+=2
    return {"ordered_exact_checks":checks,"integer_offdiagonal_pairs":offdiag,
        "zero_pairs":zeros,"boundary_ties":ties,"clipped_left_ratios":clipped,
        "six_MM_MT_channels":True,"grid_witness_domination":True,
        "powered_freeze_monotone_control":True,"integer_bucket_report":report,
        "inherited_local_tests":base.exact_tests(v,guard)}


def build_local_pairs(z,alpha,anchor,guard,deadline,precision=256,progress=None):
    """Local-only six-channel upper coefficients, API-compatible with7284."""
    z,alpha=F(z),F(alpha)
    if not 0<z<1 or not 0<alpha<=1: raise ValueError("require0<z<1,0<alpha<=1")
    if precision not in (256,384): raise ValueError("precision256/384 only")
    from flint import ctx, fmpz_poly
    ctx.prec=max(384,precision)
    emit=(lambda phase,**kw:None) if progress is None else progress
    base=load("cheap_lane_pair_local_directed.py",LOCAL_PIN,"ordered_local_base")
    rank=load("cheap_lane_residual_rank_return.py",RANK_PIN,"ordered_local_rank")
    ordered=load("cheap_lane_rho_ordered_buckets.py",ORDERED_PIN,"ordered_local_proof_interface")
    grid=ordered.load_fine().GRID
    composition=base.load("cheap_lane_pair_concavity_triple.py",base.COMPOSITION_PIN,"ordered_local_composition")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"ordered_local_triple")
    caps=rank.augment_caps(composition.build_caps(z,alpha,triple,anchor,guard,deadline),z,alpha,guard)
    v=caps["verifier"]; power=v.exact(alpha)
    tests=exact_tests(base,v,fmpz_poly,guard)
    shapes=sorted(caps["raw"]); raw=[]; powered=[]; raw_bits=0
    max_ratio=F(1); zeros=nonzeros=0; raw_hash=hashlib.sha256(); powered_hash=hashlib.sha256()
    for occ in shapes:
        guard(); raw_row=tuple(base.endpoint(v.positive(x),True) for x in caps["raw"][occ])
        assert all(x>=0 for x in raw_row)
        raw.append(raw_row); raw_bits=max(raw_bits,*(denominator_bits(x) for x in raw_row))
        frozen=[]
        for value in caps["powered"][occ]:
            upper=base.endpoint(value,True); lower=base.endpoint(value,False)
            if upper<0: raise ArithmeticError("negative powered cap")
            if value.is_zero(): integer=0; zeros+=1
            else:
                if lower<=0: raise ArithmeticError("nonzero powered cap not separated from zero")
                integer=ceil_dyadic(upper); nonzeros+=1
                max_ratio=max(max_ratio,F(integer,1<<POWER_BITS)/lower)
            assert F(integer,1<<POWER_BITS)>=upper
            frozen.append(integer)
        powered.append(tuple(frozen))
        raw_hash.update((str(occ)+":"+":".join(map(str,raw_row))+"\n").encode())
        powered_hash.update((str(occ)+":"+":".join(map(str,frozen))+"\n").encode())
    raw_int=[tuple(x.numerator<<(raw_bits-denominator_bits(x)) for x in row) for row in raw]
    del caps["raw"],caps["powered"],raw
    geometry,route_report=base.route_geometry(shapes,guard)
    emit("ordered_coherent_tables_ready",occupancies=len(shapes),working_precision=ctx.prec,
         powered_freeze_bits=POWER_BITS,raw_common_denominator_bits=raw_bits,**route_report)
    counters={"pair_polynomial_products":0}
    old_packet=aggregate_rows(powered,geometry,fmpz_poly,guard)
    old_pair=pair_products(old_packet,fmpz_poly,guard,counters); del old_packet
    raw_packet=aggregate_rows(raw_int,geometry,fmpz_poly,guard)
    raw_pair=pair_products(raw_packet,fmpz_poly,guard,counters); del raw_packet,raw_int
    correction,diagonal_denominator,bucket_report=bucket_products(powered,geometry,grid,fmpz_poly,guard,emit)
    del powered,geometry
    old_denominator=1<<(2*POWER_BITS+64)
    raw_denominator=1<<(2*raw_bits+64)
    c_lower=max(F(0),base.endpoint(2-v.exact(F(2))**power,False))
    assert 0<=c_lower<=1
    rational_coefficients={}; coefficients={}; result_hash=hashlib.sha256()
    corrected_entries=jensen_entries=0
    for f in range(65):
        guard()
        for s in range(65-f):
            index=f+STRIDE*s; mass=comb(64,f)*comb(64-f,s)
            factor=v.exact(F(mass))**v.exact(1-alpha); row=[]
            for j in range(6):
                old_value=int(old_pair[j][index]); corr=int(correction[j][index]); raw_value=int(raw_pair[j][index])
                if min(old_value,corr,raw_value)<0 or 2*corr>diagonal_denominator*old_value:
                    raise ArithmeticError("invalid exact positive correction coefficient")
                old_cap=F(old_value,old_denominator)
                lower=F(corr,diagonal_denominator*old_denominator)
                corrected=old_cap-c_lower*lower
                jensen=F(0) if raw_value==0 else base.endpoint(factor*v.exact(F(raw_value,raw_denominator))**power,True)
                final=min(corrected,jensen)
                assert 0<=final<=old_cap
                corrected_entries+=corrected<old_cap; jensen_entries+=jensen<corrected
                row.append(final)
            rational_coefficients[f,s]=tuple(row); coefficients[f,s]=tuple(v.exact(x) for x in row)
            result_hash.update((str((f,s))+":"+":".join(map(str,row))+"\n").encode())
    guard()
    report={"scope":"directed local ordered binary-pair coefficient upper bounds only",
        "directed":True,"whole_code_certificate":False,"source_factor_applied":False,
        "support_count_applied":False,"regional_marker_applied":False,
        "z":str(z),"alpha":str(alpha),"requested_precision":precision,"working_precision":ctx.prec,
        "pair_channels":["H00","H01","H10","H11","J0","J1"],
        "regular_pair":"M*M","terminal_pair":"M*T","pair_slots":64,"packing_stride":STRIDE,
        "powered_freeze":{"common_absolute_bits":POWER_BITS,"exact_zeros":zeros,
            "positive_atoms":nonzeros,"maximum_ratio_vs_Arb_lower":str(max_ratio),
            "loss_assumed_zero":False,"table_sha256":powered_hash.hexdigest()},
        "raw_table_sha256":raw_hash.hexdigest(),"raw_common_denominator_bits":raw_bits,
        "rational_coefficient_sha256":result_hash.hexdigest(),"coefficient_degrees":len(coefficients),
        "correction_constant_lower":str(c_lower),"corrected_entries":corrected_entries,
        "Jensen_further_improved_entries":jensen_entries,"ordered_buckets":bucket_report,
        "integer_pair_products":counters["pair_polynomial_products"],"route_geometry":route_report,
        "tests":tests,"field_tests":caps["tests"],"residual_rank_receipt":caps["residual_rank_receipt"],
        "nonnegative_intersections":v.intersections,"local_pin":LOCAL_PIN,"rank_pin":RANK_PIN,
        "ordered_pin":ORDERED_PIN,"composition_pin":base.COMPOSITION_PIN,
        "triple_pin":composition.TRIPLE_PIN,"generated_bytes":0}
    return {"verifier":v,"core":caps["core"],"rational_coefficients":rational_coefficients,
            "coefficients":coefficients,"report":report}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(23,100))
    parser.add_argument("--alpha",type=F,default=F(23,256))
    parser.add_argument("--precision",type=int,choices=(256,384),default=256)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    args=parser.parse_args()
    if not 0<args.z<1 or not 0<args.alpha<=1: parser.error("require0<z<1,0<alpha<=1")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    anchor=load("cheap_lane_cone_cauchy_anchors.py",ANCHOR_PIN,"ordered_local_anchor")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("ordered local deadline; no retry")
        if peak()>MEMORY: raise MemoryError("512MiB ordered local budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)
    built=build_local_pairs(args.z,args.alpha,anchor,guard,deadline,args.precision,progress)
    guard()
    print(json.dumps({**built["report"],"seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
