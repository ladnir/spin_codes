"""Exact coarse-count domination checks. Serial90s,<128MiB,no files."""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_mixed_counts as mixed
from quartet_window_transfer import gate_counts,kernel
from quartet_two_cohort import compositions,elementary,H,D
from quartet_history_screen import occupation


def cell(a,r,s,width):
    lo=(a//2)*2;return (lo,min(width,lo+1),r//2,s//2)


def phi(n):return math.prod(h**v for h,v in zip(H,n))


@lru_cache(maxsize=512)
def old_upper(nlo,nhi,slo,shi):
    result=[]
    for degree in range(nhi+1):
        feasible=[(n,s) for n in range(max(nlo,degree),nhi+1)
                  for s in range(slo,min(shi,n)+1)]
        assert feasible
        candidates=[]
        for z in (F(1,128),F(1,32),F(1,8),F(1,2),F(2)):
            l2=F(1+D[2]*z,H[2]);old=max(F(1+D[i]*z,H[i]) for i in (3,4))
            candidates.append(max(l2**s*old**(n-s) for n,s in feasible)/z**degree)
        result.append(min(candidates))
    return result


def upper_matrix(width):
    As={}
    for a in range(width+1):
        allowed=defaultdict(F)
        for (l,k),J in gate_counts(width).items():
            for v0 in range(max(0,k-a),min(k,width-a)+1):
                allowed[l,k-v0]+=F(J*math.comb(width-a,v0),math.comb(width,k))*H[1]**(l-k+v0)
        As[a]=allowed
    M={}
    for alo in range(0,width+1,2):
        ahi=min(width,alo+1)
        for r in range(ahi+1):
            nlo=max(0,alo-r);nhi=ahi-r
            for slo in range(0,nhi+1,2):
                E=old_upper(nlo,nhi,slo,min(nhi,slo+1));row=defaultdict(F)
                for l in range(width+1):
                    for u in range(min(ahi,l)+1):
                        A=As[max(alo,u)][l,u]
                        for v in range(max(0,u-nhi),min(r,u)+1):
                            row[l,(l-u)//2,v//2]+=A*math.comb(r,v)*D[1]**v*E[u-v]/H[1]**r
                source=(alo,ahi,r//2,slo//2);existing=M.setdefault(source,defaultdict(F))
                for edge,value in row.items():existing[edge]=max(existing[edge],value)
    return M


def local_tests():
    checked=0
    for width in (4,8):
        M=upper_matrix(width)
        for n in compositions(width,5):
            actual=defaultdict(F)
            for (out,l,exponent),p in kernel(n):
                actual[l,out[1]//2,out[2]//2]+=p*(1<<exponent)*F(phi(out),phi(n))
            source=cell(width-n[0],n[1],n[2],width)
            assert all(value<=M[source][e] for e,value in actual.items());checked+=len(actual)
    print(dict(exact_coarse_count_edge_bounds=checked,all_small_histograms=True),flush=True)


def region_tests():
    width=4;M=upper_matrix(width)
    direct={((width,0,0,0,0),0,0):F(1)};coarse={}
    for step in range(4):
        after=defaultdict(F)
        for (n,x,y),value in direct.items():
            for (out,l,exponent),weight in kernel(n):
                after[out,x+(l if step<2 else 0),y+(l if step>=2 else 0)]+=value*weight*(1<<exponent)*F(phi(out),phi(n))
        direct=dict(after);after=defaultdict(F)
        if step==0:
            for (n,x,y),value in direct.items():after[cell(width-n[0],n[1],n[2],width),x,y]+=value
        else:
            for (old,x,y),value in coarse.items():
                for (l,rbin,sbin),weight in M[old].items():
                    new=cell(l,2*rbin,2*sbin,width)
                    after[new,x+(l if step<2 else 0),y+(l if step>=2 else 0)]+=value*weight
        coarse=dict(after)
    exact=defaultdict(F);upper=defaultdict(F)
    for (n,x,y),value in direct.items():exact[x,y]+=value*4**(width-n[0])
    for (c,x,y),value in coarse.items():upper[x,y]+=value*4**c[1]
    assert len(exact)==81 and all(value<=upper[key] for key,value in exact.items())
    print(dict(exact_coarse_count_two_region_pairs=len(exact),true_initial_row=True),flush=True)


def interval_coefficients():
    original_h=mixed.fresh.H.copy();original_d=mixed.fresh.D.copy();checked=0
    try:
        mixed.fresh.H=np.array(H,dtype=float);mixed.fresh.D=np.array(D,dtype=float)
        for lo in range(0,9,2):
            for hi in range(lo,9):
                bins=[(s,min(hi,s+1)) for s in range(0,hi+1,2)]
                upper=mixed.old_interval_coefficients(lo,hi,bins)
                for n in range(lo,hi+1):
                    for s,t,u in compositions(n,3):
                        coeff=elementary((s,t,u),D[2:]);den=H[2]**s*H[3]**t*H[4]**u
                        for k,c in enumerate(coeff):
                            assert math.log2(c)-math.log2(den)<=upper[s//2,k]+1e-10;checked+=1
    finally:mixed.fresh.H=original_h;mixed.fresh.D=original_d
    print(dict(exact_coefficients_vs_interval_float_bound=checked,directed_full_width=False),flush=True)


def layout_tests():
    intervals=[(a,min(512,a+63)) for a in range(0,513,64)]
    cells,sources,lookup,labels,target_map=mixed.layout(intervals)
    checked=0
    for l in range(513):
        for r in range(0,l+1,64):
            for s in range(0,l-r+1,64):
                e=lookup[l,r//64,s//64];c=cells[target_map[e]]
                assert labels[e]==l and c[0]<=l<=c[1]
                assert c[2]<=r<=c[3] and c[4]<=s<=c[5];checked+=1
    logs=np.array([-np.inf,-2000.,-2.,0.]);v=mixed.positive_exp(logs,'test')
    assert v[0]==0 and v[1]==2**-1000 and v[2]==.25 and v[3]==1
    mixed.FLOOR_COUNTS.clear()
    print(dict(full_width_exact_edge_labels_checked=checked,cells=len(cells),
               upward_floor_preserves_structural_zero=True),flush=True)


def coupled_tests():
    original_h=mixed.fresh.H.copy();original_d=mixed.fresh.D.copy();checked=vertices_checked=0
    try:
        mixed.fresh.H=np.array(H,dtype=float);mixed.fresh.D=np.array(D,dtype=float)
        for width in (4,8):
            cached={}
            for n in compositions(width,5):
                a=width-n[0];lo,hi,rb,sb=cell(a,n[1],n[2],width)
                c=(lo,hi,2*rb,min(hi,2*rb+1),2*sb,min(hi-2*rb,2*sb+1))
                if c not in cached:
                    cached[c]=mixed.coupled_coefficients(c,width)
                    all_points=[(aa,r,s) for aa in range(lo,hi+1)
                                for r in range(c[2],c[3]+1) for s in range(c[4],c[5]+1) if r+s<=aa]
                    verts=mixed.cell_vertices(c)
                    for x,y,z in ((1,2,3),(-1,2,-3),(4,-2,-1),(-3,-1,2)):
                        assert max(x*aa+y*r+z*s for aa,r,s in verts)==max(x*aa+y*r+z*s for aa,r,s in all_points)
                        vertices_checked+=1
                coeff=elementary(n,D);den=phi(n)
                for k,value in enumerate(coeff):
                    assert math.log2(value)-math.log2(den)<=cached[c][k]+1e-10;checked+=1
    finally:mixed.fresh.H=original_h;mixed.fresh.D=original_d
    print(dict(coupled_full_population_coefficients_checked=checked,exact_cell_vertex_objectives=vertices_checked),flush=True)
    import itertools
    checks=0
    for caps in ((3,4,2),(1,0,3),(2,2,2)):
        for total in range(9):
            for values in ((4,1,3),(0,2,2),(1,1,1)):
                bound=min(theta*total+sum(m*max(v-theta,0) for m,v in zip(caps,values))
                          for theta in (0,)+values)
                direct=max(sum(k*v for k,v in zip(ks,values)) for ks in itertools.product(*(range(m+1) for m in caps)) if sum(ks)<=total)
                assert bound==direct;checks+=1
    print(dict(exact_row_cap_linear_programs=checks),flush=True)
    # With k<=a_lo the independently maximized population sizes admit a
    # Vandermonde coefficient for575 coordinates instead of512.
    w=512;l=408;k=208;alo=320;ahi=383
    J=next(J for kk,J,logJ in occupation(l) if kk==k)
    numerator=J*math.comb(w+ahi-alo,k)
    denominator=math.comb(w,k)*math.comb(w,l)
    bits=numerator.bit_length()-denominator.bit_length()
    while F(numerator,denominator)<F(2)**bits:bits-=1
    assert bits>0 and k<=alo
    print(dict(exact_unweighted_coarse_factor_lower_bits=bits,
               incoming_count_interval=[alo,ahi],outgoing_count=l,selected_count=k,
               actual_population=w,phantom_population=w+ahi-alo),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        local_tests();region_tests();interval_coefficients();layout_tests();coupled_tests()
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
