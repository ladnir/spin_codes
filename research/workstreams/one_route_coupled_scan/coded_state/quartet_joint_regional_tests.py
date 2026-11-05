"""Exact audit of uniform cohort buckets and quota-aware composition.

Small widths only. Serial90s,<64MiB,no files. No distance certificate.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import math
import os
from pathlib import Path
import shutil
import threading
import time
from quartet_window_transfer import kernel,gate_counts
from quartet_two_cohort import compositions,elementary,H,D

BIN=2


def label(n):return (sum(n)-n[0],n[1]//BIN,n[2]//BIN)


def phi(n):return math.prod(h**v for h,v in zip(H,n))


@lru_cache(maxsize=64)
def old_upper(count,lo,hi):
    result=[]
    for k in range(count+1):
        candidates=[]
        for z in (F(1,128),F(1,32),F(1,8),F(1,2),F(2)):
            second=F(1+D[2]*z,H[2]);old=max(F(1+D[i]*z,H[i]) for i in (3,4))
            candidates.append(max(second**s*old**(count-s) for s in (lo,hi))/z**k)
        result.append(min(candidates))
    return result


def matrix(width):
    M={}
    for a in range(width+1):
        allowed=defaultdict(F)
        for (l,k),J in gate_counts(width).items():
            for v0 in range(max(0,k-a),min(k,width-a)+1):
                allowed[l,k-v0]+=F(J*math.comb(width-a,v0)*D[0]**v0*H[1]**(l-k),math.comb(width,k))
        for r in range(a+1):
            count=a-r
            for slo in range(0,count+1,BIN):
                E=old_upper(count,slo,min(count,slo+BIN-1));row=defaultdict(F)
                for (l,u),A in allowed.items():
                    for v1 in range(max(0,u-count),min(r,u)+1):
                        row[l,(l-u)//BIN,v1//BIN]+=A*math.comb(r,v1)*D[1]**v1*E[u-v1]/H[1]**r
                source=(a,r//BIN,slo//BIN)
                existing=M.setdefault(source,defaultdict(F))
                for target,value in row.items():existing[target]=max(existing[target],value)
    return M


def local_tests():
    checked=0
    for width in (4,8):
        M=matrix(width)
        for n in compositions(width,5):
            direct=defaultdict(F)
            for (out,l,exponent),weight in kernel(n):
                direct[label(out)]+=weight*(1<<exponent)*F(phi(out),phi(n))
            assert all(value<=M[label(n)][target] for target,value in direct.items())
            checked+=len(direct)
    print(dict(exact_uniform_bucket_coefficients_checked=checked),flush=True)


def regional_tests():
    width=4;M=matrix(width)
    direct={((width,0,0,0,0),0,0):F(1)}
    coarse={((0,0,0),0,0):F(1)}
    for step in range(4):
        after=defaultdict(F)
        for (n,x,y),value in direct.items():
            for (out,l,exponent),weight in kernel(n):
                after[out,x+(l if step<2 else 0),y+(l if step>=2 else 0)]+=value*weight*(1<<exponent)*F(phi(out),phi(n))
        direct=dict(after);after=defaultdict(F)
        for (old,x,y),value in coarse.items():
            for new,weight in M[old].items():
                l=new[0]
                after[new,x+(l if step<2 else 0),y+(l if step>=2 else 0)]+=value*weight
        coarse=dict(after)
    exact=defaultdict(F);upper=defaultdict(F)
    for (n,x,y),value in direct.items():exact[x,y]+=value*4**(width-n[0])
    for (state,x,y),value in coarse.items():upper[x,y]+=value*4**state[0]
    assert len(exact)==81 and all(value<=upper[key] for key,value in exact.items())
    print(dict(exact_two_region_quota_pairs_checked=len(exact),inherited_cohorts=True,free_terminal=True),flush=True)


def numeric_coefficient_tests():
    import numpy as np
    import quartet_joint_regional as joint
    original_h=joint.fresh.H.copy();original_d=joint.fresh.D.copy();checked=0
    try:
        joint.fresh.H=np.array(H,dtype=float);joint.fresh.D=np.array(D,dtype=float)
        for count in range(9):
            bins=[(s,min(count,s+1)) for s in range(0,count+1,2)]
            upper=joint.old_coefficients(count,bins)
            for s,t,u in compositions(count,3):
                coeff=elementary((s,t,u),D[2:]);den=H[2]**s*H[3]**t*H[4]**u
                for k,c in enumerate(coeff):
                    assert math.log2(c)-math.log2(den)<=upper[s//2,k]+1e-11;checked+=1
    finally:joint.fresh.H=original_h;joint.fresh.D=original_d
    print(dict(exact_coefficients_vs_float_implementation=checked,directed_full_width=False),flush=True)


def geometry_recurrence_tests():
    checked=monotone=0
    for width in (4,8):
        for h in (F(1),F(3,2)):
            table={(l,u):F(J,math.comb(width,u))*h**(l-u)
                   for (l,u),J in gate_counts(width).items()}
            tables={}
            for a in range(width,-1,-1):
                tables[a]=table.copy()
                direct=defaultdict(F)
                for (l,k),J in gate_counts(width).items():
                    for v0 in range(max(0,k-a),min(k,width-a)+1):
                        direct[l,k-v0]+=F(J*math.comb(width-a,v0),math.comb(width,k))*h**(l-k+v0)
                for l in range(width+1):
                    for u in range(min(a,l)+1):
                        assert table.get((l,u),F(0))==direct[l,u];checked+=1
                table={(l,u):table.get((l,u),F(0))+h*table.get((l,u+1),F(0))
                       for l in range(width+1) for u in range(min(a-1,l)+1)}
            for alo in range(width+1):
                for a in range(alo,width+1):
                    for l in range(width+1):
                        for u in range(min(a,l)+1):
                            assert tables[a].get((l,u),F(0))<=tables[max(alo,u)].get((l,u),F(0))
                            monotone+=1
    print(dict(exact_inactive_coordinate_recurrence_checks=checked,
               exact_coarse_input_count_geometry_checks=monotone),flush=True)
    import numpy as np
    import quartet_joint_regional as joint
    counts=[0,1,127,392,408,424,511,512]
    actual=joint.allowed_tables(counts,counts);checked=0
    for a in counts:
        for l in counts:
            reference=joint.allowed_logs(a,l);value=actual[a][l]
            assert np.array_equal(np.isfinite(reference),np.isfinite(value))
            keep=np.isfinite(value)
            assert np.max(np.abs(reference[keep]-value[keep]))<1e-9;checked+=len(value)
    print(dict(float_geometry_recurrence_coefficients=checked,directed_full_width=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        local_tests();regional_tests();numeric_coefficient_tests();geometry_recurrence_tests()
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
