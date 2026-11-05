"""Directed obstruction to a uniform one-step product-potential envelope.

This is NOT an obstruction to the encoder or its exact closed-set moment.
Fixes current field multipliers, but covers 1<=b<=8. Uses seven homogeneous
source rows and exact balanced flows to bound every positive product potential,
every positive original-inactive marker, and every paid-cap marker in(0,1].
No data written; serial <=180s.
"""
from fractions import Fraction as F
import math
import os
from pathlib import Path
import shutil
import threading
import time
from collections import defaultdict
from quartet_closed_quotient import kernel
from quartet_two_update import log_upper


def pure_atoms():
    rows=[]
    for source in range(7):
        row=defaultdict(int)
        for (dest,u,v,a,d),ways in kernel((source,)*4):
            counts=tuple(dest.count(i) for i in range(7))
            for paid in range(v+1):
                row[counts,u+paid,707*a-1905*d-3072*paid]+=ways*math.comb(v,paid)
        rows.append([(dest,i,e,c) for (dest,i,e),c in row.items()])
    return rows


def solve(matrix,rhs):
    a=[list(map(F,row))+[F(value)] for row,value in zip(matrix,rhs)]
    n=len(a)
    for col in range(n):
        pivot=next(j for j in range(col,n) if a[j][col])
        a[col],a[pivot]=a[pivot],a[col]
        scale=a[col][col];a[col]=[v/scale for v in a[col]]
        for j in range(n):
            if j!=col:
                scale=a[j][col];a[j]=[u-scale*v for u,v in zip(a[j],a[col])]
    return [row[-1] for row in a]


def witness(rows,logx):
    # Arbitrary FLOAT proposal. Exact normalized integer weights below
    # define the witness; optimizer accuracy is irrelevant to validity.
    logs=(0.,-.04475834141004524,-.042774336422858605,-.05588058545549474,
          -1.8326154305711793,-1.4212486763304475,-1.0098819495674625)
    probs=[];transition=[]
    for source,row in enumerate(rows):
        values=[math.log2(c)+e/1024+sum(n*h for n,h in zip(dest,logs))-4*logs[source]+i*logx
                for dest,i,e,c in row]
        peak=max(values);raw=[2**(v-peak) for v in values];total=sum(raw)
        counts=[max(1,round((1<<16)*v/total)) for v in raw]
        p=[F(c,sum(counts)) for c in counts];probs.append(p)
        transition.append([sum(prob*dest[k]/4 for prob,(dest,i,e,c) in zip(p,row)) for k in range(7)])
    matrix=[[transition[i][k]-int(i==k) for i in range(7)] for k in range(6)]+[[F(1)]*7]
    alpha=solve(matrix,[0]*6+[1])
    assert min(alpha)>0 and sum(alpha)==1
    assert all(sum(alpha[i]*transition[i][k] for i in range(7))==alpha[k] for k in range(7))
    inactive=sum(a*p*i for a,row,pr in zip(alpha,rows,probs) for p,(dest,i,e,c) in zip(pr,row))
    lower=F(0)
    for a,row,pr in zip(alpha,rows,probs):
        for p,(dest,i,e,c) in zip(pr,row):
            # log coefficient = log2(c)+e/1024 exactly. Signed helper
            # selects outward endpoints; subtract an upper for -log(c/p).
            lower+=a*p*(F(e,1024)-log_upper(-1,F(c)/p))
    paid=sum(a*p*(i-dest[0]) for a,row,pr in zip(alpha,rows,probs) for p,(dest,i,e,c) in zip(pr,row))
    energy=sum(a*p*e for a,row,pr in zip(alpha,rows,probs) for p,(dest,i,e,c) in zip(pr,row))
    return inactive,lower,len([p for pr in probs for p in pr]),paid,energy


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(180,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        rows=pure_atoms();first=witness(rows,-.3);second=witness(rows,.1)
        target=F(4*2100,8192)
        print(dict(witness_inactive_means=[float(first[0]),float(second[0])],target=float(target),
                   local_atoms=first[2],seconds=time.monotonic()-start),flush=True)
        assert first[0]<target<second[0]
        theta=(second[0]-target)/(second[0]-first[0]);assert 0<theta<1
        assert theta*first[0]+(1-theta)*second[0]==target
        lower=theta*first[1]+(1-theta)*second[1]
        paid=theta*first[3]+(1-theta)*second[3]
        energy=theta*first[4]+(1-theta)*second[4]
        intercept=(1<<19)*(lower-energy/1024)-log_upper(255,math.comb(8192,2100))
        slope=(1<<19)*energy/2560-F(79174)*2/5
        assert paid <= F(255463,1<<19)
        range_lower=min(intercept,intercept+3*slope)
        assert range_lower>0
        compact_range=F(range_lower.numerator*(1<<32)//range_lower.denominator,1<<32)
        print(dict(dual_paid_per_quartet=float(paid),paid_cap_per_quartet=255463/(1<<19),
                   affine_logb_intercept_lower=float(intercept),affine_logb_slope=float(slope),
                   all_b_between_1_and_8_lower=str(compact_range),
                   all_paid_cap_markers_between_0_and_1=True),flush=True)
        global_lower=(1<<19)*lower-79174-log_upper(255,math.comb(8192,2100))
        # Compact directed endpoint for the durable checkpoint.
        scale=1<<32;compact=F(global_lower.numerator*scale//global_lower.denominator,scale)
        print(dict(all_positive_product_potentials=True,all_positive_inactive_markers=True,
                   fixed_log2_b=2.5,uniform_one_step_envelope_only=True,
                   log2_global_expression_lower=str(compact),decimal_lower=float(compact),
                   excludes_mean_atom_penalty=True,excludes_nonnegative_terminal_penalty=True,
                   seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
