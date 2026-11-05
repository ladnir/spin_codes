"""Couple both outgoing age counts with the fixed population.

Local FLOAT diagnostic and exact toy identities; not a distance certificate.
Serial Peach,90s,<256MiB, no generated files.
"""
from collections import defaultdict
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_mixed_counts as mixed
import quartet_fresh_count_screen as fresh
from quartet_history_screen import occupation,logchoose
from quartet_window_transfer import kernel,gate_counts
from quartet_two_cohort import compositions,elementary,H,D


def geometry():
    return [(np.array([k for k,J,logJ in occupation(l)]),
             np.array([logJ-logchoose(512,k)+(l-k)*math.log2(fresh.H[1])
                       for k,J,logJ in occupation(l)])) for l in range(513)]


def moment_logs(cells,x,y,geom):
    """x,y are log2 markers for selected active positions and old type1."""
    d=fresh.D.copy();d[1:]*=2**x;d[1]*=2**y
    coeff=np.array([mixed.coupled_coefficients(c,d=d) for c in cells])
    result=np.empty((len(cells),513))
    for l,(ks,logg) in enumerate(geom):
        terms=coeff[:,ks]+logg[None,:];peak=terms.max(axis=1)
        result[:,l]=peak+np.log2(mixed.positive_exp(terms-peak[:,None],'joint_moment').sum(axis=1))
    return result


def edge_caps(cells,edge_ids,labels,geom,markers,progress=False):
    """Coefficientwise bound, minimizing only after uniform source maxima."""
    rlo=np.empty(len(labels),dtype=int);rhi=rlo.copy();slo=rlo.copy();shi=rlo.copy()
    for l in range(513):
        for r in range(0,l+1,mixed.BIN):
            for s in range(0,l-r+1,mixed.BIN):
                e=edge_ids[l,r//mixed.BIN,s//mixed.BIN]
                rlo[e]=r;rhi[e]=min(l,r+mixed.BIN-1)
                slo[e]=s;shi[e]=min(l-r,s+mixed.BIN-1)
    best=np.full((len(cells),len(labels)),np.inf)
    start=time.monotonic()
    for x,y in markers:
        total=moment_logs(cells,x,y,geom)
        endpoint=x*(rhi if x>=0 else rlo)-y*(slo if y>=0 else shi)
        bound=total[:,labels]-x*labels[None,:]+endpoint[None,:]
        best=np.minimum(best,bound)
        if progress:print(dict(stage='joint_moment_edges',x=x,y=y,seconds=time.monotonic()-start),flush=True)
    return best


def exact_tests():
    checks=0;envelopes=0;coefficient_checks=0
    markers=tuple(itertools.product((F(1,2),F(1),F(2)),repeat=2))
    for width in (4,8):
        histograms=list(compositions(width,5)) if width==4 else [(2,2,2,1,1),(0,0,0,0,8),(8,0,0,0,0)]
        for n in histograms:
            denominator=math.prod(h**v for h,v in zip(H,n))
            for x,y in markers:
                selected=[F(d)*(x if i else 1)*(y if i==1 else 1) for i,d in enumerate(D)]
                poly=elementary(n,selected)
                a=width-n[0]
                bound=mixed.coupled_coefficients((a,a,n[1],n[1],n[2],n[2]),width,
                                                 d=np.array(selected,dtype=float),h=np.array(H,dtype=float))
                for degree,value in enumerate(poly):
                    assert math.log2(value)-math.log2(denominator)<=bound[degree]+1e-10
                    coefficient_checks+=1
                direct=defaultdict(F);edges=defaultdict(F)
                for (out,l,exponent),p in kernel(n):
                    value=p*(1<<exponent)*F(math.prod(h**v for h,v in zip(H,out)),denominator)
                    direct[l]+=value*x**(l-out[1])*y**out[2]
                    edges[l,out[1]//2,out[2]//2]+=value
                formula=defaultdict(F)
                for (l,k),J in gate_counts(width).items():
                    formula[l]+=F(J,math.comb(width,k))*H[1]**(l-k)*poly[k]/denominator
                assert formula==direct;checks+=len(direct)
                for (l,r,s),value in edges.items():
                    rl=2*r;rh=min(l,rl+1);sl=2*s;sh=min(l-rl,sl+1)
                    endpoint=max(x**v for v in (rl,rh))*max(y**(-v) for v in (sl,sh))
                    assert value<=formula[l]*x**(-l)*endpoint;envelopes+=1
    print(dict(exact_joint_moment_identities=checks,exact_joint_edge_envelopes=envelopes,
               joint_coefficients_vs_float_envelope=coefficient_checks),flush=True)


def diagnostic():
    intervals=[(a,min(512,a+63)) for a in range(0,513,64)]
    selected=(320,383,0,63,0,63)
    base,scale,cells,edge_ids,labels,target_map=mixed.build(intervals,only_cell=selected)
    old=base[cells.index(selected)]
    with np.errstate(divide='ignore'):oldlogs=scale+np.log2(old)
    geom=geometry();markers=tuple(itertools.product((-1.,-.5,0.,.5,1.),repeat=2))
    new=edge_caps([selected],edge_ids,labels,geom,markers)[0]
    rows=[]
    for l in (384,392,408,424,448):
        take=labels==l;a=oldlogs[take];b=np.minimum(a,new[take]);peak=max(a);after=max(b)
        rows.append(dict(input_count=l,old_edge_sum_log2=float(peak+np.log2(np.exp2(a-peak).sum())),
                         joint_edge_sum_log2=float(after+np.log2(np.exp2(b-after).sum())),
                         improved_edges=int(np.sum(b<a-1e-8)),edges=int(sum(take))))
    print(dict(incoming_cell=selected,markers=len(markers),rows=rows,
               local_only=True,directed_certificate=False,whole_code_distance=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        exact_tests();diagnostic()
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
