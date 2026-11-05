"""Small exact and numerical checks of capped regional composition."""
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
import quartet_capped_regional as capped
from quartet_mixed_count_tests import upper_matrix,cell,phi
from quartet_window_transfer import kernel
from quartet_two_cohort import compositions


def lp(m,g,f):
    return min(theta*g+sum(cap*max(value-theta,0) for cap,value in zip(m,f))
               for theta in [F(0)]+list(f))


def lp_tests():
    checked=0
    for beta in ((F(1),F(1,2),F(1,4)),(F(1,8),F(1,4),F(1,2)),(F(0),F(1),F(0))):
        values=list(itertools.product((F(0),F(1,4),F(1)),repeat=3))
        actual=capped.cap_values(np.array([beta],dtype=float),np.array(values,dtype=float))[:,0]
        for f,value in zip(values,actual):
            assert value==float(lp(beta,F(1),f));checked+=1
    print(dict(exact_dual_cap_evaluations=checked),flush=True)


def count_potential_tests():
    """An activity potential changes local weights but telescopes exactly."""
    from quartet_two_cohort import H
    checked=0
    for width in (4,8):
        histograms=compositions(width,5) if width==4 else ((2,2,2,1,1),(0,0,0,0,8),(8,0,0,0,0))
        for n in histograms:
            a=width-n[0]
            for out,l,exponent in (key for key,value in kernel(n)):
                for kappa in (F(1,2),F(2)):
                    tilted=math.prod((F(h)*(kappa if i else 1))**v
                                     for i,(h,v) in enumerate(zip(H,out)))
                    tilted/=math.prod((F(h)*(kappa if i else 1))**v
                                      for i,(h,v) in enumerate(zip(H,n)))
                    assert tilted==F(phi(out),phi(n))*kappa**(l-a)
                    checked+=1
    # Any sequence has only the initial/final potential factor remaining.
    for path in itertools.product(range(5),repeat=4):
        for kappa in (F(1,2),F(2)):
            a=0;value=F(1)
            for l in path:value*=kappa**(l-a);a=l
            assert value*kappa**(-a)==1
    print(dict(exact_activity_potential_edges=checked,telescoping_paths=1250),flush=True)


def setup(width):
    M=upper_matrix(width);cells=sorted(M);caps={c:defaultdict(F) for c in cells}
    for n in compositions(width,5):
        source=cell(width-n[0],n[1],n[2],width);row=defaultdict(F)
        for (out,l,exponent),p in kernel(n):row[l]+=p*(1<<exponent)*F(phi(out),phi(n))
        for l,value in row.items():caps[source][l]=max(caps[source][l],value)
    edges={l:sorted({e for row in M.values() for e in row if e[0]==l}) for l in range(width+1)}
    targets={l:[cell(l,2*r,2*s,width) for _,r,s in edges[l]] for l in edges}
    return M,cells,caps,edges,targets


def two_regions():
    width=4;M,cells,caps,edges,targets=setup(width)
    # Terminal upper factor for every actual histogram in a coarse cell.
    upper={(0,0):{c:F(4)**c[1] for c in cells}}
    for step in range(3,-1,-1):
        after=defaultdict(lambda:defaultdict(F))
        for (x,y),values in upper.items():
            for l in edges:
                key=(x+(l if step<2 else 0),y+(l if step>=2 else 0))
                f=[values[c] for c in targets[l]]
                for c in cells:
                    after[key][c]+=lp([M[c][e] for e in edges[l]],caps[c][l],f)
        upper=dict(after)
    direct={((width,0,0,0,0),0,0):F(1)}
    for step in range(4):
        after=defaultdict(F)
        for (n,x,y),value in direct.items():
            for (out,l,exponent),p in kernel(n):
                after[out,x+(l if step<2 else 0),y+(l if step>=2 else 0)]+=value*p*(1<<exponent)*F(phi(out),phi(n))
        direct=dict(after)
    exact=defaultdict(F)
    for (n,x,y),value in direct.items():exact[x,y]+=value*4**(width-n[0])
    initial=cell(0,0,0,width)
    assert len(exact)==81 and all(value<=upper[key][initial] for key,value in exact.items())
    print(dict(exact_capped_two_region_quota_pairs=81,inherited_cells=True,free_terminal=True),flush=True)


def implementation_counter():
    width=4;M,cells,caps,edges,targets=setup(width);index={c:i for i,c in enumerate(cells)}
    parts=[]
    for l in edges:
        g=np.array([float(caps[c][l]) for c in cells]);m=np.array([[float(M[c][e]) for e in edges[l]] for c in cells])
        parts.append(dict(l=l,targets=np.array([index[c] for c in targets[l]]),
                          beta=np.minimum(1,m/g[:,None]),cap=np.log2(g)))
    logrho=math.ceil(math.log2(max(sum(caps[c].values()) for c in cells)))+1
    terminal=np.array([4.**(c[1]-width) for c in cells])
    actual=capped.coefficients(parts,np.zeros(len(cells)),logrho,0,steps=4,width=width,center=0,terminal=terminal)
    upper={0:{c:F(4)**c[1] for c in cells}}
    for _ in range(4):
        after=defaultdict(lambda:defaultdict(F))
        for total,values in upper.items():
            for l in edges:
                f=[values[c] for c in targets[l]]
                for c in cells:after[total+l][c]+=lp([M[c][e] for e in edges[l]],caps[c][l],f)
        upper=dict(after)
    scale=2**(4*logrho)*4**width;checked=0
    for total,row in upper.items():
        for c,value in row.items():
            numeric=actual[total,index[c]]*scale
            assert abs(numeric-float(value))<=1e-10*max(1,float(value));checked+=1
    print(dict(capped_counter_vs_exact_coefficients=checked),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        lp_tests();two_regions();implementation_counter();count_potential_tests()
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
