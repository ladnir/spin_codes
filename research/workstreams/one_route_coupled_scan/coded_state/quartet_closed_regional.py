"""Bounded closed-set population transfer diagnostic, not a certificate.

Preserves original inactive labels and carries population across regions.
Maximizes run ages only inside each local coefficient. FLOAT proposals.
Serial Peach: 240 seconds, <=512 MiB, no generated files.
"""
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_closed_cells import bounds,vertices
from quartet_closed_count import coefficient_caps
from quartet_closed_transfer import TYPES,terminal


def cells(width=512,step=64):
    return [(a,min(width,a+step-1),s,min(width,s+step-1))
            for a in range(0,width+1,step) for s in range(0,a+1,step)]


def destination_cost(cell,logR,logS,inactive,width=512):
    # Every disallowed input was originally inactive. Hence a >= width-i.
    alo,ahi,slo,shi=cell
    if ahi<width-inactive:return math.inf
    clipped=(max(alo,width-inactive),ahi,slo,shi)
    ns=vertices(clipped,width)
    return min(logR*(width-n[0])+logS*n[2] for n in ns)


def initial(logh,logb,logR,logS,width=512):
    """Exact marked coefficient for the dummy all-zero incoming state."""
    h=dict(zip(TYPES,logh));values=[]
    active=logR+np.logaddexp2(h[1],h[2]+logS+707*logb/2560)
    inactive=np.logaddexp2(0.,h[2]+logR+logS+(707/2560-6/5)*logb)
    for i in range(width+1):
        values.append(math.log2(math.comb(width,i))+i*inactive+(width-i)*active)
    return np.array(values)


def terminal_caps(cs,logh,width=512):
    accepted={choice:[logh[k] for k,t in enumerate(TYPES)
                      if t%3==choice and terminal((t,))] for choice in range(3)}
    assert not accepted[1]
    out=[]
    for alo,ahi,slo,shi in cs:
        lo=max(alo,slo);hi=min(ahi,shi)
        out.append(max(-(width-a)*min(accepted[0])-a*min(accepted[2])
                       for a in (lo,hi)) if lo<=hi else -math.inf)
    return np.array(out)


def build(logh,logb=2.5,step=64):
    cs=cells(step=step);count=len(cs);start=time.monotonic()
    edges=np.full((count,513,count),np.inf);first=np.full((513,count),np.inf)
    for logR,logS in itertools.product((-.5,0.,.5),repeat=2):
        caps=coefficient_caps(logh=logh,logb=logb,logR=logR,logS=logS)
        costs=np.array([[destination_cost(c,logR,logS,i) for c in cs] for i in range(513)])
        # Impossible destinations retain -infinity, not a finite floor.
        first=np.minimum(first,initial(logh,logb,logR,logS)[:,None]-costs)
        for k,c in enumerate(cs):
            values=bounds(caps,c)
            edges[k]=np.minimum(edges[k],values[:,None]-costs)
        print(dict(marker=(logR,logS),seconds=time.monotonic()-start),flush=True)
    assert not np.any(np.isnan(edges))
    return cs,edges,first,terminal_caps(cs,logh)


def cycle_diagnostic(cs,edges,first,last,logb=2.5):
    """Lower contribution of this UPPER relaxation; NOT code failure mass.

    Each region has twelve 131 labels and four 132 labels, totaling2100.
    The first step enters the cell, then it persists through all4096 steps.
    One final terminal test. The positive relaxed transfer contains this
    term, so a positive value rules out this fixed proposal.
    """
    logden=math.log2(math.comb(8192,2100));out=[]
    for k,c in enumerate(cs):
        if not math.isfinite(last[k]):continue
        e131=edges[k,131,k];e132=edges[k,132,k]
        if not all(math.isfinite(v) for v in (e131,e132,first[131,k])):continue
        regional=12*e131+4*e132+math.log2(math.comb(16,4))
        # First inactive count131 is fixed; choose four 132 slots from15.
        leading=first[131,k]+11*e131+4*e132+math.log2(math.comb(15,4))
        value=leading+255*regional-255*logden-(713983/5-111127)*logb+last[k]
        out.append((float(value),c,float(regional-logden)))
    out.sort(reverse=True)
    print(dict(relaxation_cycle_lower=out[:5],logb=logb,
               inactive_quota=2100,updates_per_region=16,
               terminal_checked_once=True),flush=True)
    return out


def tests():
    # Exact initial polynomial has independent factors because the dummy
    # incoming state is zero. Compare against the quartet transfer itself.
    from quartet_closed_transfer import local
    for logR,logS in ((0.,0.),(-.5,.5),(.5,-.5)):
        h=np.zeros(len(TYPES));logb=2.5
        exact=np.zeros(5)
        for (dest,u,v,a,d),ways in local((0,0,0,0)):
            for paid in range(v+1):
                exact[u+paid]+=ways*math.comb(v,paid)*2**(
                    logR*sum(t%3!=0 for t in dest)+logS*a+logb*(707*a/2560-6*paid/5))
        assert np.max(abs(np.log2(exact)-initial(h,logb,logR,logS,4)))<1e-10
    checked=0
    for c in cells(8,2):
        for i in range(9):
            for r,s in itertools.product((-.5,0.,.5),repeat=2):
                cost=destination_cost(c,r,s,i,8)
                pairs=[(a,z) for a in range(c[0],c[1]+1)
                       for z in range(c[2],min(c[3],a)+1) if a>=8-i]
                assert cost==min((r*a+s*z for a,z in pairs),default=math.inf)
                checked+=1
    print(dict(initial_coefficients_checked=15,destination_cell_checks=checked),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(240,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        tests()
        # Run-age-dependent potentials must not earn a fresh artificial
        # gain at every age-maximized step. Use no hidden-age potential.
        cs,edges,first,last=build(np.zeros(len(TYPES)))
        cycle_diagnostic(cs,edges,first,last)
    finally:timer.cancel()


if __name__=='__main__':main()
