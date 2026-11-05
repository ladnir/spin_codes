"""FLOAT full-weight routing screen with actual uncovered coordinates.

Fixed construction, no files generated. Serial, 60-second watchdog,
>=20 GiB free, <16 MiB core numeric arrays for the default two kernels;
optional marked flow diagnostics have a 64 MiB numeric-array budget.
No whole-code certificate is claimed.
"""
import argparse
from fractions import Fraction as F
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np

import quartet_full_weight_kernel as full
import quartet_quiet_regions as quiet


def free_cap(p,r,m=8192):
    """Exact combinatorial expression, evaluated in FLOAT for each free count.

    Labels on uncovered positions do not enter the local integrand. For E
    such positions, maximize their remaining Bin(E,p) mass over feasible
    covered-label counts before dividing by the complete quota probability.
    """
    pi=quiet.log_choose(m,r)+r*math.log2(p)+(m-r)*math.log2(1-p)
    values=[]
    for e in range(m+1):
        lo=max(0,r-(m-e));hi=min(e,r)
        mode=min(hi,max(lo,math.floor((e+1)*p)))
        values.append(quiet.log_choose(e,mode)+mode*math.log2(p)
                      +(e-mode)*math.log2(1-p)-pi)
    return np.array(values)


def affine_cap(caps,marker):
    offsets=caps-np.arange(len(caps))*math.log2(marker)
    return float(max(offsets)),int(np.argmax(offsets))


def apply_region(kernel,vector,constant,length=16):
    for _ in range(length):
        vector=np.logaddexp2.reduce(kernel+vector[None,:],axis=1)
    return vector+constant


def screen(records,r,theta,ell,X,coverage):
    terminal=-np.arange(513)*(math.log2(X)+theta*math.log2(ell))
    potential=terminal.copy();best=math.inf
    base=(quiet.log_choose(8192,r)+theta+32*theta*(943718-1048576+128*r)
          +theta*943718*math.log2(ell)-129*(8192-r)*math.log2(coverage))
    for iteration in range(6):
        candidates=np.array([apply_region(k,potential,c) for k,c in records])
        applied=candidates.min(axis=0)
        growth=float(max(applied-potential))
        value=base+256*growth+potential[0]+float(max(terminal-potential))
        best=min(best,value)
        print(dict(iteration=iteration,log2_bound=float(value),best=float(best),
                   all_weights=True,all_deaths=True,free_terminal=True,
                   directed=False,whole_code=False),flush=True)
        potential=applied-applied[0]
    return best


def flow_diagnostic(kernel,moments):
    """Perron-weighted artificial histories, not the encoder distribution."""
    right=np.zeros(513);left=np.zeros(513)
    for _ in range(80):
        right=np.logaddexp2.reduce(kernel+right[None,:],axis=1);right-=max(right)
        left=np.logaddexp2.reduce(kernel+left[:,None],axis=0);left-=max(left)
    row=np.logaddexp2.reduce(kernel+right[None,:],axis=1)
    mass=left+right;mass-=np.logaddexp2.reduce(mass);mass=np.exp2(mass)
    # Row-normalization remains a probability even if the iterations have
    # not converged; report stationarity residual instead of assuming it.
    flow=mass[:,None]*np.exp2(kernel+right[None,:]-row[:,None])
    residual=float(np.max(np.abs(flow.sum(axis=0)-mass)))
    top=np.argsort(mass)[-5:][::-1]
    values={name:float(np.sum(flow*moment)) for name,moment in moments.items()}
    values.update(mean_state=float(mass@np.arange(513)),stationarity_residual=residual,
                  largest_weights=[(int(w),float(mass[w])) for w in top],
                  artificial_distribution_only=True)
    print(values,flush=True)


def tests():
    cells=0
    for m in (4,8,12):
        for r in range(m+1):
            p=F(2,5);pi=F(math.comb(m,r))*p**r*(1-p)**(m-r)
            caps=free_cap(float(p),r,m)
            for e in range(m+1):
                exact=max(F(math.comb(e,r-s))*p**(r-s)*(1-p)**(e-r+s)/pi
                          for s in range(max(0,r-e),min(r,m-e)+1))
                assert abs(caps[e]-math.log2(float(exact)))<1e-12
                cells+=1
            for marker in (.4,.7,1.):
                constant,_=affine_cap(caps,marker)
                assert np.all(constant+np.arange(m+1)*math.log2(marker)+1e-12>=caps)
    theta=.175;ell=.00377;X=.226;coverage=1.1;marker=.7;p=.67
    kernel=full.power(full.local(p,theta,ell,X,1,coverage,marker),2)
    full.boundary_checks(kernel,p,theta,ell,X,1,2,coverage,marker)
    # Summing every inactive label on a free position must give exactly one.
    terms=full.local(p,theta,ell,X,1,coverage,marker)
    for v in range(5):
        expected=math.comb(4,v)*((1-p)*X*coverage)**v*marker**(4-v)
        assert abs(terms[0,v]/expected-1)<1e-12
    # Marking incoming weight gives a deterministic total after coefficient
    # extraction. This tests the marked-polynomial derivative independently.
    _,moments=full.power(terms,2,{'incoming':{(a,v):a*x for (a,v),x in terms.items() if a}})
    assert np.max(np.abs(moments['incoming']-np.arange(9)[:,None]))<1e-11
    print(dict(exact_free_caps=cells,empty_quartet_checks=5),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--rows',type=int,default=2100)
        parser.add_argument('--p',type=float,default=112627/167924)
        parser.add_argument('--markers',type=float,nargs='+',default=[.6,.9])
        parser.add_argument('--ell',type=float,default=.00377)
        parser.add_argument('--coverage',type=float,default=1.)
        parser.add_argument('--flow',action='store_true')
        args=parser.parse_args()
        assert 0<args.p<1 and 0<=args.rows<=8192 and args.coverage>=1
        assert (1<<32)*args.ell>=1 and 1<=len(args.markers)<=3
        assert all(0<x<=1 for x in args.markers)
        tests()
        if args.tests_only:return
        theta=7/40;X=97905/432928;records=[]
        caps=free_cap(args.p,args.rows)
        for marker in args.markers:
            constant,touch=affine_cap(caps,marker)
            terms=full.local(args.p,theta,args.ell,X,1,args.coverage,marker)
            marks=({name:full.local(args.p,theta,args.ell,X,1,args.coverage,marker,name)
                    for name in ('free','deaths','active_coverage')} if args.flow else None)
            result=full.power(terms,128,marks)
            kernel=result[0] if args.flow else result
            full.boundary_checks(kernel,args.p,theta,args.ell,X,1,128,args.coverage,marker)
            records.append((kernel,constant))
            print(dict(marker=marker,affine_constant=constant,touch_free_count=touch,
                       seconds=time.monotonic()-start),flush=True)
            if args.flow:flow_diagnostic(kernel,result[1])
        screen(records,args.rows,theta,args.ell,X,args.coverage)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
