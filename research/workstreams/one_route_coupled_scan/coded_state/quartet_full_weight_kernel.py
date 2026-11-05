"""All-weight/all-turnover one-update kernel, FLOAT research diagnostic.

Two coefficient axes suffice once deaths are summed with their actual proof
weight. No state-band restriction, no turnover cutoff. This is a local
counting kernel, NOT an end-to-end certificate. Serial, <16 MiB core arrays
without marks; optional marked derivatives use a 32 MiB numeric-array budget.
60-second watchdog, >=20 GiB free, no generated files.
"""
import argparse
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np

import quartet_two_update as q
import quartet_quiet_regions as quiet


def local(p,theta,ell,X,score,coverage=1,free_marker=1,statistic=None,death_marker=1,covered_marker=1):
    terms={};logsigma=theta*(32+math.log2(ell))
    for a in range(5):
        before=(1<<a)-1
        for after in range(16):
            v=after.bit_count();d=(before&~after).bit_count()
            for target in range(16):
                coefficients=q.choices(4,before,after,target,True)
                if not any(coefficients):continue
                t=target.bit_count();covered=t if a else (after&target).bit_count()
                poly=sum(n*ell**j for j,n in enumerate(coefficients))
                value=math.comb(4,a)*p**t*(1-p)**(4-t)*(poly*2**(-32*covered))**theta
                value*=X**(v-a)*ell**(-theta*a)*2**(-logsigma*d)*score**(a+d)
                value*=coverage**(4-t if a else v)
                value*=free_marker**(4-v if a==0 else 0)
                value*=death_marker**d
                value*=covered_marker**covered
                if statistic=='free':value*=4-v if a==0 else 0
                elif statistic=='deaths':value*=d
                elif statistic=='active_coverage':value*=4-t if a else v
                elif statistic=='inactive_labels':value*=t
                elif statistic=='covered_inactive':value*=covered
                elif statistic is not None:raise ValueError(statistic)
                key=a,v;terms[key]=terms.get(key,0.)+value
    return {key:value for key,value in terms.items() if value>0}


def power(terms,groups,marks=None):
    size=4*groups+1;assert size*size*8*3<16*(1<<20)
    if marks is not None:
        assert len(marks)<=4 and size*size*8*14<32*(1<<20)
    current=np.full((size,size),-np.inf);current[0,0]=0
    for step in range(groups):
        nxt=np.full_like(current,-np.inf);length=4*step+1
        for (a,v),value in terms.items():
            nxt[a:a+length,v:v+length]=np.logaddexp2(nxt[a:a+length,v:v+length],
                current[:length,:length]+math.log2(value))
        if step==groups-1 and marks is not None:previous=current
        current=nxt
    moments={}
    if marks is not None:
        length=4*(groups-1)+1
        for name,marked in marks.items():
            derivative=np.full_like(current,-np.inf)
            for (a,v),value in marked.items():
                derivative[a:a+length,v:v+length]=np.logaddexp2(
                    derivative[a:a+length,v:v+length],
                    previous[:length,:length]+math.log2(groups*value))
            moment=np.zeros_like(current);valid=np.isfinite(current)
            moment[valid]=np.exp2(derivative[valid]-current[valid])
            moments[name]=moment
    current-=np.array([quiet.log_choose(4*groups,u) for u in range(size)])[:,None]
    return (current,moments) if marks is not None else current


def boundary_checks(kernel,p,theta,ell,X,score,groups,coverage=1,free_marker=1,death_marker=1,covered_marker=1):
    size=4*groups;weights=np.arange(size+1)
    birth=np.array([quiet.log_choose(size,v)+v*math.log2((1-p)*X*coverage)
                    +(size-v)*math.log2(free_marker) for v in range(size+1)])
    assert np.max(np.abs(kernel[0]-birth))<1e-8
    reset=[free_marker**4]+[sum(math.comb(4,t)*p**t*(1-p)**(4-t)*2**(-32*theta*t)*coverage**(4-t)*covered_marker**t
                    for t in range(a)) for a in range(1,5)]
    death=quiet.normalized_power(reset,groups)
    death+=weights*(-math.log2(X)-32*theta-2*theta*math.log2(ell)+2*math.log2(score))
    death+=weights*math.log2(death_marker)
    assert np.max(np.abs(kernel[:,0]-death))<1e-8
    assert abs(kernel[0,0]-size*math.log2(free_marker))<1e-10
    assert np.all(np.isfinite(kernel))
    print(dict(activation_cells=size+1,reset_cells=size+1,all_states_finite=True),flush=True)


def regional_apply(kernel,vector,p,r,length=16,chunk=32):
    """Apply the full weight-state regional envelope to a positive vector.

    Keep nonquiet count until imposing its regional quota cap. Chunking the
    incoming states bounds temporary memory; no dense regional matrix power.
    """
    n=len(vector);nonquiet=kernel.copy();nonquiet[0,0]=-np.inf
    current=vector[None,:].copy()
    for step in range(length):
        nxt=np.full((step+2,n),-np.inf)
        for lo in range(0,n,chunk):
            hi=min(n,lo+chunk)
            nxt[1:,lo:hi]=np.logaddexp2.reduce(
                nonquiet[None,lo:hi,:]+current[:,None,:],axis=2)
        nxt[:-1,0]=np.logaddexp2(nxt[:-1,0],current[:,0]+kernel[0,0])
        current=nxt
    cap=np.array([quiet.quota_cap(p,r,h) for h in range(length+1)])
    return np.logaddexp2.reduce(current+cap[:,None],axis=0)


def full_history_screen(kernels,r,theta,ell,X,score,coverage=1):
    assert score==1  # An activity tilt needs an explicit score-class factor.
    terminal=-np.arange(513)*(math.log2(X)+theta*math.log2(ell))
    potential=terminal.copy();best=math.inf
    base=(quiet.log_choose(8192,r)+theta+32*theta*(943718-1048576+128*r)
          +theta*943718*math.log2(ell)-129*(8192-r)*math.log2(coverage))
    # Diagnose the relaxed count on constant-weight histories. This counts
    # all turnovers, unlike the earlier 224/32 profile. Positive values are
    # failed upper bounds, not witnesses or lower bounds on code failure.
    fixed=np.min(np.array([base+kernel[0,1:]+4095*np.diag(kernel)[1:]
                          +terminal[1:]+256*quiet.quota_cap(p,r,16)
                          for p,kernel in kernels.items()]),axis=0)
    hardest=int(np.argmax(fixed))+1
    print(dict(constant_weight_screen=True,hardest_weight=hardest,
               log2_bound=float(fixed[hardest-1]),weight224_bound=float(fixed[223]),
               actual_code_obstruction=False),flush=True)
    for iteration in range(6):
        alternatives=np.array([regional_apply(kernel,potential,p,r) for p,kernel in kernels.items()])
        applied=alternatives.min(axis=0)
        residual=applied-potential;growth=float(max(residual))
        bound=base+256*growth+potential[0]+float(max(terminal-potential))
        best=min(best,bound)
        print(dict(inactive_rows=r,iteration=iteration,regional_log_growth=growth,
                   worst_potential_weight=int(np.argmax(residual)),log2_screen=float(bound),
                   coverage_marker=coverage,proof_p_choices=list(kernels),
                   all_weights=True,all_turnovers=True,all_quiet_lengths=True,
                   directed=False,whole_code=False),flush=True)
        potential=applied-applied[0]
    return best


def coarse_region_screen(kernels,r,theta,ell,X,coverage,width):
    """All-weight upper envelope, retaining regional endpoints before min-p.

    Each bin contains actual weights, not representative weights. Maximize
    over every incoming weight only after summing each outgoing bin.
    """
    from quartet_weight_flow import region_coefficients
    bins=[np.array([0])]+[np.arange(lo,min(lo+width,513)) for lo in range(1,513,width)]
    n=len(bins);envelope=np.full((17,n,n),np.inf)
    for p,kernel in kernels.items():
        coarse=np.empty((n,n))
        for i,rows in enumerate(bins):
            for j,cols in enumerate(bins):
                coarse[i,j]=np.max(np.logaddexp2.reduce(kernel[np.ix_(rows,cols)],axis=1))
        graded=region_coefficients(coarse)
        cap=np.array([quiet.quota_cap(p,r,h) for h in range(17)])[:,None,None]
        envelope=np.minimum(envelope,graded+cap)
    regional=np.logaddexp2.reduce(envelope,axis=0)
    terminal=np.array([max(-b*(math.log2(X)+theta*math.log2(ell))) for b in bins])
    for _ in range(256):
        terminal=np.logaddexp2.reduce(regional+terminal[None,:],axis=1)
    base=(quiet.log_choose(8192,r)+theta+32*theta*(943718-1048576+128*r)
          +theta*943718*math.log2(ell)-129*(8192-r)*math.log2(coverage))
    result=float(base+terminal[0])
    print(dict(coarse_regional_bins=n,width=width,log2_bound=result,
               incoming_max_not_representatives=True,all_weights=True,
               all_turnovers=True,all_quiet_lengths=True,directed=False,whole_code=False),flush=True)
    return result


def tests():
    from fractions import Fraction as F
    from itertools import product
    terms={(0,0):F(1,2),(1,3):F(1,3),(3,1):F(1,7),(4,4):F(1,11)}
    expected={}
    for (a,v),x in terms.items():
        for (b,w),y in terms.items():
            key=a+b,v+w;expected[key]=expected.get(key,F(0))+x*y
    actual=power(terms,2)
    for (u,v),value in expected.items():
        assert abs(actual[u,v]-math.log2(float(value/math.comb(8,u))))<1e-12
    print(dict(exact_small_power_cells=len(expected)),flush=True)

    # Independent path enumeration: retain quiet/nonquiet degree and apply
    # the quota cap only after finishing the complete region.
    raw=np.array([[1.,.2,.3],[.4,.5,.1],[.2,.3,.6]])
    terminal=np.array([1.,2.,3.]);p=.3;r=2100
    for length in (1,2,3):
        exact=np.zeros(3)
        for start in range(3):
            for suffix in product(range(3),repeat=length):
                path=(start,)+suffix;value=terminal[path[-1]];h=0
                for a,b in zip(path,path[1:]):
                    value*=raw[a,b];h+=bool(a or b)
                exact[start]+=value*2**quiet.quota_cap(p,r,h)
        got=regional_apply(np.log2(raw),np.log2(terminal),p,r,length,2)
        assert np.max(np.abs(got-np.log2(exact)))<1e-10

    # X is only a diagonal change of potential; it must not alter a path
    # count after the endpoint factors are restored.
    theta=.175;ell=.00377;gamma=1.5
    k1=power(local(p,theta,ell,.2,1,gamma),2)
    k2=power(local(p,theta,ell,.3,1,gamma),2)
    u=np.arange(9)[:,None];v=np.arange(9)[None,:]
    assert np.max(np.abs(k2-k1-(v-u)*math.log2(1.5)))<1e-11
    boundary_checks(k1,p,theta,ell,.2,1,2,gamma)
    plain=power(local(p,theta,ell,.2,1,1),2)
    assert np.all(k1+1e-12>=plain)
    print(dict(regional_path_enumerations=117,potential_conjugacy_cells=81,
               coverage_monotonicity_cells=81),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser();parser.add_argument('--compare-band',action='store_true')
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--score',type=float,default=1.)
        parser.add_argument('--regional',action='store_true')
        parser.add_argument('--coarse-regional',type=int,default=0)
        parser.add_argument('--rows',type=int,default=2100)
        parser.add_argument('--coverage',type=float,default=1.)
        parser.add_argument('--adaptive-p',action='store_true')
        parser.add_argument('--extra-p',type=float,action='append',default=[])
        parser.add_argument('--ell',type=float,default=13592/664659)
        parser.add_argument('--neutral-potential',action='store_true')
        args=parser.parse_args();assert args.score>0
        assert args.coverage>=1 and args.ell*(1<<32)>=1
        p=112627/167924;theta=7/40;ell=args.ell;X=97905/432928
        if args.neutral_potential:X=ell**(-theta)
        tests()
        if args.tests_only:return
        kernel=power(local(p,theta,ell,X,args.score,args.coverage),128)
        boundary_checks(kernel,p,theta,ell,X,args.score,128,args.coverage)
        print(dict(phase='all-weight kernel',seconds=time.monotonic()-start,
                   incoming_weights=513,outgoing_weights=513,all_deaths=True,
                   directed=False,whole_code=False),flush=True)
        if args.compare_band:
            assert args.score==args.coverage==1 and ell==13592/664659 and not args.neutral_potential
            import quartet_weight_flow as flow
            restricted=flow.step_matrix(flow.floating_box(p),0)[1:,1:]
            full=kernel[216:233,216:233]
            assert np.all(full+1e-8>=restricted)
            print(dict(three_index_comparison_cells=289,
                       largest_log2_difference=float(np.max(full-restricted)),
                       restricted_matrix_dominated=True),flush=True)
        if args.regional or args.coarse_regional:
            assert not args.compare_band and 0<=args.rows<=8192
            kernels={p:kernel}
            alternatives=args.extra_p+([args.rows/8192] if args.adaptive_p else [])
            for alternative_p in alternatives:
                assert 0<alternative_p<1
                kernels[alternative_p]=power(local(alternative_p,theta,ell,X,args.score,args.coverage),128)
                boundary_checks(kernels[alternative_p],alternative_p,theta,ell,X,args.score,128,args.coverage)
                print(dict(phase='adaptive-route kernel',seconds=time.monotonic()-start),flush=True)
            if args.coarse_regional:
                assert args.score==1 and 16<=args.coarse_regional<=512
                coarse_region_screen(kernels,args.rows,theta,ell,X,args.coverage,args.coarse_regional)
            else:
                full_history_screen(kernels,args.rows,theta,ell,X,args.score,args.coverage)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
