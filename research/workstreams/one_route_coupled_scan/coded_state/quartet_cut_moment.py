"""FLOAT ensemble bound for support-cut expansion; no distance certificate.

Serial,120s watchdog,64MiB numeric-array budget,20GiB free, no generated files.
Sums every state-support history, free terminal, and fixed-size row sets.
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
from quartet_regional_matrix import log_product
from quartet_remaining_activity import midpoint_envelope,bridge_weights


def local(p,b,reward,X,death_factor=None,canonical=False,label_marker=1):
    if death_factor is None:death_factor=1/b
    terms={}
    for a in range(5):
        before=(1<<a)-1
        for after in range(16):
            v=after.bit_count();d=(before&~after).bit_count()
            for target in range(16):
                if canonical:
                    slope,death_weight=canonical
                    if a==0 and target&after:continue
                    # Delete all a incoming selected states in this quartet.
                    # It frees at least these inactive inputs and increases
                    # deaths by at most the incoming/outgoing overlap h.
                    h=(before&after).bit_count();freed=(target&~after).bit_count()
                    if a and freed>=slope*a+death_weight*h:continue
                t=target.bit_count();c=t if a else (target&after).bit_count()
                value=math.comb(4,a)*p**t*(1-p)**(4-t)*label_marker**t*reward**v*b**(-c)*death_factor**d*X**(v-a)
                terms[a,v]=terms.get((a,v),0)+value
    return terms


def tests():
    p=F(2,5);b=F(16);reward=F(2);X=F(1,3)
    # Independent binomial formula separates overlap from inactive labels.
    checks=0
    for death_factor in (1/b,F(1,8)):
        exact=local(p,b,reward,X,death_factor)
        for a in range(5):
            for v in range(5):
                if a:
                    expected=math.comb(4,a)*(1-p+p/b)**4*reward**v*X**(v-a)*sum(
                        math.comb(a,j)*math.comb(4-a,v-j)*death_factor**(a-j)
                        for j in range(max(0,v-(4-a)),min(a,v)+1))
                else:expected=math.comb(4,v)*(reward*X*(1-p+p/b))**v
                assert exact[a,v]==expected;checks+=1
    print(dict(exact_local_cut_coefficients=checks),flush=True)
    filtered=local(p,b,reward,X,F(1,8),(F(1,4),F(3,4)))
    unfiltered=local(p,b,reward,X,F(1,8))
    assert all(0<filtered[key]<=unfiltered[key] for key in unfiltered)
    # Full deaths cannot free ANY inactive input at slope1/4, including
    # a=4: equality is excluded by minimum-cardinality tie breaking.
    for a in range(1,5):
        expected=math.comb(4,a)*(1-p)**4*F(1,8)**a*X**(-a)
        assert filtered[a,0]==expected
    for v in range(5):
        assert filtered[0,v]==math.comb(4,v)*((1-p)*reward*X)**v
    print(dict(canonical_boundary_coefficients=9,canonical_monotonic_cells=25),flush=True)
    gamma=F(4);z=1-p+p*gamma;tilted=p*gamma/z;tilt_checks=0
    for canonical in (False,(F(1,4),F(3,4))):
        marked=local(p,b,reward,X,F(1,8),canonical,label_marker=gamma)
        transformed=local(tilted,b,reward,X,F(1,8),canonical)
        for key in marked:
            assert marked[key]==z**4*transformed[key];tilt_checks+=1
    template_checks=0
    for u in range(5):
        for v in range(5):
            bits=(0,int(u>=2),int(v>=2))
            for c1 in range(5):
                for c2 in range(5):
                    counts=(c1,c2);factor=F(1)
                    for j in range(2):
                        old,new=bits[j:j+2]
                        factor*=gamma**(counts[j]*int(old or new)-new+3*old*(1-new))
                    activity=4*sum(bits[1:])
                    deaths=4*sum(old*(1-new) for old,new in zip(bits,bits[1:]))
                    coverage=sum(c*int(old or new) for c,old,new in zip(counts,bits,bits[1:]))
                    exponent=F(1,4)*activity-coverage-F(3,4)*deaths
                    assert factor==gamma**int(-exponent);template_checks+=1
    print(dict(exact_reference_label_tilt_cells=tilt_checks,
               exact_free_terminal_reference_paths=template_checks),flush=True)
    coarea_checks=0;convex_checks=0
    for u in range(5):
        for v in range(5):
            for c1 in range(5):
                for c2 in range(5):
                    total=F(0)
                    for threshold in range(1,5):
                        old=int(u>=threshold);new=int(v>=threshold)
                        total+=old+new-c1*old-c2*int(old or new)-3*old*(1-new)
                    averaged=F(u+v,4)-F(c1*u+c2*max(u,v),4)-F(3,4)*max(u-v,0)
                    assert total/4==averaged;coarea_checks+=1
    for canonical in (False,(F(1,4),F(3,4))):
        anchors=[local(p,b,reward,X,F(1,8),canonical,label_marker=F(marker)) for marker in (1,2,4,8,16)]
        for key in anchors[0]:
            assert anchors[1][key]**2<=anchors[0][key]*anchors[2][key]
            assert anchors[3][key]**2<=anchors[2][key]*anchors[4][key]
            convex_checks+=2
    print(dict(exact_averaged_threshold_paths=coarea_checks,
               exact_log_convexity_cells=convex_checks),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(120,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--logb',type=float,default=8.)
        parser.add_argument('--slope',type=float,default=.24)
        parser.add_argument('--death-weight',type=float,default=1.)
        parser.add_argument('--rows',type=int,default=2100)
        parser.add_argument('--p-levels',type=int,default=9)
        parser.add_argument('--flow',action='store_true')
        parser.add_argument('--canonical',action='store_true')
        parser.add_argument('--reference-burst',type=float,default=1.)
        parser.add_argument('--layercake',action='store_true')
        args=parser.parse_args();tests()
        if args.tests_only:return
        assert 0<args.slope<1 and 1<=args.logb<=16 and 0<args.rows<8192 and 2<=args.p_levels<=9
        assert 0<args.death_weight<=1
        assert not args.canonical or F(str(args.slope))+F(str(args.death_weight))<=1
        assert args.reference_burst>=1
        if args.reference_burst>1:
            assert args.rows==2100 and args.slope==.25 and args.death_weight==.74
        if args.layercake:
            assert args.reference_burst==1 and args.rows==2100 and args.slope==.25 and args.death_weight==.74
        r=args.rows;rho=r/8192;odds=rho/(1-rho);b=2**args.logb;reward=2**(args.slope*args.logb);X=.25
        records=[];regional=np.full((513,513),np.inf)
        canonical=(F(str(args.slope)),F(str(args.death_weight))) if args.canonical else False
        bits=(np.arange(513)>=256).astype(int)
        covered_template=np.logical_or(bits[:,None],bits[None,:])
        template_deterministic=512*(-args.slope*bits[None,:]+args.death_weight*bits[:,None]*(1-bits[None,:]))
        maxima=np.maximum(np.arange(513)[:,None],np.arange(513)[None,:])/512
        minimum_deaths=np.maximum(np.arange(513)[:,None]-np.arange(513)[None,:],0)
        def make_kernel(p):
            if args.layercake:
                anchors=[]
                for k in (0.,.5,1.):
                    gamma=b**k;z=1-p+p*gamma;tilted=p*gamma/z
                    value=full.power(local(tilted,b,1,X,b**(-args.death_weight),canonical),128)
                    anchors.append(value+512*math.log2(z))
                # log moment is convex in k. Chords lie ABOVE its value.
                fraction=2*np.minimum(maxima,.5)
                lower=(1-fraction)*anchors[0]+fraction*anchors[1]
                fraction=2*np.maximum(maxima-.5,0)
                upper=(1-fraction)*anchors[1]+fraction*anchors[2]
                kernel=np.where(maxima<=.5,lower,upper)
                return kernel+args.death_weight*minimum_deaths*args.logb
            kernel=full.power(local(p,b,reward,X,b**(-args.death_weight),canonical),128)
            if args.reference_burst>1:
                gamma=args.reference_burst;z=1-p+p*gamma;tilted=p*gamma/z
                alternative=full.power(local(tilted,b,reward,X,b**(-args.death_weight),canonical),128)
                alternative+=512*math.log2(z)
                kernel=np.where(covered_template,alternative,kernel)
                kernel+=template_deterministic*math.log2(gamma)
            return kernel
        for i in range(args.p_levels):
            exponent=i/(args.p_levels-1)
            if args.reference_burst>1 or args.layercake:exponent=2*exponent-1
            p=1/(1+1/(odds*b**exponent))
            kernel=make_kernel(p)
            direct=kernel.copy()
            pi=quiet.log_choose(8192,r)+r*math.log2(p)+(8192-r)*math.log2(1-p)
            for step in range(4):
                if step==3:records.append((direct.copy(),pi,p))
                direct,_=log_product(direct,direct)
            regional=np.minimum(regional,direct-pi)
            print(dict(phase='geometry cut moment',p=p,seconds=time.monotonic()-start),flush=True)
        refined=midpoint_envelope(records)
        assert np.max(refined-regional)<1e-7
        for name,matrix in (('endpoints',regional),('midpoint',refined)):
            terminal=-np.arange(513)*math.log2(X)
            for _ in range(256):terminal=np.logaddexp2.reduce(matrix+terminal[None,:],axis=1)
            union_log=quiet.log_choose(8192,r)+float(terminal[0])
            union_log+=3702*math.log2(args.reference_burst)
            if args.layercake:union_log+=3702*args.logb
            slack=(union_log+60)/args.logb
            print(dict(operator=name,slope=args.slope,logb=args.logb,inactive_rows=r,
                       death_weight=args.death_weight,
                       minimum_cardinality_cut_filter=args.canonical,
                       canonical_rule='quartet erasure and inactive births' if args.canonical else None,
                       reference_burst_marker=args.reference_burst,
                       averaged_threshold_reference=args.layercake,
                       separate_burst_allowance='2^-60' if args.reference_burst>1 or args.layercake else None,
                       log2_union_moment=union_log,slack_for_60bit_screen=slack,
                       all_support_histories=not args.canonical,
                       covers_a_minimizing_cut_for_every_geometry=True,all_row_sets_of_this_size=True,
                       directed=False,distance_certificate=False),flush=True)
        if args.flow:
            left=np.full(513,-np.inf);left[0]=0
            right=-np.arange(513)*math.log2(X)
            for _ in range(128):
                left=np.logaddexp2.reduce(left[:,None]+refined,axis=0);left-=max(left)
                before_right=right
                right=np.logaddexp2.reduce(refined+right[None,:],axis=1);right-=max(right)
            edges=left[:,None]+refined+before_right[None,:]
            u,v=map(int,np.unravel_index(np.argmax(edges),edges.shape))
            for u,v in ((u,v),(v,u)):
                candidates=np.array([half[u,:]+half[:,v]-pi for half,pi,p in records])
                scores=np.min(candidates,axis=0);w=int(np.argmax(scores))
                index=int(np.argmin(candidates[:,w]));p=records[index][2]
                kernel=make_kernel(p)
                trace=bridge_weights(kernel,u,w,8)+bridge_weights(kernel,w,v,8)[1:]
                print(dict(artificial_cut_bridge=[u,w,v],proof_p=p,
                           conditional_mean_weights=[round(value,3) for value in trace],
                           actual_geometry_distribution=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
