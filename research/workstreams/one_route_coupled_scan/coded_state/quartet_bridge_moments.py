"""FLOAT moments of the artificial growth/decay bridges; no code claim.

Serial,90s watchdog,64MiB numeric-array budget,20GiB free, no files.
The moments belong to a positive proof kernel, NOT actual encoder states.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np
import quartet_full_weight_kernel as full
import quartet_quiet_regions as quiet


def bridge(kernel,moments,start,end,length):
    size=len(kernel);left=[np.full(size,-np.inf)];left[0][start]=0
    for _ in range(length):left.append(np.logaddexp2.reduce(left[-1][:,None]+kernel,axis=0))
    right=[None]*(length+1);right[-1]=np.full(size,-np.inf);right[-1][end]=0
    for i in range(length-1,-1,-1):right[i]=np.logaddexp2.reduce(kernel+right[i+1][None,:],axis=1)
    z=left[-1][end];totals={name:0. for name in moments};rows=[]
    for i in range(length):
        mass=np.exp2(left[i][:,None]+kernel+right[i+1][None,:]-z)
        assert abs(float(mass.sum())-1)<1e-8
        row={name:float(np.sum(mass*value)) for name,value in moments.items()}
        row['next_activity']=float(mass.sum(axis=0)@np.arange(size))
        for name in totals:totals[name]+=row[name]
        rows.append(row)
    totals['next_activity']=sum(row['next_activity'] for row in rows)
    return float(z),totals,rows


def tests():
    weights=((F(1,3),F(2,7)),(F(1,5),F(3,8)))
    marked=np.array([[0,2],[3,1]],dtype=float)
    logz,got,_=bridge(np.log2(np.array(weights,dtype=float)),{'mark':marked},0,1,3)
    total=F(0);weighted=F(0)
    for mid in product((0,1),repeat=2):
        path=(0,)+mid+(1,);mass=F(1);stat=0
        for u,v in zip(path,path[1:]):mass*=weights[u][v];stat+=int(marked[u,v])
        total+=mass;weighted+=mass*stat
    assert abs(logz-math.log2(float(total)))<1e-12
    assert abs(got['mark']-float(weighted/total))<1e-12
    p=.4;theta=.175;ell=.00377;X=.23
    terms=full.local(p,theta,ell,X,1)
    marks={name:full.local(p,theta,ell,X,1,statistic=name)
           for name in ('free','active_coverage','covered_inactive','inactive_labels')}
    kernel,moments=full.power(terms,2,marks)
    assert np.max(abs(moments['free']+moments['active_coverage']+moments['covered_inactive']-8))<1e-10
    eps=1e-5;odds=p/(1-p)
    pp=odds*math.exp(eps)/(1+odds*math.exp(eps));pm=odds*math.exp(-eps)/(1+odds*math.exp(-eps))
    plus=full.power(full.local(pp,theta,ell,X,1),2)
    minus=full.power(full.local(pm,theta,ell,X,1),2)
    derivative=(plus-minus)*math.log(2)/(2*eps)+8*p
    assert np.max(abs(derivative-moments['inactive_labels']))<1e-7
    print(dict(exact_bridge_checks=2,coverage_partition_cells=81,label_derivative_cells=81),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--p',type=float,default=.7059455208680198)
        parser.add_argument('--theta',type=float,default=7/40)
        parser.add_argument('--ell',type=float,default=.00377)
        parser.add_argument('--profile',type=int,nargs=3,action='append')
        parser.add_argument('--trace',action='store_true')
        args=parser.parse_args();tests()
        if args.tests_only:return
        p=args.p;theta=args.theta;ell=args.ell;X=97905/432928;r=2100
        assert 0<p<1 and 0<theta<=1 and ell>0
        profiles=args.profile or ((0,348,489),(489,332,0))
        assert 1<=len(profiles)<=2 and all(0<=w<=512 for profile in profiles for w in profile)
        names=('free','deaths','covered_inactive','inactive_labels')
        terms=full.local(p,theta,ell,X,1)
        marks={name:full.local(p,theta,ell,X,1,statistic=name) for name in names}
        kernel,moments=full.power(terms,128,marks)
        results=[]
        pi=quiet.log_choose(8192,r)+r*math.log2(p)+(8192-r)*math.log2(1-p)
        for u,w,v in profiles:
            a,left,rowsa=bridge(kernel,moments,u,w,8)
            b,right,rowsb=bridge(kernel,moments,w,v,8)
            total={name:left[name]+right[name] for name in left}
            total['active_coverage']=8192-total['free']-total['covered_inactive']
            total['coverage_plus_deaths']=total['covered_inactive']+total['deaths']
            results.append(total)
            record=dict(profile=[u,w,v],proof_p=p,log2_class_bound=a+b-pi,
                       quota_required=r,quota_residual=total['inactive_labels']-r,
                       artificial_totals=total,actual_encoder_distribution=False)
            if args.trace:record['per_update']=rowsa+rowsb
            print(record,flush=True)
        pair={name:sum(result[name] for result in results) for name in results[0]}
        print(dict(repeated_profiles_artificial_totals={name:256/len(profiles)*value for name,value in pair.items()},
                   profiles_form_cycle=all(profiles[i][2]==profiles[(i+1)%len(profiles)][0]
                                           for i in range(len(profiles))),
                   rescaling_is_not_a_full_path_distribution=True,
                   actual_required_active_coverage=129*6092,actual_total_inactive_labels=256*r,
                   remaining_coverage_plus_deaths_cap=255463,
                   actual_encoder_distribution=False,files_written=0,seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
