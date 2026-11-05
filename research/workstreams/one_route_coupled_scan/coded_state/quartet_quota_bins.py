"""FLOAT routing diagnostic: retain count intervals within each region.

All weights and turnovers are covered, with explicit maxima over optional
weight bins. No representative states. Count intervals overapproximate the
exact inactive quota; no independence of actual regional labels is assumed.
Serial, 120-second watchdog, 128 MiB numeric-array budget, >=20 GiB free,
no generated files and no claim of a directed whole-code certificate.
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


def interval_ratio(p,rho,degree,lo,hi):
    slope=math.log2(rho/p)-math.log2((1-rho)/(1-p))
    endpoint=hi if slope>=0 else lo
    return degree*math.log2((1-rho)/(1-p))+endpoint*slope


def quota_apply(edges,vector,r,rho,count_width,degree=512,length=16):
    """Backward positive transfer, retaining sum of count-bin indices."""
    cap=r//count_width;current=vector[:,None]
    for step in range(length):
        size=min(cap+1,current.shape[1]+len(edges)-1)
        nxt=np.full((len(vector),size),-np.inf)
        for j,edge in enumerate(edges[:size]):
            count=min(current.shape[1],size-j)
            product,_=log_product(edge,current[:,:count])
            nxt[:,j:j+count]=np.logaddexp2(nxt[:,j:j+count],product)
        assert np.all(np.isfinite(nxt))
        current=nxt
    low=max(0,math.ceil((r-length*(count_width-1))/count_width))
    pi=quiet.log_choose(degree*length,r)+r*math.log2(rho)+(degree*length-r)*math.log2(1-rho)
    return np.logaddexp2.reduce(current[:,low:cap+1],axis=1)-pi


def tests():
    # Exact nonnegative toy polynomial matrices, independent of quartet code.
    degree=2;rho=F(2,5);length=3;n=3
    coefficients=[np.array([[F(1+i+2*j+t,17) for j in range(n)] for i in range(n)],dtype=object)
                  for t in range(degree+1)]
    exact=[np.eye(n,dtype=object)]
    for _ in range(length):
        nxt=[np.zeros((n,n),dtype=object) for _ in range(len(exact)+degree)]
        for a,left in enumerate(exact):
            for b,right in enumerate(coefficients):nxt[a+b]+=left@right
        exact=nxt
    checked=0
    for width in (1,2):
        edges=np.full((degree//width+1,n,n),np.inf)
        for p in (F(1,5),rho,F(4,5)):
            evaluated=sum((c*p**t*(1-p)**(degree-t) for t,c in enumerate(coefficients)))
            logs=np.log2(np.asarray(evaluated,dtype=float))
            for j in range(len(edges)):
                lo=j*width;hi=min(degree,lo+width-1)
                candidate=logs+interval_ratio(float(p),float(rho),degree,lo,hi)
                truth=sum(coefficients[t]*rho**t*(1-rho)**(degree-t) for t in range(lo,hi+1))
                assert np.all(candidate+1e-12>=np.log2(np.asarray(truth,dtype=float)))
                edges[j]=np.minimum(edges[j],candidate)
        for r in range(1,degree*length):
            got=quota_apply(edges,np.zeros(n),r,float(rho),width,degree,length)
            truth=np.array([float(sum(row)/math.comb(degree*length,r)) for row in exact[r]])
            assert np.all(got+1e-11>=np.log2(truth));checked+=n
    # A 16-update path carrying at most 31 labels per update cannot meet
    # r=2100. Its sole count-bin trajectory must be excluded, not repriced.
    width=32;r=2100
    assert math.ceil((r-16*(width-1))/width)>0
    print(dict(exact_toy_quota_components=checked,sparse_count_path_excluded=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(120,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--state-width',type=int,default=8)
        parser.add_argument('--count-width',type=int,default=16)
        parser.add_argument('--p-levels',type=int,default=9)
        parser.add_argument('--iterations',type=int,default=3)
        args=parser.parse_args();tests()
        if args.tests_only:return
        assert args.state_width in (1,2,4,8,16,32) and args.count_width in (8,16,32,64)
        assert 2<=args.p_levels<=9 and 1<=args.iterations<=6
        r=2100;rho=r/8192;theta=7/40;ell=.00377;X=97905/432928
        bins=[np.array([0])]+[np.arange(lo,min(513,lo+args.state_width)) for lo in range(1,513,args.state_width)]
        n=len(bins);labels=512//args.count_width+1
        assert labels*n*n*8+32*(1<<20)<128*(1<<20)
        edges=np.full((labels,n,n),np.inf)
        odds=rho/(1-rho)
        for index in range(args.p_levels):
            p=1/(1+1/(odds*2**(32*theta*index/(args.p_levels-1))))
            kernel=full.power(full.local(p,theta,ell,X,1),128)
            if args.state_width==1:
                coarse=kernel
            else:
                coarse=np.empty((n,n))
                for i,rows in enumerate(bins):
                    for j,cols in enumerate(bins):
                        coarse[i,j]=np.max(np.logaddexp2.reduce(kernel[np.ix_(rows,cols)],axis=1))
            for j in range(labels):
                lo=j*args.count_width;hi=min(512,lo+args.count_width-1)
                edges[j]=np.minimum(edges[j],coarse+interval_ratio(p,rho,512,lo,hi))
            print(dict(phase='count interval envelopes',p=p,seconds=time.monotonic()-start),flush=True)
        terminal=np.array([max(-b*(math.log2(X)+theta*math.log2(ell))) for b in bins])
        potential=terminal.copy();best=math.inf
        base=quiet.log_choose(8192,r)+theta+32*theta*(943718-1048576+128*r)+theta*943718*math.log2(ell)
        for iteration in range(args.iterations):
            applied=quota_apply(edges,potential,r,rho,args.count_width)
            growth=float(max(applied-potential))
            bound=base+256*growth+potential[0]+float(max(terminal-potential));best=min(best,bound)
            print(dict(iteration=iteration,log2_bound=float(bound),best=float(best),
                       state_bins=n,count_bins=labels,seconds=time.monotonic()-start,
                       all_weights=True,all_deaths=True,whole_code=False,directed=False),flush=True)
            potential=applied-applied[0]
        print(dict(files_written=0,seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
