"""FLOAT full-history screen retaining three quarter-region weight classes.

Positive upper envelopes, NOT a certificate. Fixed construction. Serial,
120-second watchdog, 128 MiB numeric-array cap, >=20 GiB free; no files.
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
from quartet_regional_matrix import log_product
from quartet_remaining_activity import midpoint_envelope


def partition(size,width):
    return [np.array([0])]+[np.arange(lo,min(size,lo+width)) for lo in range(1,size,width)]


def coarsen(matrix,bins,potential):
    """Max over incoming states, sum outgoing states, after conjugation."""
    transformed=matrix+potential[None,:]-potential[:,None]
    result=np.empty((len(bins),len(bins)))
    for j,dest in enumerate(bins):
        sums=np.logaddexp2.reduce(transformed[:,dest],axis=1)
        for i,source in enumerate(bins):result[i,j]=max(sums[source])
    return result


def quarter_classes(records,chunk=4):
    """Minimize the complete-region p bound before summing three classes.

    records contain four-update log matrices and complete quota log masses.
    A single p applies to all four factors for each fixed class path.
    """
    size=len(records[0][0]);result=np.full((size,size),-np.inf)
    assert chunk*size**3*8*3<48*(1<<20)
    for a in range(size):
        for lo in range(0,size,chunk):
            hi=min(size,lo+chunk)
            best=np.full((hi-lo,size,size,size),np.inf)
            for matrix,pi in records:
                term=(matrix[a,lo:hi,None,None,None]
                      +matrix[lo:hi,:,None,None]
                      +matrix[None,:,:,None]+matrix[None,None,:,:]-pi)
                np.minimum(best,term,out=best)
            block=np.logaddexp2.reduce(best,axis=0)
            block=np.logaddexp2.reduce(block,axis=0)
            block=np.logaddexp2.reduce(block,axis=0)
            result[a]=np.logaddexp2(result[a],block)
    return result


def tests():
    size=3;bins=partition(size,2);h=[F(1),F(2),F(3)]
    coeff=[np.array([[F(1+i+2*j+k,23) for j in range(size)]
                     for i in range(size)],dtype=object) for k in (0,1)]
    exact=[np.eye(size,dtype=object)]
    for _ in range(4):
        nxt=[np.zeros((size,size),dtype=object) for _ in range(len(exact)+1)]
        for i,left in enumerate(exact):
            for j,right in enumerate(coeff):nxt[i+j]+=left@right
        exact=nxt
    checks=0
    for r in (1,2,3):
        records=[];rational=[]
        for p in (F(1,5),F(2,5),F(4,5)):
            matrix=(1-p)*coeff[0]+p*coeff[1]
            small=[[max(sum(matrix[u,v]*h[v]/h[u] for v in dest)
                        for u in source) for dest in bins] for source in bins]
            pi=math.comb(4,r)*p**r*(1-p)**(4-r)
            got=coarsen(np.log2(np.array(matrix,dtype=float)),bins,np.log2(np.array(h,dtype=float)))
            assert np.max(abs(got-np.log2(np.array(small,dtype=float))))<1e-12
            records.append((got,math.log2(float(pi))));rational.append((small,pi))
        result=quarter_classes(records,1)
        assert np.max(result-quarter_classes(records[:2],1))<1e-12
        for a,source in enumerate(bins):
            for b,dest in enumerate(bins):
                bound=F(0)
                for path in product(range(len(bins)),repeat=3):
                    states=(a,)+path+(b,)
                    bound+=min(math.prod(matrix[u][v] for u,v in zip(states,states[1:]))/pi
                               for matrix,pi in rational)
                assert abs(result[a,b]-math.log2(float(bound)))<1e-11
                truth=max(sum(exact[r][u,v]*h[v]/h[u] for v in dest)
                          for u in source)/math.comb(4,r)
                assert truth<=bound;checks+=1
    print(dict(exact_quarter_class_cells=checks,additional_p_monotonicity_cells=checks,
               construction_unchanged=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(120,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--width',type=int,default=16)
        parser.add_argument('--theta',type=float,default=.13)
        parser.add_argument('--ell',type=float,default=.03)
        parser.add_argument('--tau',type=float,default=0.)
        parser.add_argument('--rows',type=int,default=2100)
        parser.add_argument('--potential',choices=('continuation','neutral'),default='continuation')
        parser.add_argument('--refine-bracket',type=int,default=-1)
        args=parser.parse_args();tests()
        if args.tests_only:return
        assert args.width in (8,16,32,64) and 0<args.theta<=1 and args.ell>0
        assert 0<args.rows<8192
        assert 0<=args.tau<=1 and (1<<32)*args.ell>=1
        assert -1<=args.refine_bracket<8
        r=args.rows;n=8192-r;theta=args.theta;ell=args.ell;X=97905/432928
        logscore=theta*args.tau*(32+math.log2(ell));score=2**logscore
        rho=r/8192;odds=rho/(1-rho);halves=[];quarters=[]
        positions=list(range(9))
        if args.refine_bracket>=0:positions+= [args.refine_bracket+j/9 for j in range(1,9)]
        positions.sort()
        for completed,i in enumerate(positions,1):
            p=1/(1+1/(odds*2**(32*theta*i/8)))
            kernel=full.power(full.local(p,theta,ell,X,score),128)
            square,_=log_product(kernel,kernel)
            quarter,_=log_product(square,square)
            half,_=log_product(quarter,quarter)
            pi=quiet.log_choose(8192,r)+r*math.log2(p)+n*math.log2(1-p)
            halves.append((half,pi,p));quarters.append(quarter)
            print(dict(phase='quarter matrices',completed=completed,proof_p=p,seconds=time.monotonic()-start),flush=True)
        envelope=midpoint_envelope(halves)
        terminal=-np.arange(513)*(math.log2(X)+theta*math.log2(ell)-logscore)
        potential=terminal.copy()
        if args.potential=='continuation':
            for _ in range(64):
                potential=np.logaddexp2.reduce(envelope+potential[None,:],axis=1)
                potential-=max(potential)
        potential-=potential[0]
        bins=partition(513,args.width)
        records=[(coarsen(matrix,bins,potential),record[1]) for matrix,record in zip(quarters,halves)]
        # Full matrices plus chunk tensors remain below the stated cap.
        assert (len(quarters)+len(halves)+10)*513*513*8+4*len(bins)**3*8*3<128*(1<<20)
        matrix=quarter_classes(records)
        tail=np.array([max((terminal-potential)[part]) for part in bins])
        for _ in range(256):tail=np.logaddexp2.reduce(matrix+tail[None,:],axis=1)
        base=quiet.log_choose(8192,r)+theta+32*theta*(943718-1048576+128*r)+theta*943718*math.log2(ell)
        base-=943718*logscore
        original=terminal.copy()
        for _ in range(256):original=np.logaddexp2.reduce(envelope+original[None,:],axis=1)
        print(dict(inactive_rows=r,bin_width=args.width,bin_count=len(bins),potential=args.potential,
                   refined_bracket=args.refine_bracket,proof_p_count=len(positions),
                   theta=theta,ell=ell,tau=args.tau,midpoint_log2_bound=float(base+original[0]),
                   quarter_class_log2_bound=float(base+tail[0]+potential[0]),
                   directed=False,whole_code=False,seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
