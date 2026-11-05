"""FLOAT all-weight screen retaining both regional endpoints before min-p.

The regional quota is bounded by its full conditioning factor. The local
count, all 513 weights, and all turnovers are unchanged. No directed or
whole-code certificate. Serial, 60s watchdog, >=20 GiB free, no files written.
Block products have a 32 MiB numeric-array budget, independent of p-grid size.
"""
import argparse
import math
import os
from pathlib import Path
import shutil
import threading
import time

import numpy as np

import quartet_full_weight_kernel as full
import quartet_quiet_regions as quiet


def log_product(a,b,chunk=16):
    """Positive log matrix product with explicit tiny-product fallback.

    Row/column scaling permits ordinary BLAS in most blocks. If the scaled
    product is tiny, recompute that entry in log space rather than accepting
    underflow as zero. This is a FLOAT screen, not directed arithmetic.
    """
    result=np.full((a.shape[0],b.shape[1]),-np.inf);fallback=0
    for lo in range(0,a.shape[1],chunk):
        left=a[:,lo:lo+chunk];right=b[lo:lo+chunk,:]
        row=np.max(left,axis=1);col=np.max(right,axis=0)
        assert np.all(np.isfinite(row)) and np.all(np.isfinite(col))
        product=np.exp2(left-row[:,None])@np.exp2(right-col[None,:])
        bad=product<1e-200
        block=row[:,None]+col[None,:]+np.log2(np.maximum(product,1e-300))
        ii,jj=np.nonzero(bad);fallback+=len(ii)
        for start in range(0,len(ii),2048):
            i=ii[start:start+2048];j=jj[start:start+2048]
            block[i,j]=np.logaddexp2.reduce(left[i,:]+right[:,j].T,axis=1)
        result=np.logaddexp2(result,block)
    assert np.all(np.isfinite(result))
    return result,fallback


def tests():
    rng=np.random.default_rng(92741);cells=0;fallbacks=0
    for scale in (1.,100.,2000.):
        a=rng.uniform(-scale,scale,(9,13));b=rng.uniform(-scale,scale,(13,7))
        expected=np.logaddexp2.reduce(a[:,:,None]+b[None,:,:],axis=1)
        for chunk in (1,4,16):
            got,falls=log_product(a,b,chunk)
            assert np.max(np.abs(got-expected))<2e-11
            cells+=got.size;fallbacks+=falls
    print(dict(log_product_cells=cells,tiny_product_fallbacks=fallbacks),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--p',nargs='+',type=float,default=[2100/8192,112627/167924,.94])
        parser.add_argument('--rows',type=int,default=2100)
        parser.add_argument('--ell',type=float,default=.00377)
        parser.add_argument('--coverage',type=float,default=1.)
        parser.add_argument('--theta',type=float,default=7/40)
        parser.add_argument('--p-levels',type=int,default=0)
        args=parser.parse_args();tests()
        if args.tests_only:return
        assert 0<=args.rows<=8192 and args.coverage>=1 and args.ell*(1<<32)>=1
        assert 0<args.theta<=1
        theta=args.theta
        if args.p_levels:
            assert 2<=args.p_levels<=12 and 0<args.rows<8192
            odds=args.rows/(8192-args.rows)
            args.p=[1/(1+1/(odds*2**(32*theta*i/(args.p_levels-1)))) for i in range(args.p_levels)]
        assert 1<=len(args.p)<=12 and all(0<p<1 for p in args.p)
        X=97905/432928;regional=np.full((513,513),np.inf)
        for p in args.p:
            kernel=full.power(full.local(p,theta,args.ell,X,1,args.coverage),128)
            full.boundary_checks(kernel,p,theta,args.ell,X,1,128,args.coverage)
            direct=kernel.copy();falls=0
            for _ in range(4):direct,n=log_product(direct,direct);falls+=n
            # Independent vector application cross-checks the complete 16
            # updates, including free endpoints and very small entries.
            vector=-np.arange(513)*(math.log2(X)+theta*math.log2(args.ell))
            reference=vector.copy()
            for _ in range(16):reference=np.logaddexp2.reduce(kernel+reference[None,:],axis=1)
            check=np.logaddexp2.reduce(direct+vector[None,:],axis=1)
            error=float(np.max(np.abs(check-reference)));assert error<1e-7
            pi=quiet.log_choose(8192,args.rows)+args.rows*math.log2(p)+(8192-args.rows)*math.log2(1-p)
            regional=np.minimum(regional,direct-pi)
            print(dict(p=p,region_vector_error=error,tiny_product_fallbacks=falls,
                       seconds=time.monotonic()-start),flush=True)
        terminal=-np.arange(513)*(math.log2(X)+theta*math.log2(args.ell))
        for _ in range(256):terminal=np.logaddexp2.reduce(regional+terminal[None,:],axis=1)
        base=(quiet.log_choose(8192,args.rows)+theta+32*theta*(943718-1048576+128*args.rows)
              +theta*943718*math.log2(args.ell)-129*(8192-args.rows)*math.log2(args.coverage))
        print(dict(log2_bound=float(base+terminal[0]),theta=theta,all_weights=True,all_deaths=True,
                   regional_endpoint_minimum=True,directed=False,whole_code=False,
                   seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
