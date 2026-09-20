"""Adaptive weighted Bernoulli majorant, simultaneously bounding all Q.

The adaptive process is only an upper bound on fixed group assignments.
Gamma^(1/256) distributes outer counting cost across the 256 regions.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.special import gammaln
import bridge as base
import general_occupancy as general


def normalized_regions(epoch,maximum,t=128,length=8192):
    current=np.eye(3)[None,:,:]
    for completed in range(0,length,t):
        outsize=min(maximum,completed+t)+1
        updated=np.zeros((outsize,3,3))
        for a in range(min(t,maximum)+1):
            lo=max(a,0);hi=min(outsize,a+len(current))
            if hi<=lo:continue
            degree=np.arange(lo,hi);old=degree-a
            logweight=(gammaln(t+1)-gammaln(a+1)-gammaln(t-a+1)
                       +gammaln(completed+1)-gammaln(old+1)-gammaln(completed-old+1)
                       -gammaln(completed+t+1)+gammaln(degree+1)+gammaln(completed+t-degree+1))
            updated[lo:hi]+=np.exp(logweight)[:,None,None]*(current[old]@epoch[a])
        current=updated
    return current


def adaptive_logs(region,ps,roots):
    maximum=len(region)-1;current=region.copy();scale=0.;values=[]
    for q in range(1,maximum+1):
        a=current[:-1];b=current[1:]
        updated=roots[0]*((1-ps[0])*a+ps[0]*b)
        for p,r in zip(ps[1:],roots[1:]):
            np.maximum(updated,r*((1-p)*a+p*b),out=updated)
        normalizer=float(updated.max());assert normalizer>0
        current=updated/normalizer;scale+=math.log(normalizer)
        values.append(general.log_power(current[0].copy())+256*scale)
    return np.array(values)


def screen(maximum):
    t,s,spectrum=base.load_map(general.NAME);kernel=general.kernel_spectrum();caps=general.q2.deterministic_caps()
    weights=[np.array(band) for band in general.BANDS]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in general.BANDS]
    mid=[(band[0]+band[-1])/512 for band in general.BANDS]
    configurations=[]
    for shift_index in range(-2,25):
        factor=math.exp(shift_index/4)
        ps=[p*factor/(1-p+p*factor) for p in mid]
        gamma=[float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)]
        configurations.append((shift_index,ps,np.exp(np.array(gamma)/256)))
    qvalues=np.arange(1,maximum+1)
    placements=np.array([math.log(math.comb(8192,int(q)))+int(q)*math.log(5) for q in qvalues])
    best=np.full(maximum,np.inf);witness=[None]*maximum
    for tenth in range(-80,-24,5):
        lam=math.exp(tenth/10)
        epoch=np.array(general.epoch_matrices(t,s,spectrum,kernel,math.exp(-lam),float,maximum)).reshape(-1,3,3)
        region=normalized_regions(epoch,maximum)
        for shift_index,ps,roots in configurations:
            value=adaptive_logs(region,ps,roots)+209716*lam+placements
            improved=value<best
            for index in np.flatnonzero(improved):witness[int(index)]=(tenth,shift_index,ps)
            best=np.minimum(best,value)
        print('Adaptive range tilt',tenth,'margin at',maximum,-best[-1]/math.log(2),flush=True)
    rows=[dict(occupation=int(q),margin_bits_diagnostic=float(-value/math.log(2)),witness_tenth=w[0],
               shift_index=w[1],p=[base.encode(F.from_float(p)) for p in w[2]]) for q,value,w in zip(qvalues,best,witness)]
    base.write_new(base.HERE/'generated'/f'adaptive_q1_q{maximum}_screen.json',
                   dict(status='ADAPTIVE_RANGE_SCREEN_ONLY',rows=rows,source_sha256={'adaptive_range.py':base.sha(Path(__file__))}))
    print('Selected margins',[(r['occupation'],r['margin_bits_diagnostic']) for r in rows if r['occupation'] in (1,4,8,16,17,24,32,48,64,96,128)],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--maximum',type=int,default=64)
    args=parser.parse_args();assert 1<=args.maximum<=8192;screen(args.maximum)
