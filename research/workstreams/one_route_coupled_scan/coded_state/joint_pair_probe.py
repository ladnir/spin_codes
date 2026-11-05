"""FLOAT only: retain the full fresh-pair rank in a weighted output count.

Small bounded arrays (<32 MiB), serial run, no generated files. This is a
diagnostic of an upper bound, not directed arithmetic or a certificate.
"""
import math
import time
import numpy as np
from changing_pair_probe import P
from feedback_support_probe import K,W,logchoose,route_log


def pair_dimensions(phase1=False):
    u = np.arange(129,dtype=np.int16)[:,None]
    v = np.arange(129,dtype=np.int16)[None,:]
    ru,rv = np.minimum(64,u),np.minimum(64,v)
    total = (u+v).ravel()
    result = np.empty((64,64,257),dtype=np.int16)
    for h0 in range(64):
        for h1 in range(64):
            separate = np.maximum(0,h0+ru-64)+np.maximum(0,h1+rv-64)
            if phase1:
                # Each 64-row coset of A1 B_W is MDS for |W|<=64;
                # for larger W retain an arbitrary 64 columns.
                first = np.minimum((h0+1)//2,ru)
                separate = np.maximum(separate,first+np.maximum(0,h1+rv-64))
            joint = np.maximum(0,min(64,h0+h1)+ru-64)+np.maximum(0,h0+h1+rv-128)
            ranks = np.maximum(separate,joint)
            ceiling = np.full(257,-1,dtype=np.int16)
            np.maximum.at(ceiling,total,(u+v-ranks).ravel())
            result[h0,h1,:] = ceiling
    return result


def logsum(values,axis=None):
    return np.logaddexp2.reduce(values,axis=axis)


def screen(phase1=False):
    started = time.monotonic()
    dimensions = pair_dimensions(phase1)
    h = np.arange(64)
    pair_h = h[:,None]+h[None,:]
    c128 = np.array([logchoose(128,int(j)) for j in h])
    c256 = np.array([logchoose(256,int(j)) for j in h])
    # Positive output marker z=q^-beta. Bounding each local support sum
    # by 2^(number of coordinates) times its largest q-exponent is valid.
    proposals = []
    for beta in (0.48,0.49,0.50,64/127,0.51,0.52,0.54,0.56,0.60):
        same = np.maximum((64-h)*(1-beta),128*(1-beta)-h)
        joint = np.max(dimensions-beta*np.arange(257),axis=2)
        for tau in (4,8,9,10,11,12,16):
            bits = math.log2(tau)
            for xbits in (-4,-2,-1,0,1,2,3,4,6,8):
                a = float(logsum(c256+bits*same+h*xbits))
                b = float(logsum(c128[:,None]+c128[None,:]+bits*joint+pair_h*xbits,axis=None))
                c = float(logsum(c128+bits*same+h*xbits))
                proposals.append((beta,tau,xbits,a,b,c))
    for d in (0,32,64,128,256,1024,8191):
        top = 63*(8192+d)
        for H in (top-2048,top):
            def score(p):
                beta,tau,xbits,a,b,c = p
                return (logchoose(P,d)+math.log2(tau)*(beta*W-K+8193+d+128*(1-beta))
                        +(P-d)*a+d*b+c-H*xbits+route_log(H,math.log2(tau)))
            winner = min(proposals,key=score)
            print(f'd={d}, H={H}, FLOAT log2={score(winner):.3f}, '
                  f'(beta,tau,log2x)={winner[:3]}',flush=True)
        if time.monotonic()-started>45:
            raise RuntimeError('45-second diagnostic cap')
    print(f'Elapsed {time.monotonic()-started:.3f}s; selected FLOAT points only.')


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase1',action='store_true')
    args = parser.parse_args()
    screen(args.phase1)
