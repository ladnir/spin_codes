"""Binary64 one-group screen for all ranks using worst local column shapes.

Preserves exact union-support size. This is NOT an outward certificate.
"""
from math import comb,exp,log,log2
import numpy as np

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS
from rank_two_moment import census,region_transfers


def worst_regions(transfers):
    shapes,empty,active=transfers
    return [empty]+[np.max(active[[i for i,shape in enumerate(shapes)
                                  if sum(w!=0 for w in shape)==b]],axis=0) for b in range(1,5)]


def probabilities(all_regions):
    best=np.ones((64,4,257))
    for tilt,region in all_regions:
        current=np.ones((1,7))
        first=np.stack([r[0] for r in region[1:]],axis=1)
        transposed=[r.T*comb(4,b) for b,r in enumerate(region)]
        for remaining in range(64):
            values=current@first
            for b in range(1,5):
                for u in range(b,min(257,4*remaining+b+1)):
                    # A screen only. A certified replay must use outward
                    # arithmetic instead of binary64 and this underflow floor.
                    moment=max(float(values[u-b,b-1]),1e-300)
                    logbound=log(moment)+tilt*209715-log(comb(4*remaining,u-b))
                    best[remaining,b-1,u]=min(best[remaining,b-1,u],exp(min(0,logbound)))
            if remaining==63:
                break
            shifted=[current@r for r in transposed]
            count=current.shape[0]
            current=np.zeros((count+4,7))
            for b,row in enumerate(shifted):
                current[b:b+count]+=row
        print('Worst-shape diagnostic tilt',tilt,flush=True)
    result=[1.0]+[sum(comb(4,b)*comb(4*l,u-b)/comb(256,u)*best[l,b-1,u]
                     for l in range(64) for b in range(1,5) if 0<=u-b<=4*l) for u in range(1,257)]
    # CDF summation needs decreasing weights: take a decreasing majorant,
    # not differences of upper CDFs interpreted as true shell multiplicities.
    for u in range(255,-1,-1):
        result[u]=max(result[u],result[u+1])
    return result


def main():
    spectrum=authenticated_caps()
    caps=support_caps(spectrum,g=4,dimensions=dimension_caps())
    data=census()
    weights=probabilities([(float(t),worst_regions(region_transfers(data,float(t)))) for t in TILTS])
    for h,row in enumerate(caps,1):
        terms=[row[u]*(weights[u]-weights[u+1]) for u in range(256)]+[row[-1]*weights[-1]]
        total=2048*sum(terms)
        dominant=sorted(enumerate(terms),key=lambda p:p[1],reverse=True)[:4]
        print('rank',h,'log2 one-group union',log2(total),
              'dominant support CDF terms',[(u,log2(2048*v)) for u,v in dominant if v],flush=True)
    print('Binary64 only; no full certificate. Multiple active groups not included.')


if __name__=='__main__':
    main()
