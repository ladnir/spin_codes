"""Binary64 screen: maximize complete shape actions, not matrix entries.

Exact support coefficients, seven-state envelopes, one active group only.
No outward numerical certificate is claimed by this script.
"""
from math import comb,exp,log,log2
import numpy as np

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS
from rank_two_moment import census,region_transfers


def probabilities(data):
    best=np.ones((64,4,257))
    for tilt in map(float,TILTS):
        shapes,empty,active=region_transfers(data,tilt)
        candidates=[empty[None,:,:]]+[active[[i for i,s in enumerate(shapes) if sum(w!=0 for w in s)==b]]
                                     for b in range(1,5)]
        current=np.ones((1,7))
        for remaining in range(64):
            # Each candidate acts as a whole on the continuation vector.
            # Maximizing separately for each entering coordinate is still an
            # upper envelope, but does not splice entries from different shapes.
            transformed=[np.max(current@rows.transpose(0,2,1),axis=0) for rows in candidates]
            for b in range(1,5):
                for u in range(b,min(257,4*remaining+b+1)):
                    moment=max(float(transformed[b][u-b,0]),1e-300)
                    logbound=log(moment)+tilt*209715-log(comb(4*remaining,u-b))
                    best[remaining,b-1,u]=min(best[remaining,b-1,u],exp(min(0,logbound)))
            if remaining==63:
                break
            count=current.shape[0]
            current=np.zeros((count+4,7))
            for b,row in enumerate(transformed):
                current[b:b+count]+=comb(4,b)*row
        print('Backward shape-action diagnostic tilt',tilt,flush=True)
    result=[1.0]+[sum(comb(4,b)*comb(4*l,u-b)/comb(256,u)*best[l,b-1,u]
                     for l in range(64) for b in range(1,5) if 0<=u-b<=4*l) for u in range(1,257)]
    for u in range(255,-1,-1):
        result[u]=max(result[u],result[u+1])
    return result


def main():
    spectrum=authenticated_caps()
    caps=support_caps(spectrum,g=4,dimensions=dimension_caps())
    weights=probabilities(census())
    for h,row in enumerate(caps,1):
        terms=[row[u]*(weights[u]-weights[u+1]) for u in range(256)]+[row[-1]*weights[-1]]
        total=2048*sum(terms)
        dominant=sorted(enumerate(terms),key=lambda p:p[1],reverse=True)[:4]
        print('rank',h,'log2 one-group union',log2(total),
              'dominant CDF terms',[(u,log2(2048*v)) for u,v in dominant if v],flush=True)
    print('Binary64 diagnostic only. Multiple active groups not included.')


if __name__=='__main__':
    main()
