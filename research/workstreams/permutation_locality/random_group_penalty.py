"""Binary64 screen retaining a rank-dependent cap on all-one columns.

For rank h, at most u-d_(h-1) union columns equal 1111. Use an exponential
weight on those columns before taking a worst-shape envelope. No cert here.
"""
import argparse
from itertools import product
from math import comb,log,log2
import numpy as np

from bch_joint_support import authenticated_caps,support_caps
from basis_lattice import improve_caps
from shortened_bound import dimension_caps
from joint_support import gf2_rank,griesmer,span
from group_rank_one_verify import TILTS
from rank_two_moment import census,region_transfers


def support_test():
    for rows,n in (([15,51,85],7),([0x97,0x4b,0x2d,0x1e],8),([1,2,4],3)):
        words=span(rows)
        d=min(w.bit_count() for w in words[1:])
        for tup in product(words,repeat=4):
            h=gf2_rank(tup)
            if h<2:
                continue
            union=tup[0]|tup[1]|tup[2]|tup[3]
            all_one=tup[0]&tup[1]&tup[2]&tup[3]
            assert all_one.bit_count()<=union.bit_count()-griesmer(d,h-1)
    print('All-one column cap passed exhaustive small-code checks',flush=True)


def probabilities(data,penalties=('1','0.75','0.5','0.25'),tilts=TILTS):
    # The current verified lower generalized weights for ranks 1,2,3.
    previous=(38,57,67)
    best=np.ones((3,64,4,257))
    for tilt in map(float,tilts):
        shapes,empty,active=region_transfers(data,tilt)
        for penalty_text in penalties:
            penalty=float(penalty_text)
            penalized=active*np.array([penalty**s.count(4) for s in shapes])[:,None,None]
            region=[empty]+[np.max(penalized[[i for i,s in enumerate(shapes) if sum(w!=0 for w in s)==b]],axis=0)
                            for b in range(1,5)]
            current=np.ones((1,7))
            first=np.stack([r[0] for r in region[1:]],axis=1)
            transposed=[r.T*comb(4,b) for b,r in enumerate(region)]
            for l in range(64):
                values=current@first
                for b in range(1,5):
                    for u in range(max(b,38),min(257,4*l+b+1)):
                        moment=max(float(values[u-b,b-1]),1e-300)
                        base=log(moment)+tilt*209715-log(comb(4*l,u-b))
                        for h,d in enumerate(previous):
                            if u>=d:
                                logbound=base-(u-d)*log(penalty)
                                best[h,l,b-1,u]=min(best[h,l,b-1,u],np.exp(min(0,logbound)))
                if l==63:
                    break
                count=current.shape[0]
                shifted=[current@r for r in transposed]
                current=np.zeros((count+4,7))
                for b,row in enumerate(shifted):
                    current[b:b+count]+=row
        print('Column-penalty diagnostic tilt',tilt,flush=True)
    result=[]
    for h in range(3):
        weights=[1.0]+[sum(comb(4,b)*comb(4*l,u-b)/comb(256,u)*best[h,l,b-1,u]
                          for l in range(64) for b in range(1,5) if 0<=u-b<=4*l) for u in range(1,257)]
        for u in range(255,-1,-1):
            weights[u]=max(weights[u],weights[u+1])
        result.append(weights)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fine',action='store_true')
    args=parser.parse_args()
    support_test()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    parameters={}
    if args.fine:
        parameters=dict(penalties=('1','0.875','0.75','0.625','0.5','0.375','0.25','0.125'),
                        tilts=('0.00016','0.0002','0.00025','0.00032','0.0004','0.0005',
                               '0.00064','0.0008','0.001','0.00125','0.0016','0.002','0.0025'))
    weights=probabilities(census(),**parameters)
    for h,(row,prob) in enumerate(zip(caps[1:],weights),2):
        terms=[row[u]*(prob[u]-prob[u+1]) for u in range(256)]+[row[-1]*prob[-1]]
        total=2048*sum(terms)
        dominant=sorted(enumerate(terms),key=lambda p:p[1],reverse=True)[:4]
        print('rank',h,'log2 one-group union',log2(total),
              'dominant CDF terms',[(u,log2(2048*v)) for u,v in dominant if v],flush=True)
    print('Binary64 diagnostic only. Multiple active groups not included.')


if __name__=='__main__':
    main()
