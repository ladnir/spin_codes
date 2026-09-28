"""Binary64 rank-four screen conditioned on the all-one annihilator support.

The outer counts use integers. Moment arithmetic is diagnostic, not outward.
The three-space shell count is bounded by a CDF, never a difference of caps.
"""
from collections import Counter
from itertools import combinations, product
from math import comb, log, log2
import argparse
import numpy as np

from bch_joint_support import authenticated_caps, support_caps
from basis_lattice import improve_caps
from shortened_bound import dimension_caps
from joint_support import span, gf2_rank
from rank_two_moment import census, region_transfers
from rank_three_flags import TILTS, PENALTIES


def self_test():
    for rows in ([1,2,4,8], [0x97,0x4b,0x2d,0x1e], [0x5555,0x3333,0x0f0f,0x00ff]):
        words=span(rows)
        whole=0
        for w in words:
            whole|=w
        u=whole.bit_count()
        hyperplanes={tuple(sorted(span(b))) for b in combinations(words[1:],3) if gf2_rank(b)==3}
        assert len(hyperplanes)==15
        hcounts=Counter()
        for h in hyperplanes:
            support=0
            for w in h:
                support|=w
            v=support.bit_count()
            hcounts[v]+=1
            coset=set(words)-set(h)
            assert len(coset)==8 and sum(w.bit_count() for w in coset)==8*u-4*v
            assert min(w.bit_count() for w in coset)<=(2*u-v)//2
        actual=Counter()
        for tup in product(words[1:],repeat=4):
            if gf2_rank(tup)!=4:
                continue
            all_one=(tup[0]&tup[1]&tup[2]&tup[3]).bit_count()
            v=((tup[0]^tup[1])|(tup[0]^tup[2])|(tup[0]^tup[3])).bit_count()
            assert all_one==u-v
            actual[v]+=1
        assert sum(actual.values())==20160
        for v,hcount in hcounts.items():
            assert actual[v]==1344*hcount
            eligible=sum(0<w.bit_count()<=(2*u-v)//2 for w in words)
            assert actual[v]<=1344*hcount*eligible
    print('Rank-four coset identity and 1344 assignments checked exhaustively',flush=True)
    for rows in ([1,2,4,8,16], [0x97,0x4b,0x2d,0x1e], [0x5555,0x3333,0x0f0f,0x00ff]):
        words=span(rows)
        n=max(words).bit_length()
        hyperplanes={tuple(sorted(span(b))) for b in combinations(words[1:],3) if gf2_rank(b)==3}
        for h in hyperplanes:
            support=0
            for w in h:
                support|=w
            v=support.bit_count()
            fiber=sum(not(w&~support) for w in words)
            cosets={tuple(sorted(w^x for x in h)) for w in words if w not in h}
            counts=Counter(((support|coset[0]).bit_count()-v) for coset in cosets)
            for m,actual in counts.items():
                upper=(comb(n-v,m)*fiber-(8 if m==0 else 0))//8
                assert actual<=upper
    print('Exact exterior-weight/fiber count checked on all three-subspaces of small codes',flush=True)


def flag_probabilities(data, maximum=256, tilts=TILTS, penalties=PENALTIES):
    values=np.arange(67,maximum+1)
    best=np.ones((len(values),64,4,maximum+1))
    for tilt in map(float,tilts):
        shapes,empty,active=region_transfers(data,tilt)
        for penalty in penalties:
            penalized=active*np.array([penalty**s.count(4) for s in shapes])[:,None,None]
            region=[empty]+[np.max(penalized[[i for i,s in enumerate(shapes) if sum(w!=0 for w in s)==b]],axis=0)
                            for b in range(1,5)]
            current=np.ones((1,7))
            first=np.stack([r[0] for r in region[1:]],axis=1)
            transposed=[r.T*comb(4,b) for b,r in enumerate(region)]
            for l in range(64):
                moments=current@first
                for b in range(1,5):
                    last=min(maximum,4*l+b)
                    if last<72:
                        continue
                    u=np.arange(72,last+1)
                    base=np.log(np.maximum(moments[u-b,b-1],1e-300))+tilt*209715-np.array([log(comb(4*l,int(w-b))) for w in u])
                    bounds=np.exp(np.minimum(0,base[None,:]-(u[None,:]-values[:,None])*log(penalty)))
                    best[:,l,b-1,72:last+1]=np.minimum(best[:,l,b-1,72:last+1],bounds)
                if l==63:
                    break
                shifted=[current@r for r in transposed]
                count=current.shape[0]
                current=np.zeros((min(maximum+1,count+4),7))
                for b,row in enumerate(shifted):
                    size=min(count,len(current)-b)
                    if size>0:
                        current[b:b+size]+=row[:size]
        print('Rank-four flag tilt',tilt,flush=True)
    probabilities=np.ones((len(values),maximum+1))
    for u in range(72,maximum+1):
        probabilities[:,u]=sum(comb(4,b)*comb(4*l,u-b)/comb(256,u)*best[:,l,b-1,u]
                              for l in range(64) for b in range(1,5) if 0<=u-b<=4*l)
    return values,probabilities


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse',action='store_true')
    args=parser.parse_args()
    self_test()
    spectrum=authenticated_caps()
    dimensions=dimension_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    # Every three-space has 2520 spanning four-row assignments. Its exact
    # shell count is at most its CDF cap divided by this exact multiplicity.
    h_shell=[cap//2520 for cap in caps[2]]
    options=dict(tilts=TILTS[::2],penalties=(1,.75,.5,.25)) if args.coarse else {}
    values,probabilities=flag_probabilities(census(),**options)
    outer_cdf=np.cumsum(np.array(spectrum,dtype=object))
    results=[]
    dominant=[]
    shell_results=[]
    shell_dominant=[]
    for index,v in enumerate(values):
        counts=[0 if u<max(72,v) else min(caps[3][u],1344*h_shell[v]*int(outer_cdf[(2*u-v)//2]))
                for u in range(257)]
        raw=probabilities[index]
        shell_terms=[]
        for u in range(max(72,v),257):
            # For fixed H with support S, restriction outside S has fibers
            # of size 2^dim(C_S), bounded by the authenticated shortening cap.
            # Each nonzero coset x+H has eight words with the same exterior.
            extensions=(comb(256-v,u-v)*(1<<dimensions[v])-(8 if u==v else 0))//8
            count=min(counts[u],1344*h_shell[v]*extensions)
            term=2048*count*raw[u]
            shell_terms.append(term)
            if term:
                shell_dominant.append((term,u,int(v),log2(count),log2(raw[u])))
        shell_results.append((int(v),sum(shell_terms)))
        p=raw.copy()
        for u in range(255,-1,-1):
            p[u]=max(p[u],p[u+1])
        terms=[2048*counts[u]*(p[u]-p[u+1]) for u in range(256)]+[2048*counts[-1]*p[-1]]
        results.append((int(v),sum(terms)))
        dominant.extend((x,u,int(v)) for u,x in enumerate(terms) if x)
    total=sum(x for _,x in results)
    print('Rank-four flags log2 one-group union',log2(total),flush=True)
    print('Dominant annihilator supports',[(v,log2(x)) for v,x in sorted(results,key=lambda p:p[1],reverse=True)[:10] if x],flush=True)
    print('Dominant (union support, annihilator support, log2 term)',[(u,v,log2(x)) for x,u,v in sorted(dominant,reverse=True)[:10]],flush=True)
    print('Exterior-shell rank-four log2 one-group union',log2(sum(x for _,x in shell_results)),flush=True)
    print('Exterior-shell dominant (u,v,log2 term,log2 count,log2 probability)',[(u,v,log2(x),c,p) for x,u,v,c,p in sorted(shell_dominant,reverse=True)[:12]],flush=True)
    print('Binary64 diagnostic only; rank four not certified and multiple groups not covered.')


if __name__=='__main__':
    main()
