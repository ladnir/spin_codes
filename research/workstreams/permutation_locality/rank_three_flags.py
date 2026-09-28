"""Binary64 sparse rank-three screen using a distinguished rank-two subspace.

For the all-one column pattern, its annihilating subspace H has dimension 2.
If supp(U)=u and supp(H)=v, the pattern count is u-v. Its four coset words
average weight u-v/2. This couples dangerous columns to the outer counts.
"""
import argparse
from collections import Counter
from itertools import combinations_with_replacement,combinations,product
from math import comb,log,log2
import numpy as np

from bch_joint_support import authenticated_caps,support_caps
from basis_lattice import improve_caps
from shortened_bound import dimension_caps
from joint_support import span,gf2_rank
from rank_two_types import subspace_cap
from rank_two_moment import census,region_transfers
from random_group_moment import probabilities as unrestricted_probabilities,worst_regions

TILTS=('0.00016','0.0002','0.00025','0.00032','0.00034','0.00036','0.00038','0.0004',
       '0.0005','0.00064','0.0008','0.001','0.00125','0.0016','0.002','0.0025')
PENALTIES=(1,.875,.75,.625,.5,.375,.25,.125)


def flag_test():
    for rows,n in (([15,51,85],7),([0x97,0x4b,0x2d,0x1e],8),([1,2,4],3)):
        words=span(rows)
        triples={tuple(sorted(span((a,b,c)))) for a,b,c in combinations(words[1:],3) if c not in (a,b,a^b)}
        for space in triples:
            union=0
            for w in space:
                union|=w
            u=union.bit_count()
            hyperplanes={tuple(sorted((0,a,b,a^b))) for a,b in combinations(space[1:],2)}
            assert len(hyperplanes)==7
            for h in hyperplanes:
                v=(h[1]|h[2]|h[3]).bit_count()
                coset=set(space)-set(h)
                assert sum(w.bit_count() for w in coset)==4*u-2*v
                assert min(w.bit_count() for w in coset)<=(2*u-v)//2
        spectrum=[sum(w.bit_count()==i for w in words) for i in range(n+1)]
        hspaces={tuple(sorted((0,a,b,a^b))) for a,b in combinations(words[1:],2)}
        hcounts=Counter((h[1]|h[2]|h[3]).bit_count() for h in hspaces)
        actual=Counter()
        absent=total=0
        for tup in product(words,repeat=4):
            if gf2_rank(tup)!=3:
                continue
            total+=1
            normal=None
            for mask in range(1,16):
                value=0
                for j in range(4):
                    if mask>>j&1:
                        value^=tup[j]
                if not value:
                    normal=mask
                    break
            assert normal is not None
            if normal.bit_count()%2:
                absent+=1
                continue
            union=(tup[0]|tup[1]|tup[2]|tup[3]).bit_count()
            all_one=(tup[0]&tup[1]&tup[2]&tup[3]).bit_count()
            actual[union,union-all_one]+=1
        assert absent*15==8*total
        for v,hcount in hcounts.items():
            for u in range(v,n+1):
                count=sum(c for (uu,vv),c in actual.items() if uu<=u and vv==v)
                assert count<=168*hcount*sum(spectrum[1:(2*u-v)//2+1])
    print('Rank-three hyperplane/coset mean identities checked on small codes',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,default=100)
    args=parser.parse_args()
    flag_test()
    maximum=args.maximum
    assert 67<=maximum<=200
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)[2]
    shells=[0]*(maximum+1)
    allowed=[w for w in range(1,257) if spectrum[w]]
    for weights in combinations_with_replacement(allowed,3):
        if sum(weights)%2==0 and sum(weights)//2<=maximum:
            shells[sum(weights)//2]+=subspace_cap(weights,spectrum)
    values=np.arange(57,maximum+1)
    # One additional class for hyperplanes in F2^4 excluding the all-one vector.
    best=np.ones((len(values)+1,64,4,257))
    data=census()
    full_transfers=[]
    for tilt in map(float,TILTS):
        transfers=region_transfers(data,tilt)
        full_transfers.append((tilt,worst_regions(transfers)))
        shapes,empty,active=transfers
        for penalty in PENALTIES+(0,):
            penalized=active*np.array([penalty**s.count(4) for s in shapes])[:,None,None]
            region=[empty]+[np.max(penalized[[i for i,s in enumerate(shapes) if sum(w!=0 for w in s)==b]],axis=0)
                            for b in range(1,5)]
            current=np.ones((1,7))
            first=np.stack([r[0] for r in region[1:]],axis=1)
            transposed=[r.T*comb(4,b) for b,r in enumerate(region)]
            for l in range(64):
                moments=current@first
                for b in range(1,5):
                    last=min(maximum if penalty else 256,4*l+b)
                    if last<67:
                        continue
                    u=np.arange(67,last+1)
                    base=np.log(np.maximum(moments[u-b,b-1],1e-300))+tilt*209715-np.array([log(comb(4*l,int(w-b))) for w in u])
                    if penalty:
                        bounds=np.exp(np.minimum(0,base[None,:]-(u[None,:]-values[:,None])*log(penalty)))
                        best[:-1,l,b-1,67:last+1]=np.minimum(best[:-1,l,b-1,67:last+1],bounds)
                    else:
                        best[-1,l,b-1,67:last+1]=np.minimum(best[-1,l,b-1,67:last+1],np.exp(np.minimum(0,base)))
                if l==63:
                    break
                shifted=[current@r for r in transposed]
                count=current.shape[0]
                current=np.zeros((min(maximum+1 if penalty else 257,count+4),7))
                for b,row in enumerate(shifted):
                    size=min(count,len(current)-b)
                    if size>0:
                        current[b:b+size]+=row[:size]
        print('Flag diagnostic tilt',tilt,flush=True)
    probabilities=np.ones((len(values)+1,maximum+1))
    for u in range(67,maximum+1):
        probabilities[:,u]=sum(comb(4,b)*comb(4*l,u-b)/comb(256,u)*best[:,l,b-1,u]
                              for l in range(64) for b in range(1,5) if 0<=u-b<=4*l)
    for u in range(maximum-1,-1,-1):
        probabilities[:,u]=np.maximum(probabilities[:,u],probabilities[:,u+1])
    outer_cdf=np.cumsum(np.array(spectrum,dtype=object))
    # Every 3-space has 2520 spanning row tuples. For a fixed hyperplane H,
    # exactly 168 row maps send its annihilator to the all-one pattern.
    # Eight of the 15 image hyperplanes exclude that pattern altogether.
    results=[]
    for index,v in enumerate(values):
        counts=[0 if u<max(67,v) else min(caps[u],168*shells[v]*int(outer_cdf[(2*u-v)//2]))
                for u in range(maximum+1)]
        p=probabilities[index]
        union=2048*(counts[-1]*p[-1]+sum(counts[u]*(p[u]-p[u+1]) for u in range(maximum)))
        results.append((int(v),union))
    absent=[(caps[u]*8)//15 for u in range(257)]
    # caps counts all embeddings; the exact unknown CDF is divisible by 15.
    p=[1.0]*67+[sum(comb(4,b)*comb(4*l,u-b)/comb(256,u)*best[-1,l,b-1,u]
                    for l in range(64) for b in range(1,5) if 0<=u-b<=4*l) for u in range(67,257)]
    for u in range(255,-1,-1):
        p[u]=max(p[u],p[u+1])
    no_four=2048*(absent[-1]*p[-1]+sum(absent[u]*(p[u]-p[u+1]) for u in range(256)))
    sparse=sum(v for _,v in results)+no_four
    unrestricted=unrestricted_probabilities(full_transfers)
    present=[(cap*7)//15 for cap in caps]
    dense=2048*(present[-1]*unrestricted[-1]+sum(present[u]*(unrestricted[u]-unrestricted[u+1]) for u in range(maximum+1,256)))
    print('Sparse flags plus all no-all-one tuples log2',log2(sparse),'no-all-one component',log2(no_four),flush=True)
    print('Dominant hyperplane supports',[(v,log2(x)) for v,x in sorted(results,key=lambda p:p[1],reverse=True)[:8] if x],flush=True)
    print('Dense remainder log2',log2(dense),'combined log2',log2(sparse+dense),flush=True)
    print('Binary64 diagnostic only; no new certified margin.')


if __name__=='__main__':
    main()
