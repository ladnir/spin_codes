"""Binary64 diagnostic for rank-two types, using exact local orbit counts.

The iid-column experiment is conditioned on its complete type composition.
This is a screening calculation, NOT an outward numerical certificate.
"""
import argparse
from collections import Counter,defaultdict
from itertools import permutations
from math import comb,exp,factorial,lgamma,log,log2,prod

import numpy as np

from group_moment import maps
from group_rank_one_verify import TILTS
from random_rank_one import shapes_and_masks
from rank_two_types import composition,embeddings,subspace_cap,self_test
from bch_joint_support import authenticated_caps


def intersection_counts(state_shape,input_shape):
    result=Counter()
    for columns in set(permutations(input_shape)):
        counts={0:1}
        for v,a in zip(state_shape,columns):
            updated=Counter()
            for old,count in counts.items():
                for overlap in range(max(0,v+a-4),min(v,a)+1):
                    updated[old+overlap]+=count*comb(v,overlap)*comb(4-v,a-overlap)
            counts=updated
        result.update(counts)
    return result


def census():
    images,columns,spectrum=maps()
    shapes,_=shapes_and_masks()
    allowed=defaultdict(list)
    for mask,shape in enumerate(shapes):
        if mask:
            allowed[shape].append(mask)
    windows=defaultdict(Counter)
    for image in images[1:]:
        for start in range(0,128,16):
            windows[image.bit_count()][shapes[(image>>start)&65535]]+=1
    moments={shape:defaultdict(Counter) for shape in allowed}
    for state_shape in set(shapes):
        for input_shape in allowed:
            counts=intersection_counts(state_shape,input_shape)
            assert sum(counts.values())==len(allowed[input_shape])
            for v in spectrum:
                frequency=windows[v][state_shape]
                for overlap,count in counts.items():
                    moments[input_shape][v][v+sum(input_shape)-2*overlap]+=frequency*count
    atoms=defaultdict(Counter)
    cancel={shape:defaultdict(Counter) for shape in allowed}
    for start in range(0,128,16):
        syndromes=[0]*65536
        for mask in range(1,65536):
            bit=mask&-mask
            syndrome=syndromes[mask^bit]^columns[start+bit.bit_length()-1]
            syndromes[mask]=syndrome
            shape=shapes[mask]
            atoms[shape][syndrome]+=1
            if syndrome:
                image=images[syndrome]
                cancel[shape][image.bit_count()][(image^(mask<<start)).bit_count()]+=1
    for shape in allowed:
        choices=8*len(allowed[shape])
        assert sum(atoms[shape].values())==choices
        assert sum(sum(row.values()) for row in cancel[shape].values())==choices-atoms[shape][0]
        for v,count in spectrum.items():
            assert sum(moments[shape][v].values())==choices*count
    assert {shape:row[0] for shape,row in atoms.items() if row[0]}=={(0,2,3,3):1,(1,2,2,3):1}
    # Independent direct overlap checks, not just conservation of total counts.
    for mask in (0x1234,0xabcd,0xaaaa,0x000f):
        for shape in allowed:
            expected=Counter((mask&x).bit_count() for x in allowed[shape])
            assert expected==intersection_counts(shapes[mask],shape)
    print('All 69 nonempty orbit moments and syndrome counts reconstructed and checked',flush=True)
    return spectrum,allowed,moments,atoms,cancel


def region_transfers(data,tilt,windows=8,epochs=256):
    spectrum,allowed,moments,atoms,cancel=data
    levels=sorted(spectrum)
    m=(1<<19)-1
    powers=np.exp(-tilt*np.arange(145))
    empty=np.zeros((7,7))
    empty[0,0]=1
    empty[1,1]=powers[48]/2
    for j,w in enumerate(levels):
        empty[1,j+2]=powers[48]*spectrum[w]/(2*m)
    for i,v in enumerate(levels):
        empty[i+2,i+2]=powers[v]/2
        for j,w in enumerate(levels):
            empty[i+2,j+2]+=powers[v]*spectrum[w]/(2*m)
    shapes=sorted(allowed)
    active=np.zeros((len(shapes),7,7))
    for k,shape in enumerate(shapes):
        choices=windows*len(allowed[shape])
        q=atoms[shape][0]/choices
        active[k,0,0]=q*powers[sum(shape)]
        active[k,0,1]=(1-q)*powers[sum(shape)]
        moment=powers[48-sum(shape)]
        maxatom=max(count for syndrome,count in atoms[shape].items() if syndrome)
        least=min(w for row in cancel[shape].values() for w in row)
        cancellation=min(moment,maxatom/choices*powers[least])
        for i in range(1,7):
            if i>=2:
                v=levels[i-2]
                moment=sum(count*powers[w] for w,count in moments[shape][v].items())/(choices*spectrum[v])
                cancellation=sum(count*powers[w] for w,count in cancel[shape][v].items())/(choices*spectrum[v])
            active[k,i,0]=cancellation/2+moment/(2*m)
            active[k,i,1]=moment/2
            for j,w in enumerate(levels):
                active[k,i,j+2]=moment*spectrum[w]/(2*m)
    rz=np.eye(7)
    ra=np.zeros_like(active)
    for _ in range(epochs):
        ra=ra@empty+rz@active
        rz=rz@empty
    return shapes,rz,ra/epochs


def conditioning_log_probability(counts):
    n=sum(counts)
    return lgamma(n+1)-sum(lgamma(v+1) for v in counts)+sum(v*log(v/n) for v in counts if v)


def type_probability(counts,row_weights,transfers,tilt):
    shapes,empty,active=transfers
    p=np.zeros(5)
    p[0]=counts[0]/256
    for count,w in zip(counts[1:],row_weights):
        p[w]+=count/256
    probabilities=np.array([factorial(4)/prod(factorial(v) for v in Counter(shape).values())
                            *prod(p[w] for w in shape) for shape in shapes])
    assert abs(sum(probabilities)+p[0]**4-1)<1e-12
    region=empty*p[0]**4+np.einsum('i,ijk->jk',probabilities,active)
    moment=np.linalg.matrix_power(region,64)[0].sum()
    if moment<=0:
        raise ArithmeticError('moment underflow: use scaled or outward arithmetic')
    return min(0,(log(moment)+tilt*209715-conditioning_log_probability(counts))/log(2))


def fixed_support_regions(counts,row_weights,transfers):
    shapes,empty,active=transfers
    u=sum(counts[1:])
    p=np.zeros(5)
    for count,w in zip(counts[1:],row_weights):
        p[w]+=count/u
    result=[empty]
    for b in range(1,5):
        probabilities=np.array([factorial(4)/prod(factorial(v) for v in Counter(shape).values())
                                *prod(p[w] for w in shape if w)/comb(4,b)
                                if sum(w!=0 for w in shape)==b else 0 for shape in shapes])
        assert abs(sum(probabilities)-1)<1e-12
        result.append(np.einsum('i,ijk->jk',probabilities,active))
    return result


def fixed_support_update(best,counts,row_weights,transfers,tilt):
    """Condition exactly on support size and first active macroregion.

    Nonzero label counts still use an iid-conditioning upper bound. Binary64
    only; floor tiny moments conservatively for this diagnostic, not a cert.
    """
    u=sum(counts[1:])
    region=fixed_support_regions(counts,row_weights,transfers)
    current=np.ones((1,7))
    first=np.stack([r[0] for r in region[1:]],axis=1)
    transposed=[r.T*comb(4,b) for b,r in enumerate(region)]
    conditioning=conditioning_log_probability(counts[1:])
    for remaining in range(64):
        values=current@first
        for b in range(1,5):
            if 0<=u-b<current.shape[0]:
                moment=max(float(values[u-b,b-1]),1e-300)
                logbound=log(moment)+tilt*209715-log(comb(4*remaining,u-b))-conditioning
                best[remaining,b-1]=min(best[remaining,b-1],exp(min(0,logbound)))
        if remaining==63:
            break
        shifted=[current@r for r in transposed]
        count=current.shape[0]
        current=np.zeros((min(u+1,count+4),7))
        for b,row in enumerate(shifted):
            size=min(count,current.shape[0]-b)
            if size>0:
                current[b:b+size]+=row[:size]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weights',nargs=3,type=int,action='append')
    parser.add_argument('--fixed-support',action='store_true')
    args=parser.parse_args()
    self_test()
    spectrum=authenticated_caps()
    data=census()
    triples=args.weights or [(38,38,38),(38,38,40),(40,40,40),(48,48,48),(64,64,64),(128,128,128)]
    best={tuple(w):{v:0.0 for v in embeddings()} for w in triples}
    fixed={w:{v:np.ones((64,4)) for v in embeddings()} for w in best} if args.fixed_support else None
    for tilt in map(float,TILTS):
        transfers=region_transfers(data,tilt)
        for weights in best:
            counts=composition(weights,256)
            assert counts is not None and subspace_cap(weights,spectrum)>0
            for row_weights in best[weights]:
                value=type_probability(counts,row_weights,transfers,tilt)
                best[weights][row_weights]=min(best[weights][row_weights],value)
                if fixed is not None:
                    fixed_support_update(fixed[weights][row_weights],counts,row_weights,transfers,tilt)
        print('Diagnostic tilt completed',tilt,flush=True)
    for weights,bounds in best.items():
        counts=composition(weights,256)
        union=sum(embeddings()[v]*2**p for v,p in bounds.items())
        total=11+log2(subspace_cap(weights,spectrum))+log2(min(1,union))
        print('weights',weights,'composition',counts,'conditioning loss bits',-conditioning_log_probability(counts)/log(2),
              'single-type union log2',total,'worst row map',max(bounds,key=bounds.get),flush=True)
        if fixed is not None:
            u=sum(counts[1:])
            unions=sum(embeddings()[v]*p for v,p in fixed[weights].items())
            probability=sum(comb(4,b)*comb(4*l,u-b)/comb(256,u)*min(1,unions[l,b-1])
                            for l in range(64) for b in range(1,5) if 0<=u-b<=4*l)
            total=11+log2(subspace_cap(weights,spectrum))+log2(probability)
            print('  EXACT-SUPPORT diagnostic: single-type union log2',total,
                  'remaining conditioning loss bits',-conditioning_log_probability(counts[1:])/log(2),flush=True)
    print('Binary64 diagnostic only. Unscreened types, higher ranks, and multiple groups remain.')


if __name__=='__main__':
    main()
