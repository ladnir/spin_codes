"""Exact XOR convolutions for two-window low-expansion-weight probabilities."""
from fractions import Fraction
from itertools import combinations_with_replacement
from math import comb
import numpy as np
from fresh_collision import pair_completions


def hadamard(values):
    out=np.array(values,dtype=np.int64,copy=True)
    width=1
    while width<len(out):
        blocks=out.reshape(-1,2*width)
        left=blocks[:,:width].copy();right=blocks[:,width:].copy()
        blocks[:,:width]=left+right;blocks[:,width:]=left-right
        width*=2
    return out


def census(inputs,prepared):
    low,high,words,spectrum=inputs
    weights=np.bitwise_count(low)+np.bitwise_count(high)
    n=len(low)
    transforms={0:np.full(n,32,dtype=np.int64)}
    for b,entries in words.items():
        counts=np.zeros(n,dtype=np.int64)
        for syndrome,_,_ in entries:counts[syndrome]+=1
        transforms[b]=hadamard(counts)
        assert np.array_equal(hadamard(transforms[b]),counts*n)
    indicators={cut:hadamard(((weights>0)&(weights<=cut)).astype(np.int64)) for cut in (48,56)}
    result={}
    for a,b in combinations_with_replacement(range(1,5),2):
        transformed=transforms[a]*transforms[b]
        for d in range(5):transformed-=pair_completions(a,b,d)*transforms[d]
        denominator=32*31*comb(4,a)*comb(4,b)
        assert int(np.abs(transformed).max())<=denominator
        counts=hadamard(transformed)
        assert not np.any(counts%n)
        counts//=n
        reference=prepared[0][1][1][(a,),(b,)]
        assert int(counts.sum())==denominator and int(counts.min())>=0
        assert int(counts[0])==reference[1] and int(counts[1:].max())==reference[2]
        # Independent direct enumeration of every ordered distinct-window pair.
        sa=np.array([s for s,_,_ in words[a]],dtype=np.int64)
        sb=np.array([s for s,_,_ in words[b]],dtype=np.int64)
        wa=np.array([((lo|(hi<<64)).bit_length()-1)//4 for _,lo,hi in words[a]])
        wb=np.array([((lo|(hi<<64)).bit_length()-1)//4 for _,lo,hi in words[b]])
        direct=np.bincount((sa[:,None]^sb[None,:])[wa[:,None]!=wb[None,:]],minlength=n)
        assert np.array_equal(counts,direct)
        records={}
        for cutoff,indicator in indicators.items():
            cardinality=sum(size for weight,size in spectrum.items() if weight<=cutoff)
            # This also bounds every intermediate of the inverse transform.
            assert n*denominator*cardinality<1<<63
            hits=hadamard(transformed*indicator)
            assert not np.any(hits%n)
            hits//=n
            assert 0<=int(hits.min())<=int(hits.max())<=denominator
            assert int(hits.sum())==denominator*cardinality
            assert int(hits[0])==sum(sum(hist.values()) for v,hist in reference[3].items() if v<=cutoff)
            for state in (0,1,17,int(hits.argmax())):
                selected=(weights[np.arange(n)^state]>0)&(weights[np.arange(n)^state]<=cutoff)
                assert int(counts[selected].sum())==int(hits[state])
            by_weight={v:Fraction(int(hits[weights==v].max()),denominator) for v in spectrum}
            density=Fraction(int(counts[weights==cutoff].max()),denominator)
            records[cutoff]=(Fraction(int(hits.max()),denominator),Fraction(int(hits[0]),denominator),by_weight,density)
        result[a,b]=records
        print('Pair-tail exact XOR counts',a,b,[(cut,str(p[0])) for cut,p in records.items()],flush=True)
    return result
