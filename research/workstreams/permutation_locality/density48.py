"""Separate an auxiliary pointwise density cap on expansion-weight-48 states."""
from collections import Counter
from itertools import combinations_with_replacement
import numpy as np
from flint import arb,arb_mat
from occupancy_memory import Z,F,C,U,rounded
from occupancy_window_average import dyadic_powers
from group_rank_one_verify import up

D48=11


def census(inputs,prepared):
    low,high,words,spectrum=inputs
    weights=np.bitwise_count(low)+np.bitwise_count(high)
    atoms=prepared[0][0][3];records={}
    for a in range(1,5):
        for b,entries in words.items():
            hist={}
            for state,count in atoms[a,].items():
                if not count:continue
                for syndrome,lo,hi in entries:
                    target=state^syndrome
                    if weights[target]!=48:continue
                    output=(int(low[state])^lo).bit_count()+(int(high[state])^hi).bit_count()
                    hist.setdefault(target,Counter())[output]+=count
            denominator=sum(atoms[a,].values())*len(entries)
            records[a,b]=hist,denominator
    print('Fresh-to-density48 exact weighted fibers:',sum(len(h) for h,_ in records.values()),flush=True)
    return records


def averages(inputs,records,tilt):
    low,high,words,spectrum=inputs
    weights=np.bitwise_count(low)+np.bitwise_count(high)
    selected=np.flatnonzero(weights==48)
    powers,scale=dyadic_powers(tilt)
    result={}
    for b,entries in words.items():
        sums=np.zeros(len(selected),dtype=np.uint64)
        uniform={v:np.zeros(len(selected),dtype=np.uint64) for v in spectrum}
        for syndrome,lo,hi in entries:
            incoming=selected^syndrome
            output=np.bitwise_count(low[incoming]^np.uint64(lo))+np.bitwise_count(high[incoming]^np.uint64(hi))
            values=powers[output];sums+=values
            for v in spectrum:uniform[v]+=values*(weights[incoming]==v)
        denominator=len(entries)*scale
        maximum=int(sums.max())
        winner=int(selected[int(sums.argmax())])
        check=sum(int(powers[(int(low[winner^s])^lo).bit_count()+(int(high[winner^s])^hi).bit_count()]) for s,lo,hi in entries)
        assert maximum==check
        fresh=[]
        for a in range(1,5):
            fibers,count=records[a,b]
            numerator=max((sum(c*int(powers[w]) for w,c in h.items()) for h in fibers.values()),default=0)
            fresh.append(up(arb(numerator)/(count*scale)))
        result[b]=(up(arb(maximum)/denominator),max(fresh),
                   {v:up(arb(int(row.max()))/(denominator*spectrum[v])) for v,row in uniform.items()})
    print('Density48 restricted maxima and independent maximizing-state replay passed',flush=True)
    return result


def lift(base,values,tilt,penalty,pairs=None):
    result=[]
    for j,old in enumerate(base):
        assert old.nrows()==11
        t=arb_mat([[old[a,b] if a<11 and b<11 else arb(0) for b in range(12)] for a in range(12)])
        if not j:
            f48=(-arb(tilt)*48).exp();f56=(-arb(tilt)*56).exp()
            t[C,C]=up(f56/2);t[D48,C]=up((f48-f56)/2)
            t[D48,D48]=up(f48/2)
        else:
            for source in range(11):t[source,D48]=old[source,C]
            if j==1:
                for source,index in ((C,0),(F,1)):
                    t[source,D48]=min(t[source,D48],max(up(value[index]*arb(penalty)**int(b==4)/2) for b,value in values.items()))
                for k,v in enumerate((48,56,64,72,80)):
                    t[U+k,D48]=min(t[U+k,D48],max(up(value[2][v]*arb(penalty)**int(b==4)/2) for b,value in values.items()))
            elif j==2 and pairs is not None:
                candidates=[]
                for a,b in combinations_with_replacement(range(1,5),2):
                    density=pairs[a,b][48][3]
                    f=(-arb(tilt)*(a+b)).exp()*arb(penalty)**((a==4)+(b==4))
                    candidates.append(up(f*arb(density.numerator)/density.denominator))
                t[Z,D48]=min(t[Z,D48],max(candidates))
        result.append(rounded(t))
    return result
