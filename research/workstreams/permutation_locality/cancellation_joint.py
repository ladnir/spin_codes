"""Joint feedback/output enumeration for exact lazy-return coefficients.

Enumerates every input using at most three distinct four-bit windows.
No floating-point counts or sampled states enter the enumeration.
"""
from collections import Counter
from itertools import combinations,combinations_with_replacement,product
from math import comb,factorial,prod
import argparse
import numpy as np
from flint import arb,ctx
from group_moment import maps
from group_rank_one_verify import up
from occupancy_memory import C,F,U,Z


def census(maximum=3,check_old=True):
    assert maximum in (1,2,3)
    images,columns,spectrum=maps();levels=[0]+sorted(spectrum)
    low=np.array([w&((1<<64)-1) for w in images],dtype=np.uint64)
    high=np.array([w>>64 for w in images],dtype=np.uint64)
    class_index=np.searchsorted(levels,np.bitwise_count(low)+np.bitwise_count(high))
    feedback=np.zeros((32,15),dtype=np.uint32)
    inputs_lo=np.zeros((32,15),dtype=np.uint64);inputs_hi=inputs_lo.copy()
    for window in range(32):
        for bits in range(1,16):
            value=0
            for bit in range(4):
                if bits>>bit&1:value^=columns[4*window+bit]
            feedback[window,bits-1]=value
            word=bits<<(4*window)
            inputs_lo[window,bits-1]=word&((1<<64)-1);inputs_hi[window,bits-1]=word>>64
    fresh=np.zeros((4,1<<19),dtype=np.uint8)
    for window in range(32):
        for bits in range(1,16):fresh[bits.bit_count()-1,feedback[window,bits-1]]+=1
    assert int(fresh.max())==1 and not fresh[:,0].any()
    assert list(map(int,fresh.sum(axis=1)))==[32*comb(4,b) for b in range(1,5)]
    result={}
    for j in range(1,maximum+1):
        shapes=list(combinations_with_replacement(range(1,5),j))
        indices=np.indices((15,)*j).reshape(j,-1)
        lane_weights=np.array([b.bit_count() for b in range(1,16)])
        shape_ids=np.array([shapes.index(tuple(sorted(lane_weights[indices[:,i]]))) for i in range(indices.shape[1])])
        joint=np.zeros(len(shapes)*len(levels)*129,dtype=np.int64)
        fresh_hist=np.zeros((4,len(shapes)*129),dtype=np.int64)
        assert comb(32,j)*15**j<1<<63
        for windows in combinations(range(32),j):
            syndrome=np.zeros(indices.shape[1],dtype=np.uint32)
            word_lo=np.zeros(indices.shape[1],dtype=np.uint64);word_hi=word_lo.copy()
            for slot,window in enumerate(windows):
                syndrome^=feedback[window,indices[slot]]
                word_lo^=inputs_lo[window,indices[slot]];word_hi^=inputs_hi[window,indices[slot]]
            weight=np.bitwise_count(low[syndrome]^word_lo)+np.bitwise_count(high[syndrome]^word_hi)
            ids=(shape_ids*len(levels)+class_index[syndrome])*129+weight
            local=np.bincount(ids,minlength=len(joint));joint+=local
            fresh_ids=shape_ids*129+weight
            for a in range(4):
                fresh_hist[a]+=np.bincount(fresh_ids[fresh[a,syndrome]!=0],minlength=fresh_hist.shape[1])
            if windows in (tuple(range(j)),tuple(range(32-j,32))):
                reference=Counter()
                for masks in product(range(1,16),repeat=j):
                    word=0;s=0
                    for window,bits in zip(windows,masks):
                        word^=bits<<(4*window)
                        for bit in range(4):
                            if bits>>bit&1:s^=columns[4*window+bit]
                    shape=tuple(sorted(b.bit_count() for b in masks))
                    index=(shapes.index(shape)*len(levels)+levels.index(images[s].bit_count()))*129+(images[s]^word).bit_count()
                    reference[index]+=1
                assert all(int(n)==reference[i] for i,n in enumerate(local))
        for i,shape in enumerate(shapes):
            den=comb(32,j)*factorial(j)//prod(factorial(n) for n in Counter(shape).values())*prod(comb(4,b) for b in shape)
            rows={v:[int(x) for x in joint[(i*len(levels)+k)*129:(i*len(levels)+k+1)*129]] for k,v in enumerate(levels)}
            ff=[[int(x) for x in fresh_hist[a,i*129:(i+1)*129]] for a in range(4)]
            assert sum(sum(row) for row in rows.values())==den
            assert all(sum(row)<=den for row in ff)
            assert all(not n or abs(v-sum(shape))<=w<=v+sum(shape)
                       for v,row in rows.items() for w,n in enumerate(row))
            result[shape]=rows,ff,den
        print('Joint cancellation census:',j,'windows;',len(shapes),'shapes; all input totals and direct samples checked',flush=True)
    if check_old and maximum>=2:
        from zero_moment import census as prior_census
        prior=prior_census(maximum)
        for shape,(rows,_,den) in result.items():
            if len(shape)<2:continue
            hist,old_den=prior[len(shape)][shape]
            assert all(sum(row[w] for row in rows.values())*old_den==hist[w]*den for w in range(129))
        print('Every joint output marginal matches the earlier independent cancellation census',flush=True)
    return result,spectrum


def lazy_coefficients(data,tilt,penalty,rounds=2,input_penalty=1,odd_penalty=1):
    rows,spectrum=data;levels=sorted(spectrum)
    powers=[(-arb(tilt)*w).exp() for w in range(129)]
    alpha=arb(2)**(-rounds);maxima={}
    for shape,(joint,fresh,den) in rows.items():
        scale=alpha*arb(penalty)**shape.count(4)*arb(input_penalty)**sum(shape)*arb(odd_penalty)**sum(b%2 for b in shape)
        moments={v:sum((n*p for n,p in zip(hist,powers)),arb(0))/den for v,hist in joint.items() if v}
        values={C:sum(moments.values(),arb(0)),
                F:max(sum((n*p for n,p in zip(hist,powers)),arb(0))/(den*32*comb(4,a)) for a,hist in enumerate(fresh,1))}
        values.update({U+i:moments[v]/spectrum[v] for i,v in enumerate(levels)})
        for source,value in values.items():
            key=len(shape),source
            maxima[key]=max(maxima.get(key,arb(0)),up(scale*value))
    return maxima


def check_fresh(data,prepared):
    from fractions import Fraction
    from fresh_collision import probability
    tables={1:{},2:{},**prepared[0][2]};checks=0
    for shape,(_,fresh,den) in data[0].items():
        for a,hist in enumerate(fresh,1):
            observed=Fraction(sum(hist),den*32*comb(4,a))
            expected=(Fraction(int(a==shape[0]),32*comb(4,a)) if len(shape)==1
                      else probability(tables,a,shape))
            assert observed==expected,(a,shape,observed,expected)
            checks+=1
    print('Joint fresh marginals:',checks,'independent overlap/kernel identities passed',flush=True)


def check_feedback(data,feedback):
    checks=0
    for shape,(joint,_,den) in data[0].items():
        zero,_,old_den,classes=feedback[shape]
        assert sum(joint[0])*old_den==zero*den
        assert all(sum(joint[v])*old_den==n*den for v,n in classes.items())
        checks+=1
    print('Joint feedback marginals:',checks,'independent Walsh-distribution identities passed',flush=True)


def refine(base,data,tilt,penalty,rounds=2,input_penalty=1,odd_penalty=1):
    _,spectrum=data;changes=0
    for (j,source),lazy in lazy_coefficients(data,tilt,penalty,rounds,input_penalty,odd_penalty).items():
        if j>=len(base):continue
        value=lazy
        if source!=C:
            refresh=min(up(base[j][source,U+i]/spectrum[v]) for i,v in enumerate(sorted(spectrum)))
            value=up(value+refresh)
        if value<base[j][source,Z]:changes+=1
        base[j][source,Z]=min(base[j][source,Z],value)
    print('Joint cancellation refinement:',changes,'zero-return coefficients tightened',flush=True)
    return base


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,choices=(1,2,3),default=3)
    parser.add_argument('--precision',type=int,default=192)
    args=parser.parse_args();data=census(args.maximum);ctx.prec=args.precision
    from occupancy_model import local_data
    from occupancy_memory import prepare
    check_fresh(data,prepare(local_data(4)))
    values=lazy_coefficients(data,'.032','.75')
    for (j,source),value in sorted(values.items()):print('LAZY',j,source,value,flush=True)
