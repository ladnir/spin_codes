"""Joint feedback/output census for lazy returns to zero.

If the incoming state is s=CX, its output is (A C + I)X. We count feedback
and this output weight jointly, never multiply independent marginal bounds.
"""
from itertools import combinations,product
from math import comb
import numpy as np
from flint import arb
from feedback_exact import annihilators,distribution


def census(images,columns,bits,through=2,verbose=False):
    windows=len(columns)//4;size=1<<bits;width=len(columns);radix=width+1
    if (not 1<=windows<=32 or len(columns)%4 or len(images)!=size or images[0]!=0
            or type(through) is not int or not 0<=through<=min(4,windows)
            or any(not 0<=x<1<<width for x in images)
            or any(images[s]!=images[s^(s&-s)]^images[s&-s] for s in range(1,size))):
        raise ValueError('linear expansion, complete packets, and census degree at most four required')
    chars=annihilators(columns,bits)
    syndromes=np.zeros((windows,15),dtype=np.uint32)
    low=np.zeros((windows,15),dtype=np.uint64);high=low.copy()
    mask=(1<<64)-1
    for w in range(windows):
        for v in range(1,16):
            s=0
            for b,c in enumerate(columns[4*w:4*w+4]):
                if v>>b&1:s^=c
            y=images[s]^(v<<(4*w));syndromes[w,v-1]=s;low[w,v-1]=y&mask;high[w,v-1]=y>>64
    result=[]
    for j in range(through+1):
        denominator=comb(windows,j)*15**j
        if denominator>np.iinfo(np.uint32).max:raise OverflowError('joint census exceeds uint32 count bound')
        histogram=np.zeros(size*radix,dtype=np.uint32)
        choices=np.array(list(product(range(15),repeat=j)),dtype=np.intp).reshape(15**j,j)
        for positions in combinations(range(windows),j):
            s=np.zeros(15**j,dtype=np.uint32);lo=np.zeros(15**j,dtype=np.uint64);hi=lo.copy()
            for i,w in enumerate(positions):
                s^=syndromes[w,choices[:,i]];lo^=low[w,choices[:,i]];hi^=high[w,choices[:,i]]
            weights=np.bitwise_count(lo)+np.bitwise_count(hi)
            indices=s.astype(np.int64)*radix+weights
            np.add.at(histogram,indices,1)
        indices=np.flatnonzero(histogram);counts=histogram[indices].astype(np.uint64)
        targets,starts=np.unique(indices//radix,return_index=True)
        totals=np.add.reduceat(counts,starts)
        reference,total=distribution(chars,windows,j)
        if denominator!=total or sum(map(int,counts))!=total or not np.array_equal(np.flatnonzero(reference),targets):
            raise ArithmeticError('joint census support/mass disagrees with independent Walsh census')
        if not np.array_equal(reference[targets],totals):
            raise ArithmeticError('joint census feedback marginal disagrees with Walsh inversion')
        levels=np.array([images[int(s)].bit_count() for s in targets],dtype=np.int16)
        result.append(dict(occupancy=j,denominator=denominator,targets=targets,starts=starts,
                           weights=indices%radix,counts=counts,levels=levels,width=width))
        if verbose:print('JOINT RETURN census checked occupancy',j,'inputs',denominator,flush=True)
    return result


def outward(row,z):
    """Exact dyadic upper endpoints; every uint64 sum is overflow-checked."""
    if not 0<z<=1:raise ValueError('output weight in (0,1] required')
    denominator=row['denominator'];bits=min(44,62-denominator.bit_length());scale=1<<bits
    if denominator*scale>=1<<63:raise OverflowError('weighted count accumulation bound exceeded')
    powers=np.array([scale if w==0 else int((z**w*scale).upper().ceil().unique_fmpz())
                     for w in range(row['width']+1)],dtype=np.uint64)
    sums=np.add.reduceat(row['counts']*powers[row['weights']],row['starts'])
    active=row['targets']!=0;divisor=denominator*scale
    zero=int(sums[0]) if len(sums) and row['targets'][0]==0 else 0
    maximum=int(sums[active].max()) if active.any() else 0
    total=sum(map(int,sums[active]))
    by_level={int(level):sum(map(int,sums[(row['levels']==level)&active])) for level in set(row['levels'][active])}
    level_maxima={int(level):int(sums[(row['levels']==level)&active].max()) for level in set(row['levels'][active])}
    return dict(zero=arb(zero)/divisor,maximum=arb(maximum)/divisor,total=arb(total)/divisor,
                classes={level:arb(value)/divisor for level,value in by_level.items()},
                class_maxima={level:arb(value)/divisor for level,value in level_maxima.items()},
                numerator=maximum,denominator=divisor)


def actual_census(data,through):
    """Regenerate and cross-check the integer census for the actual maps."""
    from group_moment import maps
    images,columns,_=maps()
    if columns!=data['columns'] or len(images)!=(1<<data['bits']):
        raise ArithmeticError('joint-return census map mismatch')
    return census(images,columns,data['bits'],through,verbose=True)


def refine_class_returns(data,local,checked_census,z):
    """Tighten zero entries; apply before changing the refreshed U entries."""
    from occupancy_kernel import up
    alpha=arb(2)**-data['updates'];L=(1<<data['bits'])-1
    result=[matrix*arb(1) for matrix in local]
    for row in checked_census:
        j=row['occupancy']
        if type(j) is not int or not 0<=j<len(result):raise ValueError('joint census occupancy outside local geometry')
        matrix=result[j];moments=outward(row,z)
        caps=[moments['maximum'],moments['total']/L,
              *(moments['class_maxima'].get(int(level),arb(0)) for level in data['birth_class_levels'])]
        for source,cap in enumerate(caps,1):
            matrix[source,0]=min(matrix[source,0],up(alpha*cap+matrix[source,2]/L))
    return result


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--through',type=int,default=2)
    parser.add_argument('--tilt',default='.064')
    args=parser.parse_args()
    from packet_census import maps
    images,columns,_=maps();data=census(images,columns,19,args.through)
    for row in data:
        moment=outward(row,(-arb(args.tilt)).exp())
        print('GF16 joint return',row['occupancy'],'inputs',row['denominator'],'weighted maximum',moment['maximum'],
              'weighted sum',moment['total'],'zero',moment['zero'],flush=True)
