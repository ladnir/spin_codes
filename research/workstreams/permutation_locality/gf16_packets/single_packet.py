"""Exact dyadic upper census of one GF16 packet and shifted feedback."""
import numpy as np
from flint import arb


def census(images,columns,z,precision=44):
    if len(columns)%4 or not 0<z<=1 or precision<1:
        raise ValueError('complete packets, positive output weight, and precision required')
    scale=1<<precision;W=len(columns)//4;choices=15*W
    if choices*scale>=1<<64:raise ValueError('dyadic sums must fit uint64')
    powers=np.array([int((z**w*scale).upper().ceil().unique_fmpz()) for w in range(len(columns)+1)],dtype=np.uint64)
    low=np.array([v&((1<<64)-1) for v in images],dtype=np.uint64)
    high=np.array([v>>64 for v in images],dtype=np.uint64)
    states=np.arange(len(images),dtype=np.uint64)
    density=np.zeros(len(images),dtype=np.uint64);feedback=np.zeros(len(images),dtype=np.uint64)
    inputs=[]
    for w in range(W):
        for mask in range(1,16):
            syndrome=0
            for b in range(4):
                if mask>>b&1:syndrome^=columns[4*w+b]
            word=mask<<(4*w);lo=word&((1<<64)-1);hi=word>>64
            index=states^syndrome
            weights=np.bitwise_count(low[index]^lo)+np.bitwise_count(high[index]^hi)
            values=powers[weights].copy();values[syndrome]=0 # exclude entering state zero
            density+=values;feedback[syndrome]+=powers[mask.bit_count()]
            inputs.append((syndrome,word))
    denominator=choices*scale
    winner=int(density.argmax())
    for target in set((0,1,winner)):
        direct=sum(int(powers[(images[target^syndrome]^word).bit_count()])
                   for syndrome,word in inputs if target^syndrome)
        if direct!=int(density[target]):raise ArithmeticError('shifted density census mismatch')
    return dict(density=int(density[1:].max()),global_density=int(density.max()),
                cancellation=int(density[0]),feedback=int(feedback[1:].max()),
                zero=int(feedback[0]),denominator=denominator)


def actual(z):
    from group_moment import maps
    images,columns,_=maps()
    return census(images,columns,z)
