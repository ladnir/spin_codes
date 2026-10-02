"""Exact feedback distribution for uniform nonzero GF16 packets.

Choose j distinct packet positions uniformly, then choose each value
uniformly among 15 nonzero packets. Integer Walsh inversion counts every
feedback target; an explicit l1 bound prevents machine-integer overflow.
"""
import argparse
from math import comb
import numpy as np


def fwht(values):
    out=np.array(values,dtype=np.int64,copy=True)
    if out.ndim!=1 or not len(out) or len(out)&(len(out)-1):
        raise ValueError('power-of-two vector required')
    if sum(abs(int(x)) for x in out)>np.iinfo(np.int64).max:
        raise OverflowError('Walsh intermediate bound exceeds int64')
    step=1
    while step<len(out):
        block=out.reshape(-1,2*step);left=block[:,:step].copy();right=block[:,step:]
        block[:,:step]=left+right;block[:,step:]=left-right
        step*=2
    return out


def annihilators(columns,bits):
    if (type(bits) is not int or not 1<=bits<=24 or not columns or len(columns)%4
            or any(type(c) is not int or not 0<=c<1<<bits for c in columns)):
        raise ValueError('binary feedback columns and complete packet geometry required')
    chars=np.arange(1<<bits,dtype=np.uint32);result=np.zeros(len(chars),dtype=np.int64)
    for start in range(0,len(columns),4):
        mask=np.zeros(len(chars),dtype=np.uint8)
        for col in columns[start:start+4]:mask|=np.bitwise_count(chars&col)&1
        result+=mask==0
    return result


def distribution(annihilator_counts,windows,j):
    counts=np.asarray(annihilator_counts)
    if (type(windows) is not int or not 1<=windows<=32 or type(j) is not int or not 0<=j<=windows
            or counts.ndim!=1 or not len(counts) or len(counts)&(len(counts)-1)
            or not np.issubdtype(counts.dtype,np.integer) or np.any(counts<0) or np.any(counts>windows)):
        raise ValueError('valid character counts, window count, and occupancy required')
    # Fourier coefficient: [x^j] (1+15x)^c (1-x)^(windows-c).
    values=[sum(comb(c,k)*15**k*comb(windows-c,j-k)*(-1)**(j-k)
                for k in range(max(0,j-windows+c),min(c,j)+1)) for c in range(windows+1)]
    multiplicities=np.bincount(counts,minlength=windows+1)
    bound=sum(abs(x)*int(n) for x,n in zip(values,multiplicities))
    if bound>np.iinfo(np.int64).max:raise OverflowError('exact feedback census needs larger integers')
    transformed=fwht(np.array(values,dtype=np.int64)[counts])
    size=len(counts);denominator=comb(windows,j)*15**j
    if np.any(transformed%size):raise ArithmeticError('nonintegral inverse Walsh count')
    result=transformed//size
    if np.any(result<0) or sum(map(int,result))!=denominator:
        raise ArithmeticError('invalid feedback counts')
    return result,denominator


def census(columns,bits,max_packets=6):
    windows=len(columns)//4
    if type(max_packets) is not int or not 0<=max_packets<=windows:
        raise ValueError('occupancy census must fit the packet geometry')
    counts=annihilators(columns,bits);rows=[]
    for j in range(max_packets+1):
        values,denominator=distribution(counts,windows,j)
        target=int(values[1:].argmax())+1
        rows.append(dict(occupancy=j,zero=int(values[0]),nonzero_peak=int(values[target]),
                         peak_target=target,denominator=denominator))
    return rows


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--through',type=int,default=6)
    args=parser.parse_args()
    from packet_census import maps
    _,columns,_=maps()
    for row in census(columns,19,args.through):print(row,flush=True)
