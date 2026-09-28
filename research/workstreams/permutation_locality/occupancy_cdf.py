"""Cumulative-count bounds for coordinatewise decreasing support weights.

CDF increments are a dominating measure, not shell upper bounds. For each
coordinate, summation by parts bounds a nonnegative decreasing function by
nonnegative multiples of the upper CDF. Repeating along the coordinates is
valid without a sign condition on mixed differences.
"""
from itertools import product
import numpy as np


def differences(caps):
    values=[int(v) for v in caps]
    assert values[0]>=0
    assert all(a<=b for a,b in zip(values,values[1:]))
    return [values[0]]+[b-a for a,b in zip(values,values[1:])]


def majorant(values):
    """Least coordinatewise decreasing majorant on a rectangular grid."""
    result=values.copy()
    for axis in range(result.ndim):
        result=np.flip(np.maximum.accumulate(np.flip(result,axis),axis=axis),axis)
    return result


def tensor_sum(weights,measures):
    return sum(weights[index]*np.prod([measures[j][u] for j,u in enumerate(index)],dtype=object)
               for index in product(*(range(len(m)) for m in measures)))


def self_test():
    checks=0
    # Every four-point integer distribution with mass <=8; move positive
    # mass to earlier positions to generate dominating CDFs.
    for counts in product(range(3),repeat=4):
        dominating=list(counts)
        for j in range(1,4):
            if dominating[j]:
                dominating[j]-=1
                dominating[j-1]+=1
        caps=list(np.cumsum(dominating))
        assert differences(caps)==dominating
        assert all(a<=b for a,b in zip(np.cumsum(counts),caps))
        for q in (1,2,3):
            raw=np.array([((sum((j+1)*u*u for j,u in enumerate(index))+3)%11)
                          for index in product(range(4),repeat=q)],dtype=object).reshape((4,)*q)
            decreasing=majorant(raw)
            assert np.all(decreasing>=raw)
            assert np.array_equal(majorant(decreasing),decreasing)
            assert tensor_sum(decreasing,[counts]*q)<=tensor_sum(decreasing,[dominating]*q)
            checks+=1
    # Decreasing does not require nonnegative mixed differences.
    values=np.array([[1,1],[1,0]],dtype=object)
    assert values[0,0]-values[0,1]-values[1,0]+values[1,1]<0
    assert tensor_sum(values,[[1,1],[1,1]])<=tensor_sum(values,[[2,0],[2,0]])
    # A truncated lower cube is handled by extending the function by zero.
    values=np.zeros((4,4),dtype=object)
    values[:2,:2]=[[4,2],[3,1]]
    assert np.array_equal(majorant(values),values)
    print('CDF domination tests passed:',checks,'integer tensor inequalities',flush=True)


if __name__=='__main__':
    self_test()
