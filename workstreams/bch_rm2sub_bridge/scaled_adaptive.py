"""Per-degree binary exponents avoid dense-range underflow.

Each 3x3 matrix has its own shared power-of-two exponent. All operations
have a directed mode for certificates and a nearest mode for discovery.
"""
import math
import numpy as np
from flint import arb


def initial(region,outward=True):
    mantissas=[];exponents=[]
    for row in region:
        maximum=max(row)
        m,e=maximum.man_exp()
        exponent=int(e)+int(m).bit_length() if m else 0
        values=[float(v*arb(2)**(-exponent)) for v in row]
        mantissas.append(values);exponents.append(exponent)
    array=np.array(mantissas).reshape(-1,3,3)
    return (np.nextafter(array,np.inf) if outward else array),np.array(exponents,dtype=np.int64)


def matrices(mantissas,exponents,left,right,outward=True):
    current=mantissas.copy();powers=exponents.copy()
    def rounded(x):return np.nextafter(x,np.inf) if outward else x
    for q in range(1,len(current)):
        common=np.maximum(powers[:-1],powers[1:])
        a=rounded(np.ldexp(current[:-1],(powers[:-1]-common)[:,None,None]))
        b=rounded(np.ldexp(current[1:],(powers[1:]-common)[:,None,None]))
        updated=rounded(rounded(left[0]*a)+rounded(right[0]*b))
        for x,y in zip(left[1:],right[1:]):
            np.maximum(updated,rounded(rounded(x*a)+rounded(y*b)),out=updated)
        maximum=updated.max(axis=(1,2));assert (maximum>0).all() and np.isfinite(maximum).all()
        _,shift=np.frexp(maximum)
        current=rounded(np.ldexp(updated,-shift[:,None,None]));powers=common+shift
        yield q,current[0].copy(),int(powers[0])


def terminal_log(mantissas,exponents,ps,roots):
    import general_occupancy as general
    result=None
    for result in matrices(mantissas,exponents,roots*(1-ps),roots*ps,outward=False):pass
    assert result is not None
    _,matrix,exponent=result
    return general.log_power(matrix)+256*exponent*math.log(2)
