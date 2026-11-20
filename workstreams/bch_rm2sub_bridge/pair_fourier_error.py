"""Point-dependent error enclosure for the fixed degree-128 pair evaluator.

Let R_j bound both exact and computed linear-form magnitudes. A is the
absolute monomial sum at R, and D is its sum of first derivatives. The
enclosure 2^-40*A + 2^-49*D + 2^-1000 covers arithmetic, input/root error,
and underflow respectively. All nonnegative enclosure operations round up.
"""
import numpy as np
from screen_pair_type_bound import SIGNS


def up(x):return np.nextafter(x,np.inf)


def evaluate(values,powers,counts):
    linear=SIGNS@values;squared=linear*linear;fourth=squared*squared
    table=np.empty((4,33,values.shape[1]),dtype=np.complex128);table[:,0]=1
    for n in range(1,33):table[:,n]=table[:,n-1]*fourth
    magnitude=up(np.sqrt(up(up(linear.real*linear.real)+up(linear.imag*linear.imag))))
    radius=up(magnitude+2**-49)
    fourth_upper=up(radius*radius);fourth_upper=up(fourth_upper*fourth_upper)
    positive=np.empty((4,33,values.shape[1]),dtype=np.float64);positive[:,0]=1
    for n in range(1,33):positive[:,n]=up(positive[:,n-1]*fourth_upper)
    value=np.zeros(values.shape[1],dtype=np.complex128)
    absolute=np.zeros(values.shape[1]);derivative=np.zeros(values.shape[1])
    for exponents,count in zip(powers,counts):
        value+=count*table[0,exponents[0]]*table[1,exponents[1]]*table[2,exponents[2]]*table[3,exponents[3]]
        term=np.full(values.shape[1],count);factor=np.zeros(values.shape[1])
        for j,e in enumerate(exponents):
            term=up(term*positive[j,e])
            if e:factor=up(factor+up((4*int(e))/radius[j]))
        absolute=up(absolute+term);derivative=up(derivative+up(term*factor))
    error=up(up(up(absolute*2**-40)+up(derivative*2**-49))+2**-1000)
    return value,error
