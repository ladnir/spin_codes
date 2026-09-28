"""Pointwise Bernoulli-mixture domination of a shuffled BCH row.

Zero and all-one rows are kept as separate deterministic components.
Every other component is only a positive comparison measure, not a new
code distribution. Rational arithmetic checks all weight shells.
"""
from fractions import Fraction as Q
from math import comb,log


def envelope(caps,central,theta):
    central,theta=Q(central),Q(theta);n=len(caps)-1
    if n<1 or caps[0]!=1 or caps[-1]!=1 or central<=0 or not 0<theta<Q(1,2):
        raise ValueError('one zero/all-one row, positive central mass and a bias below one half required')
    density=[Q(c,comb(n,w)) for w,c in enumerate(caps)]
    flat=central/Q(2)**n
    tail=max(Q(0),max((density[w]-flat)/(theta**w*(1-theta)**(n-w)+(1-theta)**w*theta**(n-w))
                      for w in range(1,n)))
    components=[(Q(1),Q(0)),(Q(1),Q(1)),(central,Q(1,2)),(tail,theta),(tail,1-theta)]
    for w,a in enumerate(density):
        bound=sum((c*p**w*(1-p)**(n-w) for c,p in components),Q(0))
        if bound<a:raise ArithmeticError('row-mixture shell domination failed')
    return components


def pair_components(rows):
    """Return (name, coefficient, p00, p_single, p_double, active_label).

    The inactive pair is the zero/zero row component. All other terms
    dominate nonzero pair messages; even a comparison input that happens
    to be zero retains its active label in the occupancy accounting.
    """
    result=[]
    for i,(c,p) in enumerate(rows):
        for j in range(i,len(rows)):
            d,r=rows[j];coefficient=c*d*(1 if i==j else 2)
            if coefficient:
                result.append((f'{i},{j}',coefficient,(1-p)*(1-r),p*(1-r)+(1-p)*r,p*r,int(i+j!=0)))
    return result


if __name__=='__main__':
    import model
    from bch_joint_support import authenticated_caps
    caps=authenticated_caps()
    for exponent in (129,130,131,132):
        for theta in (Q(1,5),Q(1,4),Q(3,10),Q(1,3)):
            rows=envelope(caps,1<<exponent,theta);tail=rows[3][0]
            print('ROW ENVELOPE central',exponent,'theta',theta,'tail log2',
                  (log(tail.numerator)-log(tail.denominator))/log(2),flush=True)
