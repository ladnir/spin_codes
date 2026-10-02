"""A fixed-operator chord bound for an entire mean-activity interval.

Let T_j be nonnegative upper transition matrices conditional on exactly j
active packets among W positions, and let E be the number of steps. With
P=WE, the weighted transition moment equals

    (1-m)^P H(m/(x(1-m))),

where H has nonnegative coefficients and degree at most P. Here x is the
base activity tilt, not an additional randomness parameter. The function
log H(exp(u)) is convex. Its endpoint chord bounds the whole interval.
Combining that chord with the outer linear penalty gives a concave scalar
maximization, without separately maximizing the two input probabilities.

This module constructs a new candidate bound only. It does not modify the
existing certificate engine or claim coverage beyond the supplied cell.
"""
from fractions import Fraction as Q
from math import comb

from flint import arb, arb_mat, ctx


def aq(value):
    value=Q(value)
    return arb(value.numerator)/arb(value.denominator)


def upper_fraction(value):
    m,e=value.upper().man_exp()
    return Q(int(m))*Q(2)**int(e)


def compact_endpoint(value):
    """Keep outward dyadics compact even when their exponent is enormous."""
    value=Q(value)
    n,d=value.numerator,value.denominator
    if d & (d-1):
        # Only toy callers supply nondyadic endpoints; production endpoints
        # come directly from upper_fraction and always take the first form.
        return dict(rational=[str(n),str(d)])
    exponent=1-d.bit_length()
    if n:
        shift=(abs(n)&-abs(n)).bit_length()-1
        n>>=shift;exponent+=shift
    return dict(dyadic=[n,exponent])


def concave_maximum(a, degree, penalty, cell, steps=80):
    """Upper-bound max a log(m)+(P-a)log(1-m)-b*m exactly in scope."""
    a,penalty=Q(a),Q(penalty)
    lo,hi=map(Q,cell)
    if (type(degree) is not int or degree<1 or not 0<=a<=degree
            or not 0<lo<=hi<1 or type(steps) is not int or steps<1):
        raise ValueError('valid degree, slope, interior mean cell, and bisection budget required')
    def derivative(m):
        return a/m-(degree-a)/(1-m)-penalty
    def function(m):
        return aq(a)*aq(m).log()+aq(degree-a)*aq(1-m).log()-aq(penalty*m)
    if derivative(lo)<=0:
        bracket=(lo,lo)
    elif derivative(hi)>=0:
        bracket=(hi,hi)
    else:
        left,right=lo,hi
        for _ in range(steps):
            mid=(left+right)/2
            if derivative(mid)>0:left=mid
            else:right=mid
        bracket=(left,right)
    mid=sum(bracket)/2
    # A tangent is an upper bound for this concave function. The unique
    # maximizer is inside the exact derivative bracket, including endpoints.
    correction=max(derivative(mid)*(bracket[0]-mid),derivative(mid)*(bracket[1]-mid))
    value=function(mid)+aq(correction)
    return upper_fraction(value),dict(bracket=list(map(str,bracket)),at=str(mid),
                                     tangent_correction=str(correction))


def chord_from_endpoints(hlo, hhi, degree, cell, base_tilt, penalty=0):
    """Bound a fixed nonnegative polynomial from valid endpoint upper bounds.

    hlo and hhi must bound the SAME polynomial H at the declared odds.
    They are supplied by transition_endpoints in the production wrapper.
    Clipping the outward chord slope to P is safe by the independent degree
    bound H(v)<=H(vlo)*(v/vlo)^P. Clipping a negative slope to zero weakens
    the chord because v>=vlo. Thus the scalar majorant is always concave.
    """
    hlo,hhi,base_tilt,penalty=map(Q,(hlo,hhi,base_tilt,penalty))
    lo,hi=map(Q,cell)
    if (hlo<=0 or hhi<=0 or type(degree) is not int or degree<1
            or not 0<lo<=hi<1 or base_tilt<=0):
        raise ValueError('positive endpoints and tilt, degree, and interior mean cell required')
    vlo=lo/(base_tilt*(1-lo));vhi=hi/(base_tilt*(1-hi))
    if lo==hi:
        value=aq(hlo).log()+degree*aq(1-lo).log()-aq(penalty*lo)
        return upper_fraction(value),dict(slope=None,point=str(lo))
    slope_ball=(aq(hhi).log()-aq(hlo).log())/aq(vhi/vlo).log()
    slope=min(Q(degree),max(Q(0),upper_fraction(slope_ball)))
    peak,witness=concave_maximum(slope,degree,penalty,(lo,hi))
    # x*vlo=lo/(1-lo), so this anchor avoids an unnecessary cancellation.
    value=aq(hlo).log()-aq(slope)*aq(lo/(1-lo)).log()+aq(peak)
    return upper_fraction(value),dict(slope=str(slope),maximum=witness,
        endpoint_uppers=[compact_endpoint(hlo),compact_endpoint(hhi)],odds=[str(vlo),str(vhi)])


def transition_endpoints(local, epochs, cell, base_tilt):
    """Evaluate one fixed set of conditional operators at both endpoint odds."""
    lo,hi=map(Q,cell);base_tilt=Q(base_tilt)
    if (not isinstance(local,(list,tuple)) or len(local)<2
            or type(epochs) is not int or epochs<1 or not 0<lo<=hi<1 or base_tilt<=0):
        raise ValueError('fixed conditional operators, epochs, and interior cell required')
    size=local[0].nrows();windows=len(local)-1
    if (size<1 or any(m.nrows()!=size or m.ncols()!=size for m in local)
            or any(not m[i,j]>=0 for m in local for i in range(size) for j in range(size))):
        raise ValueError('equal nonnegative square conditional matrices required')
    # Freeze each operator to exact outward dyadic endpoints before the two
    # evaluations. This creates one well-defined nonnegative polynomial.
    fixed=[arb_mat([[aq(upper_fraction(m[i,j])) for j in range(size)]
                    for i in range(size)]) for m in local]
    values=[]
    for mean in (lo,hi):
        odds=mean/(base_tilt*(1-mean))
        matrix=sum((m*aq(comb(windows,j)*odds**j) for j,m in enumerate(fixed)),arb_mat(size,size))
        power=matrix**epochs
        value=sum((power[0,j] for j in range(size)),arb(0))
        if not value>0:raise ArithmeticError('strictly positive endpoint moment required')
        values.append(upper_fraction(value))
    return tuple(values)


def outward(model, cell, witness):
    """New cell bound using existing outer and shuffle bounds with fixed T_j."""
    import occupancy_birth_classes
    import scalar_cover as sc
    if (not model.variance_shuffle or model.inner.__name__!='birth_classes'
            or Q(witness['tilt'])!=model.tilt or len(witness['parameters'])!=3
            or 'variance_dual' not in witness):
        raise ValueError('base-tilt birth-class witness with variance dual required')
    lam,eta,mu=map(Q,witness['parameters'])
    if lam<=0 or mu<0:raise ValueError('positive output tilt and nonnegative occupancy dual required')
    precision=ctx.prec
    local=occupancy_birth_classes.outward(model.data,lam)
    hlo,hhi=transition_endpoints(local,sc.EPOCHS,cell,model.tilt)
    log_inner,details=chord_from_endpoints(hlo,hhi,sc.PACKETS,cell,model.tilt,sc.G*eta)
    cs,_,_=model.family(model.tilt)
    mass=sum((aq(c)*aq(eta*f+mu*a).exp()
              for c,f,a in zip(cs,model.features,model.active)),arb(0))
    loss=model.comparison_loss(cell,model.tilt,witness['variance_dual'])
    value=(aq(log_inner)+sc.G*mass.log()-aq(mu*model.q_min)
           +sc.REGIONS*aq(loss).log()+aq(lam*model.threshold)).exp()
    if ctx.prec!=precision:raise ArithmeticError('precision changed during chord verification')
    return upper_fraction(value),details
