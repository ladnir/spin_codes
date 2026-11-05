"""Directed switch sensitivity on a Gaussian ball, not a distance certificate.

The global constant remains necessary in the far tail. Here tilted pair
probabilities retain their actual three outcomes and common tilt variables.
All computations are file-free; the optional integration uses a fixed center.
"""
import math
import shutil
import sys
import time
from decimal import Decimal as D, getcontext
from fractions import Fraction as F
from itertools import product

from random_pairing_flow_envelope_verify import I, arithmetic_tests
from random_pairing_flow_quadratic import (
    TYPES, LAW, OUTPUT, LOSS, RETUNED_LOG_WEIGHT_NUMERATORS,
)
from random_pairing_switch_lipschitz import switch_keys


def threshold_mean_bounds(outcomes, probabilities):
    """E[v_i]=Pr[v_i>=1]+Pr[v_i>=2], with the probability simplex retained."""
    result=[]
    for i in range(5):
        value=I()
        for threshold in (1,2):
            inside=[j for j in range(3) if outcomes[j][i]>=threshold]
            outside=[j for j in range(3) if j not in inside]
            a=sum((probabilities[j] for j in inside),I())
            b=1-sum((probabilities[j] for j in outside),I())
            value+=I(max(D(0),a.lo,b.lo),min(D(1),a.hi,b.hi))
        result.append(value)
    return result


def pair_mean_bounds(radius, log_weights=RETUNED_LOG_WEIGHT_NUMERATORS, center=None):
    """Uniform for ||Z||<=radius and either one center or all t>=0,sum(t)=512.

    Pair outcome j has unnormalized weight
      LAW[j] exp(.1 LOSS_j+.2j + <u+.16t/512+beta Z, OUTPUT_j>).
    The four-dimensional Gaussian embeds isometrically in the zero-sum
    subspace of R^5. Outcome differences have zero sum. Cauchy--Schwarz
    therefore bounds every log weight ratio by beta*radius*||difference||.
    With no fixed center, its contribution is between .16 min(d) and
    .16 max(d). We bound each probability by its ratio to all other weights.
    """
    radius=I(radius);assert radius.lo>=0
    beta=(I(4)/25/512).sqrt();u=[I(x)/10**9 for x in log_weights]
    if center is not None:
        center=[I(x) for x in center]
        assert len(center)==5 and min(x.lo for x in center)>=0
        total=sum(center,I());assert total.lo<=512<=total.hi
    result={}
    for a,b in TYPES:
        outcomes=[OUTPUT[a,b,j] for j in range(3)]
        logs=[(I(LAW[j].numerator)/LAW[j].denominator).ln()+I(LOSS[a,b,j])/10+I(j)/5
              for j in range(3)]
        probabilities=[]
        for j in range(3):
            denominator=I(1)
            for k in range(3):
                if k==j:continue
                difference=[x-y for x,y in zip(outcomes[k],outcomes[j])]
                assert sum(difference)==0
                nominal=logs[k]-logs[j]+sum((x*y for x,y in zip(u,difference)),I())
                if center is None:nominal+=(I(4)/25)*I(min(difference),max(difference))
                else:nominal+=(I(4)/25/512)*sum((x*y for x,y in zip(center,difference)),I())
                deviation=beta*radius*I(sum(x*x for x in difference)).sqrt()
                denominator+=(nominal+I(-deviation.hi,deviation.hi)).exp()
            probabilities.append(1/denominator)
        result[a,b]=threshold_mean_bounds(outcomes,probabilities)
    return result


def switch_squared_bounds(radius, log_weights=RETUNED_LOG_WEIGHT_NUMERATORS, center=None):
    means=pair_mean_bounds(radius,log_weights,center);result={}
    for key in switch_keys():
        gradient=[sum((coefficient*means[pair][i] for pair,coefficient in key),I()) for i in range(5)]
        # Projection onto the contrast space cannot increase the norm.
        result[key]=sum((x.sq() for x in gradient),I()).upper()
    return result


def subset_squared_bound(bounds, subset):
    allowed=set(subset)
    return I(max((value.hi for key,value in bounds.items()
                  if all(a in allowed and b in allowed for (a,b),_ in key)),default=D(0)))


def exact_tests(deadline=math.inf):
    checks=arithmetic_tests()
    # Exact rational probability boxes include all their simplex points.
    for a,b in TYPES:
        assert time.monotonic()<deadline,'switch-ball test deadline'
        outcomes=[OUTPUT[a,b,j] for j in range(3)]
        for numerators in product(range(1,4),repeat=3):
            total=sum(numerators);p=[F(x,total) for x in numerators]
            probabilities=[I(max(D(0),(I(x.numerator)/x.denominator-I(1)/10).lo),
                             min(D(1),(I(x.numerator)/x.denominator+I(1)/10).hi)) for x in p]
            bounds=threshold_mean_bounds(outcomes,probabilities)
            for i,value in enumerate(bounds):
                actual=sum(p[j]*outcomes[j][i] for j in range(3))
                # Do not convert outward-rounded zero's subnormal Decimal
                # endpoint into a million-digit Fraction.
                if actual==0:assert value.lo<=0<=value.hi
                else:
                    assert value.lo<=0 or F(value.lo)<=actual
                    assert actual<=F(value.hi)
                checks+=1
    # Radius zero at a fixed center must enclose direct evaluation.
    for center in ((0,0,512,0,0),(100,100,112,100,100)):
        means=pair_mean_bounds(0,center=center)
        u=[I(x)/10**9+(I(4)/25)*c/512 for x,c in zip(RETUNED_LOG_WEIGHT_NUMERATORS,center)]
        for a,b in TYPES:
            weights=[(I(LAW[j].numerator)/LAW[j].denominator)*
                     (I(LOSS[a,b,j])/10+I(j)/5+sum((x*y for x,y in zip(u,OUTPUT[a,b,j])),I())).exp()
                     for j in range(3)]
            for i,actual in enumerate(means[a,b]):
                direct=sum((w*OUTPUT[a,b,j][i] for j,w in enumerate(weights)),I())/sum(weights,I())
                assert actual.lo<=direct.hi and direct.lo<=actual.hi
                assert actual.hi-actual.lo<D('1e-30');checks+=1
    assert len(switch_keys())==55
    return checks


def main():
    assert shutil.disk_usage('.').free>=20*(1<<30),'disk guard'
    start=time.monotonic();getcontext().prec=80 if '--precision80' in sys.argv else 40
    checks=exact_tests(start+40);print(dict(exact_and_interval_checks=checks,precision=getcontext().prec),flush=True)
    for radius in (0,4,8,12,16):
        assert time.monotonic()-start<40,'switch-ball deadline'
        bounds=switch_squared_bounds(radius);key=max(bounds,key=lambda k:bounds[k].hi)
        assert bounds[key].hi<32
        print(dict(radius=radius,squared_gradient_upper=str(bounds[key].hi),worst_switch=key,
                   all_centers=True,all_switches=len(bounds),whole_code_certificate=False),flush=True)
    print(dict(seconds=time.monotonic()-start,far_tail_squared_bound=32),flush=True)
    assert time.monotonic()-start<45


if __name__=='__main__':main()
