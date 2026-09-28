"""Exact first moments of input/expansion overlap for fixed packet shapes.

This independent candidate does not enter the production verifier or the
cached local-family builder. See ../../OVERLAP_MEAN.md for the argument.
All polynomial arithmetic and overlap allocations use Python integers and
fractions. Only the final exponential moment uses outward Arb arithmetic.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement,product
from math import comb,factorial,prod
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from mass_density_screen import baseline
from occupancy_memory import C,Z
from flint import arb


def kraw(n,h,b):
    return sum((-1)**r*comb(h,r)*comb(n-h,b-r)
               for r in range(max(0,b-(n-h)),min(b,h)+1))


def shape_tuples(maximum,width):
    return sorted((s for s in product(range(maximum+1),repeat=width) if sum(s)<=maximum),
                  key=lambda s:(sum(s),s))


def polynomial(histogram,width,maximum):
    """Coefficients of the character polynomial on the other windows."""
    shapes=shape_tuples(maximum,width)
    index={s:i for i,s in enumerate(shapes)}
    previous=[[(index[tuple(c-int(b==i) for i,c in enumerate(s))],b)
               for b in range(width) if s[b]] for s in shapes]
    coefficients=[0]*len(shapes); coefficients[0]=1
    processed=0
    for h,copies in enumerate(histogram):
        factors=[kraw(width,h,b) for b in range(1,width+1)]
        for _ in range(copies):
            processed+=1
            for i in range(len(shapes)-1,0,-1):
                if sum(shapes[i])>processed: continue
                coefficients[i]+=sum(coefficients[p]*factors[b] for p,b in previous[i])
    return dict(zip(shapes,coefficients))


def moments(expansion,feedback,maximum,*,width=4):
    """shape -> (denominator, exact mean overlap, universal translated upper).

expansion lists the output word of each state-basis vector. feedback lists
the state image of each input-basis vector. A shape samples distinct windows
and uniform masks of its prescribed positive weights, including every
assignment of those weights to the windows.
"""
    if (type(width) is not int or width<1 or type(maximum) is not int
            or not feedback or len(feedback)%width or not 1<=maximum<=len(feedback)//width
            or not expansion or any(type(x) is not int or x<0 for x in (*expansion,*feedback))
            or any(x>=1<<len(feedback) for x in expansion)
            or any(x>=1<<len(expansion) for x in feedback)):
        raise ValueError('matching binary linear maps and complete packet windows required')
    windows=len(feedback)//width
    shapes=shape_tuples(maximum,width)[1:]
    signed=[0]*len(shapes); absolute=[0]*len(shapes)
    denominators=[comb(windows,sum(s))*factorial(sum(s))//prod(factorial(c) for c in s)
                  *prod(comb(width,b)**s[b-1] for b in range(1,width+1)) for s in shapes]
    weights=[sum((b+1)*c for b,c in enumerate(s)) for s in shapes]
    polynomials={}
    for i in range(len(feedback)):
        character=sum(((image>>i)&1)<<b for b,image in enumerate(expansion))
        signs=[(character&column).bit_count()%2 for column in feedback]
        window_weights=[sum(signs[width*w:width*(w+1)]) for w in range(windows)]
        own=window_weights[i//width]
        hist=Counter(window_weights); hist[own]-=1
        key=tuple(hist[h] for h in range(width+1))
        if key not in polynomials:
            polynomials[key]=polynomial(key,width,maximum-1)
        coefficients=polynomials[key]
        derivative=[(-1)**signs[i]*kraw(width-1,own-signs[i],b-1) for b in range(1,width+1)]
        for j,s in enumerate(shapes):
            value=sum(derivative[b]*coefficients[tuple(c-int(b==a) for a,c in enumerate(s))]
                      for b in range(width) if s[b])
            # Each coordinate is 1 in exactly W*D/output_length inputs.
            assert abs(value)*len(feedback)<=weights[j]*denominators[j]
            signed[j]+=value; absolute[j]+=abs(value)
    result={}
    for s,W,D,S,A in zip(shapes,weights,denominators,signed,absolute):
        shape=tuple(b+1 for b,c in enumerate(s) for _ in range(c))
        mean=Q(W*D-S,2*D); upper=Q(W*D+A,2*D)
        assert 0<=mean<=upper<=W
        result[shape]=D,mean,upper
    return result


def build(maximum=10):
    from group_moment import maps
    images,feedback,_=maps()
    expansion=[images[1<<b] for b in range(19)]
    result=moments(expansion,feedback,maximum)
    for j in range(1,maximum+1):
        D,mean,upper=result[(4,)*j]
        print('EXACT OVERLAP all-four',j,'mean',mean,'translated upper',upper,flush=True)
    print('EXACT OVERLAP complete shapes',len(result),flush=True)
    return result


def chord_histogram(classes,weight,mean):
    """Positive output-weight envelope from class masses and an overlap cap.

Greedily allocate the overlap budget to the smallest expansion weights.
The allocation is independent of z in (0,1]. It maximizes the chord upper,
not the actual output moment. An overlarge budget can simply remain unused.
"""
    if (type(weight) is not int or weight<1 or not isinstance(mean,Q)
            or not 0<=mean<=weight
            or any(type(v) is not int or v<=0 or not isinstance(p,Q) or p<0 for v,p in classes.items())
            or sum(classes.values(),Q(0))>1):
        raise ValueError('positive expansion classes, exact masses, and overlap cap required')
    remaining=mean; histogram={}
    for v,p in sorted(classes.items()):
        overlap=min(remaining,weight*p); remaining-=overlap
        high=overlap/weight
        for w,mass in ((v-weight,high),(v+weight,p-high)):
            histogram[w]=histogram.get(w,Q(0))+mass
    assert all(p>=0 for p in histogram.values())
    assert sum(histogram.values(),Q(0))==sum(classes.values(),Q(0))
    return {w:p for w,p in histogram.items() if p}


def zero_refine(base,overlap,feedback,tilt,penalty,*,minimum=5,maximum=10,rounds=2):
    """Tighten only the scalar mature-density to zero coefficient.

Feedback zero and expansion-class counts must be exact. The feedback peak
need not be exact here. Apply before coupled mass-column replacement.
"""
    if not 1<=minimum<=maximum<len(base) or not 0<Q(penalty)<=1 or Q(tilt)<=0 or rounds<1:
        raise ValueError('valid complete shape range and positive witnesses required')
    required={s for j in range(minimum,maximum+1) for s in combinations_with_replacement(range(1,5),j)}
    if not required<=set(overlap) or not required<=set(feedback):
        raise ValueError('complete overlap and feedback censuses required')
    bounds={j:arb(0) for j in range(minimum,maximum+1)}
    for shape in required:
        zero,_,D,counts=feedback[shape]
        if D!=overlap[shape][0] or zero+sum(counts.values())!=D:
            raise ValueError('feedback and overlap counting measures disagree')
        hist=chord_histogram({v:Q(int(n),int(D)) for v,n in counts.items()},sum(shape),overlap[shape][1])
        value=sum((arb(p.numerator)/p.denominator*(-arb(tilt)*w).exp() for w,p in hist.items()),arb(0))
        scale=Q(penalty)**shape.count(4)/2**rounds
        value=baseline.up(value*scale.numerator/scale.denominator)
        bounds[len(shape)]=max(bounds[len(shape)],value)
    result=[t*1 for t in base]; changes=0
    for j,value in bounds.items():
        changes+=int(value<result[j][C,Z])
        result[j][C,Z]=min(value,result[j][C,Z])
    print('OVERLAP CHORD tightened zero-return entries',changes,flush=True)
    return result


if __name__=='__main__':
    build()
