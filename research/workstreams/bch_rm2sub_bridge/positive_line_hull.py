"""Exact pruning of redundant nonnegative Bernoulli coefficient pairs.

Pairs (a,b) represent a*A+b*B for A,B>=0. Their upper line envelope at
x=B/A is determined by rational breakpoints. No floating hull decision is
trusted. Endpoint and kept-line intersection checks independently verify
that every discarded line stays below the retained envelope on x>=0.
"""
from fractions import Fraction as F


def indices(left,right):
    pairs=[(F.from_float(float(a)),F.from_float(float(b))) for a,b in zip(left,right)]
    lines=sorted((b,a,i) for i,(a,b) in enumerate(pairs));hull=[];starts=[]
    for slope,intercept,index in lines:
        if hull and slope==hull[-1][0]:
            if intercept<=hull[-1][1]:continue
            hull.pop();starts.pop()
        start=None
        while hull:
            start=(hull[-1][1]-intercept)/(slope-hull[-1][0])
            if starts[-1] is None or start>starts[-1]:break
            hull.pop();starts.pop()
        if not hull:start=None
        hull.append((slope,intercept,index));starts.append(start)
    keep=[line[2] for i,line in enumerate(hull) if i+1==len(hull) or starts[i+1]>=0]
    assert keep
    retained=[pairs[i] for i in keep]
    assert max(b for a,b in retained)>=max(b for a,b in pairs)
    points={F(0)}
    for a,b in retained:
        for c,d in retained:
            if b!=d:
                x=(c-a)/(b-d)
                if x>=0:points.add(x)
    for x in points:
        upper=max(a+b*x for a,b in retained)
        assert all(a+b*x<=upper for a,b in pairs)
    return keep
