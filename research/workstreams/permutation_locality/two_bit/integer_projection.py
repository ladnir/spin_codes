"""Conservative exact reachability pruning near a mixture-coordinate boundary.

Return 'possibly reachable' when the bounded enumeration would be costly.
Only exhaustive rational infeasibility is evidence for discarding a cell.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import gcd,lcm


@lru_cache(maxsize=8192)
def possibly_reachable(values,lo,hi,slots,term_limit=12,node_limit=20000):
    values=tuple(sorted(set(map(Q,values))))
    lo,hi=Q(lo),Q(hi)
    if (not 0<=lo<=hi or type(slots) is not int or slots<0 or min(term_limit,node_limit)<1
            or any(v<=0 for v in values)):
        raise ValueError('positive increments, bounded interval, and nonnegative slot count required')
    if lo==0:return True
    if not values or not slots:return False
    # Scale once; the search itself uses small integers, not Fraction objects.
    denominator=lcm(*(v.denominator for v in values))
    increments=tuple(v.numerator*(denominator//v.denominator) for v in values)
    divisor=gcd(*increments);increments=tuple(v//divisor for v in increments)
    low_scaled=lo*denominator/divisor;high_scaled=hi*denominator/divisor
    lower=-(-low_scaled.numerator//low_scaled.denominator)
    upper=high_scaled.numerator//high_scaled.denominator
    if lower>upper:return False
    terms=min(slots,upper//increments[0])
    if terms>term_limit:return True
    nodes=0
    def search(index,total,left):
        nonlocal nodes
        nodes+=1
        if nodes>node_limit:return True  # Inconclusive, never an emptiness proof.
        if lower<=total<=upper:return True
        if index==len(increments) or not left or total+left*increments[-1]<lower:return False
        value=increments[index]
        for count in range(min(left,(upper-total)//value)+1):
            if search(index+1,total+count*value,left-count):return True
        return False
    return search(0,0,terms)


def integer_projection_empty(features,cell,slots):
    """Every actual composition has at most slots nonnegative increments."""
    if (type(slots) is not int or slots<1 or not cell or len(cell)%2
            or any(not 0<=lo<=hi for lo,hi in zip(cell[::2],cell[1::2]))
            or any(len(f)*2!=len(cell) or any(x<0 for x in f) for f in features)):
        raise ValueError('nonnegative features, box and positive slot count required')
    # A component that alone exceeds any upper coordinate cannot occur.
    eligible=[f for f in features if all(x<=slots*hi for x,hi in zip(f,cell[1::2]))]
    for j,(lo,hi) in enumerate(zip(cell[::2],cell[1::2])):
        if lo==0:continue
        increments=tuple(sorted({f[j] for f in eligible if f[j]>0}))
        if not possibly_reachable(increments,slots*lo,slots*hi,slots):return True
    return False


def exceptional_count_cap(exceptions,budget,central_cap,slots):
    if not exceptions:return 0
    bound=min(slots,int(budget/min(cost for cost,_ in exceptions)))
    central=[f[2] for _,f in exceptions if f[2]>0]
    noncentral=[cost for cost,f in exceptions if not f[2]]
    # Separate upper bounds remain valid even though both classes consume
    # the same excess budget. This avoids charging unlimited central types
    # in a cell that permits only one central row.
    central_count=min(slots,int(central_cap/min(central))) if central else 0
    other_count=min(slots,int(budget/min(noncentral))) if noncentral else 0
    return min(bound,central_count+other_count)


def boundary_composition_empty(features,active,cell,slots,q_min,term_limit=16,node_limit=20000):
    """Enumerate exceptional types near the least-active mixture vertex.

    If a is the least active feature sum f0+f1, every active type pays
    nonnegative excess f0+f1-a. A small upper excess allows enumeration
    of all nonminimal types. The remaining minimal-type count is checked
    by exact integer intervals in every coordinate simultaneously.
    """
    if (len(features)!=len(active) or not 1<=q_min<=slots or len(cell)!=6
            or any(b not in (0,1) for b in active)
            or any(len(f)!=3 or any(x<0 for x in f) for f in features)
            or any(not 0<=lo<=hi for lo,hi in zip(cell[::2],cell[1::2]))):
        raise ValueError('nonnegative three-coordinate features and bounded occupancy required')
    if any(any(f) for f,b in zip(features,active) if not b):return False
    positive=[tuple(map(Q,f)) for f,b in zip(features,active) if b]
    if not positive:return True
    minimum=min(f[0]+f[1] for f in positive)
    anchors={f for f in positive if f[0]+f[1]==minimum}
    if len(anchors)!=1 or minimum<=0:return False
    anchor=next(iter(anchors));budget=slots*(cell[1]+cell[3])-q_min*minimum
    if budget<0:return True
    exceptional=sorted(set((f[0]+f[1]-minimum,f) for f in positive
                           if f!=anchor and all(x<=slots*hi for x,hi in zip(f,cell[1::2]))))
    exceptional=[(cost,f) for cost,f in exceptional if cost<=budget]
    max_terms=exceptional_count_cap(exceptional,budget,slots*cell[5],slots)
    if max_terms>term_limit:return False
    denominator=lcm(*(x.denominator for f in [anchor,*(f for _,f in exceptional)] for x in f),
                    minimum.denominator)
    anchor_int=tuple(int(x*denominator) for x in anchor)
    kinds=[(int(cost*denominator),tuple(int(x*denominator) for x in f)) for cost,f in exceptional]
    scaled=budget*denominator;remaining_budget=scaled.numerator//scaled.denominator
    lower=[];upper=[]
    for lo,hi in zip(cell[::2],cell[1::2]):
        x,y=slots*lo*denominator,slots*hi*denominator
        lower.append(-(-x.numerator//x.denominator));upper.append(y.numerator//y.denominator)
    if any(lo>hi for lo,hi in zip(lower,upper)):return True
    nodes=0
    def reachable(index,total,count,remaining):
        nonlocal nodes
        nodes+=1
        if nodes>node_limit:return True
        if any(x>hi for x,hi in zip(total,upper)):return False
        first,last=max(0,q_min-count),slots-count
        for x,a,lo,hi in zip(total,anchor_int,lower,upper):
            if a:
                first=max(first,-(-(lo-x)//a));last=min(last,(hi-x)//a)
            elif not lo<=x<=hi:first,last=1,0;break
        if first<=last:return True
        if index==len(kinds):return False
        cost,feature=kinds[index]
        for number in range(min(slots-count,remaining//cost)+1):
            if reachable(index+1,tuple(x+number*y for x,y in zip(total,feature)),
                         count+number,remaining-number*cost):return True
        return False
    return not reachable(0,(0,0,0),0,remaining_budget)
