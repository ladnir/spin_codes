"""Complete integer composition lists near the least-active mixture vertex.

None means the bounded enumeration was inconclusive. A returned tuple is
the entire list, including label-type distinctions. No partial list may
be used as a proof of coverage.
"""
from fractions import Fraction as Q
from math import lcm
from integer_projection import exceptional_count_cap


def boundary_compositions(features,active,cell,slots,q_min,*,term_limit=16,node_limit=20000,list_limit=32):
    if (type(slots) is not int or not 1<=q_min<=slots or len(features)!=len(active)
            or len(cell)!=6 or min(term_limit,node_limit,list_limit)<1
            or any(len(f)!=3 or any(x<0 for x in f) for f in features)
            or any(b not in (0,1) for b in active)
            or any(not 0<=lo<=hi for lo,hi in zip(cell[::2],cell[1::2]))):
        raise ValueError('nonnegative features, box and bounded integer occupancy required')
    features=[tuple(map(Q,f)) for f in features]
    inactive=[i for i,b in enumerate(active) if not b]
    if len(inactive)!=1 or any(features[inactive[0]]):return None
    selected=[i for i,b in enumerate(active) if b]
    if not selected:return ()
    minimum=min(sum(features[i][:2]) for i in selected)
    anchors=[i for i in selected if sum(features[i][:2])==minimum]
    if len(anchors)!=1 or minimum<=0:return None
    anchor=anchors[0];budget=slots*(cell[1]+cell[3])-q_min*minimum
    if budget<0:return ()
    exceptions=sorted((sum(features[i][:2])-minimum,i) for i in selected if i!=anchor
                      and all(x<=slots*hi for x,hi in zip(features[i],cell[1::2])))
    exceptions=[(cost,i) for cost,i in exceptions if cost<=budget]
    if exceptional_count_cap([(cost,features[i]) for cost,i in exceptions],budget,slots*cell[5],slots)>term_limit:return None
    denominator=lcm(*(x.denominator for f in features for x in f),minimum.denominator)
    vectors=[tuple(int(x*denominator) for x in f) for f in features]
    kinds=[(int(cost*denominator),i) for cost,i in exceptions]
    remaining=int(budget*denominator);lower=[];upper=[]
    for lo,hi in zip(cell[::2],cell[1::2]):
        x,y=slots*lo*denominator,slots*hi*denominator
        lower.append(-(-x.numerator//x.denominator));upper.append(y.numerator//y.denominator)
    if any(lo>hi for lo,hi in zip(lower,upper)):return ()
    nodes=0;result=[];counts=[0]*len(features)
    class Exhausted(Exception):pass
    def visit(index,total,count,left):
        nonlocal nodes
        nodes+=1
        if nodes>node_limit:raise Exhausted
        if any(x>hi for x,hi in zip(total,upper)):return
        if index==len(kinds):
            first,last=max(0,q_min-count),slots-count
            for x,a,lo,hi in zip(total,vectors[anchor],lower,upper):
                if a:
                    first=max(first,-(-(lo-x)//a));last=min(last,(hi-x)//a)
                elif not lo<=x<=hi:return
            if len(result)+max(0,last-first+1)>list_limit:raise Exhausted
            for n in range(first,last+1):
                counts[anchor]=n;counts[inactive[0]]=slots-count-n
                result.append(tuple(counts))
            return
        cost,i=kinds[index]
        for n in range(min(slots-count,left//cost)+1):
            counts[i]=n
            visit(index+1,tuple(x+n*y for x,y in zip(total,vectors[i])),count+n,left-n*cost)
        counts[i]=0
    try:visit(0,(0,0,0),0,remaining)
    except Exhausted:return None
    return tuple(sorted(result))
