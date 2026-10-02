"""Complete integer composition enumeration near the least-dense vertex.

None means the work limit was reached, never a partial coverage claim.
"""
from fractions import Fraction as Q
from math import lcm


def compositions(features,active,cell,slots,q_min,*,term_limit=8,node_limit=12000,list_limit=32):
    dimension=len(cell)//2
    if (not dimension or len(cell)!=2*dimension or not 1<=q_min<=slots
            or len(features)!=len(active) or any(len(f)!=dimension or min(f)<0 for f in features)
            or any(b not in (0,1) for b in active)
            or any(not 0<=lo<=hi for lo,hi in zip(cell[::2],cell[1::2]))):
        raise ValueError('bounded nonnegative feature box required')
    features=[tuple(map(Q,f)) for f in features]
    inactive=[i for i,a in enumerate(active) if not a]
    selected=[i for i,a in enumerate(active) if a]
    if len(inactive)!=1 or any(features[inactive[0]]) or not selected:return None
    minimum=min(sum(features[i]) for i in selected)
    anchors=[i for i in selected if sum(features[i])==minimum]
    if minimum<=0 or len(anchors)!=1:return None
    anchor=anchors[0];budget=slots*sum(cell[1::2])-q_min*minimum
    if budget<0:return ()
    kinds=sorted((sum(features[i])-minimum,i) for i in selected if i!=anchor
                 and all(x<=slots*hi for x,hi in zip(features[i],cell[1::2])))
    kinds=[(cost,i) for cost,i in kinds if cost<=budget]
    if kinds and budget//kinds[0][0]>term_limit:return None
    denominator=lcm(*(x.denominator for f in features for x in f))
    vectors=[tuple(int(x*denominator) for x in f) for f in features]
    kinds=[(int(cost*denominator),i) for cost,i in kinds]
    lower=[];upper=[]
    for lo,hi in zip(cell[::2],cell[1::2]):
        x,y=slots*lo*denominator,slots*hi*denominator
        lower.append(-(-x.numerator//x.denominator));upper.append(y.numerator//y.denominator)
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
                counts[anchor]=n;counts[inactive[0]]=slots-count-n;result.append(tuple(counts))
            return
        cost,i=kinds[index]
        for n in range(min(slots-count,left//cost)+1):
            counts[i]=n
            visit(index+1,tuple(x+n*y for x,y in zip(total,vectors[i])),count+n,left-n*cost)
        counts[i]=0
    try:visit(0,(0,)*dimension,0,int(budget*denominator))
    except Exhausted:return None
    return tuple(sorted(result))
