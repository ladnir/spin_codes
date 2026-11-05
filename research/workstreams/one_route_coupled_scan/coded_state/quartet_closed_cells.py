"""Uniform source-population cells with exact original-inactive labels.

Positive coefficient formulas; FLOAT proposals only. No distance cert.
"""
import math
import numpy as np
from quartet_closed_population import multinomial


def vertices(cell,width=512):
    alo,ahi,slo,shi=cell
    pairs={(a,s) for a in (alo,ahi) for s in (slo,shi) if 0<=s<=a<=width}
    pairs.update((a,a) for a in (alo,ahi,slo,shi) if alo<=a<=ahi and slo<=a<=shi)
    assert pairs
    return np.array([(width-a,a-s,s) for a,s in sorted(pairs)],dtype=int)


def bounds(caps,cell,width=512):
    """All inactive counts; Cauchy in source counts, exact univariate extraction."""
    ns=vertices(cell,width);groups=width//4;active=np.flatnonzero(ns.max(axis=0)>0)
    omitted=int(active[np.argmax(ns.mean(axis=0)[active])]);kept=[int(i) for i in active if i!=omitted]
    base=[];powers=[];js=[]
    for m,values in caps.items():
        if any(m[i] for i in range(3) if i not in active):continue
        for j,value in enumerate(values):
            base.append(value+math.log2(multinomial(m)));powers.append([m[i] for i in kept]+[j]);js.append(j)
    base=np.array(base);powers=np.array(powers,dtype=float);js=np.array(js)
    target=np.tile(np.r_[ns.mean(axis=0)[kept]/groups,0.],(width+1,1));target[:,-1]=np.arange(width+1)/groups
    x=np.zeros_like(target);dimension=x.shape[1]
    def evaluate(x):
        terms=base[None,:]+x@powers.T;peak=terms.max(axis=1);raw=np.exp2(terms-peak[:,None]);p=raw/raw.sum(axis=1)[:,None]
        total=peak+np.log2(raw.sum(axis=1));mean=p@powers
        covariance=np.einsum('bt,ti,tj->bij',p,powers,powers)-np.einsum('bi,bj->bij',mean,mean)
        return total,p,mean,covariance
    for _ in range(32):
        total,p,mean,cov=evaluate(x)
        direction=np.linalg.solve(math.log(2)*cov+1e-10*np.eye(dimension)[None,:,:],(target-mean)[:,:,None])[:,:,0]
        direction/=np.maximum(1.,np.max(abs(direction),axis=1)/3)[:,None]
        old=total-np.sum(target*x,axis=1);rate=np.ones(width+1)
        for _ in range(20):
            proposal=np.clip(x+rate[:,None]*direction,-40,40)
            accept=evaluate(proposal)[0]-np.sum(target*proposal,axis=1)<=old+1e-12
            if np.all(accept):break
            rate=np.where(accept,rate,rate*.5)
        x=proposal
    total,p,mean,cov=evaluate(x)
    q=np.column_stack([p[:,js==j].sum(axis=1) for j in range(5)])
    current=np.zeros((width+1,width+1));current[:,0]=1.
    for _ in range(groups):
        after=np.zeros_like(current)
        for j in range(5):after[:,j:]+=q[:,j,None]*current[:,:width+1-j]
        current=after
    mass=np.diag(current);assert np.all(mass>0)
    lognormal=np.array([math.log2(multinomial(n)) for n in ns])
    source=-x[:,:len(kept)]@ns[:,kept].T-lognormal[None,:]
    return groups*total-np.arange(width+1)*x[:,-1]+source.max(axis=1)+np.log2(mass)


def tests():
    from quartet_closed_count import bound
    from quartet_two_cohort import compositions
    caps={m:np.log2([1+3*m[1]+7*m[2]+2*m[1]*m[2]+j*(m[0]+1) for j in range(5)]) for m in compositions(4,3)}
    checked=0
    for cell in ((0,4,0,4),(4,8,0,4),(4,8,4,8),(8,8,0,8),(0,0,0,0)):
        upper=bounds(caps,cell,8)
        for a in range(cell[0],cell[1]+1):
            for s in range(cell[2],min(a,cell[3])+1):
                n=(8-a,a-s,s)
                # Exact two-quartet coefficient, not a comparison to a
                # separately optimized upper bound.
                for inactive in range(9):
                    total=0.
                    for m,values in caps.items():
                        other=tuple(ni-mi for ni,mi in zip(n,m))
                        if min(other)<0:continue
                        for j in range(5):
                            if 0<=inactive-j<=4:
                                total+=multinomial(m)*multinomial(other)*2**(values[j]+caps[other][inactive-j])
                    exact=total/multinomial(n)
                    assert math.log2(exact)<=upper[inactive]+1e-8;checked+=1
    print(dict(uniform_population_cell_comparisons=checked),flush=True)
