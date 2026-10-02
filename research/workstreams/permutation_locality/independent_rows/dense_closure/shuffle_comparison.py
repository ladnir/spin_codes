"""Shuffle comparisons that retain randomness in the packet-type counts.

The exact homogeneous bound is special-purpose. The posterior bound accepts
arbitrary mixtures: guaranteed constant packets give a lower bound on the
number of equal observations, whose posterior type counts remain random.
"""
from fractions import Fraction as Q
from math import comb,factorial,lgamma,log,log1p


def structure(components,counts):
    law=None;q=0;fixed={};slots=sum(counts)
    if len(components)!=len(counts) or any(type(c) is not int or c<0 for c in counts):
        raise ValueError('nonnegative integer component counts required')
    for row,n in zip(components,counts):
        if not n:continue
        p=row[2]
        if p[0]==1:continue
        constant=next((j for j,x in enumerate(p) if x==1),None)
        if constant is not None:fixed[constant]=fixed.get(constant,0)+n
        elif law is None:law=p;q=n
        elif p==law:q+=n
        else:raise ValueError('more than one stochastic packet law')
    if law is not None and any(law[j] for j in fixed):
        raise ValueError('constant categories must be disjoint from stochastic support')
    return slots,q,law,tuple(sorted(fixed.items()))


def binomial_loss(slots,q,p):
    """Exact maximum of Bin(q,p)/Bin(slots,q*p/slots) on their counts."""
    p=Q(p)
    if type(slots) is not int or type(q) is not int or not 0<=q<=slots or not 0<=p<=1:
        raise ValueError('valid binomial geometry required')
    if not q or q==slots or not p:return Q(1)
    mean=q*p
    k=min(q,mean.numerator//mean.denominator+1)
    reference=Q(q,slots)*p
    return (Q(comb(q,k),comb(slots,k))*p**k*(1-p)**(q-k)
            / (reference**k*(1-reference)**(slots-k)))


def binomial_logloss(slots,q,p):
    if not q or q==slots or p<=0:return 0.
    k=min(q,int(q*p)+1);reference=q*p/slots
    value=lgamma(q+1)-lgamma(q-k+1)-lgamma(slots+1)+lgamma(slots-k+1)
    if k:value+=k*(log(p)-log(reference))
    if k<q:value+=(q-k)*log1p(-p)
    return value-(slots-k)*log1p(-reference)


def fixed_loss(slots,fixed):
    remaining=slots-sum(n for _,n in fixed)
    multiplicity=factorial(slots)//factorial(remaining)
    mass=Q(remaining,slots)**remaining
    for _,n in fixed:
        multiplicity//=factorial(n);mass*=Q(n,slots)**n
    return 1/(multiplicity*mass)


def exact(geometry,tilt):
    slots,q,law,fixed=geometry
    remaining=slots-sum(n for _,n in fixed)
    p=Q(0) if law is None else sum(p*t for p,t in zip(law[1:],tilt[1:]))/sum(p*t for p,t in zip(law,tilt))
    return fixed_loss(slots,fixed)*binomial_loss(remaining,q,p)


def floating(geometry,tilt):
    slots,q,law,fixed=geometry
    remaining=slots-sum(n for _,n in fixed)
    value=lgamma(remaining+1)-lgamma(slots+1)
    if remaining:value-=remaining*log(remaining/slots)
    for _,n in fixed:value+=lgamma(n+1)-n*log(n/slots)
    p=0. if law is None else sum(float(p)*t for p,t in zip(law[1:],tilt[1:]))/sum(float(p)*t for p,t in zip(law,tilt))
    return value+binomial_logloss(remaining,q,p)


def multinomial_mass(counts,probabilities):
    coefficient=factorial(sum(counts));mass=Q(1)
    for n,p in zip(counts,probabilities):
        coefficient//=factorial(n);mass*=p**n
    return coefficient*mass


def multinomial_mode(total,probabilities,caps):
    """Exact mode subject to coordinate caps, by decreasing marginal gains.

    Starting at min(cap_i,floor(total*p_i)) takes only gains >=1/total;
    all remaining gains are <=1/total. Greedy completion is therefore
    exact. Zero-probability coordinates cannot carry positive mass.
    """
    probabilities=tuple(map(Q,probabilities))
    if (type(total) is not int or total<0 or len(caps)!=len(probabilities)
            or any(type(c) is not int or c<0 for c in caps)
            or not probabilities or min(probabilities)<0 or sum(probabilities)!=1):
        raise ValueError('normalized probabilities and nonnegative integer caps required')
    if total>sum(c for c,p in zip(caps,probabilities) if p):return None
    counts=[min(c,int(total*p)) for c,p in zip(caps,probabilities)]
    for _ in range(total-sum(counts)):
        i=max((i for i,p in enumerate(probabilities) if p and counts[i]<caps[i]),
              key=lambda i:probabilities[i]/(counts[i]+1))
        counts[i]+=1
    return tuple(counts)


def posterior_structure(components,counts):
    if (len(components)!=len(counts) or any(type(c) is not int or c<0 for c in counts)
            or not sum(counts)):
        raise ValueError('positive total of nonnegative integer counts required')
    merged={}
    for row,n in zip(components,counts):
        if not n:continue
        law=tuple(map(Q,row[2]))
        if len(law)!=5 or min(law)<0 or sum(law)!=1:raise ValueError('five-category probability laws required')
        merged[law]=merged.get(law,0)+n
    laws=tuple(merged);ns=tuple(merged.values());slots=sum(ns)
    fixed=tuple(sum(n for n,law in zip(ns,laws) if law[j]==1) for j in range(5))
    prior=multinomial_mass(ns,[Q(n,slots) for n in ns])
    logprior=lgamma(slots+1)-sum(lgamma(n+1) for n in ns)+sum(n*log(n/slots) for n in ns)
    return ns,laws,fixed,prior,logprior


def posterior_exact(geometry,tilt):
    """Pointwise shuffle/iid bound; exact rationals, not an exact supremum."""
    ns,laws,fixed,prior,_=geometry
    tilt=tuple(map(Q,tilt))
    if len(tilt)!=5 or min(tilt)<=0:raise ValueError('five positive tilts required')
    tilted=[tuple(p*t for p,t in zip(law,tilt)) for law in laws]
    tilted=[tuple(p/sum(law) for p in law) for law in tilted]
    peak=Q(1)
    for j,m in enumerate(fixed):
        if not m:continue
        posterior=[n*law[j] for n,law in zip(ns,tilted)]
        total=sum(posterior);posterior=[p/total for p in posterior]
        mode=multinomial_mode(m,posterior,ns)
        mass=Q(0) if mode is None else multinomial_mass(mode,posterior)
        peak=min(peak,mass)
    return peak/prior


def posterior_floating(geometry,tilt):
    ns,laws,fixed,_,logprior=geometry
    tilted=[[float(p)*t for p,t in zip(law,tilt)] for law in laws]
    tilted=[[p/sum(law) for p in law] for law in tilted]
    logpeak=0.
    for j,m in enumerate(fixed):
        if not m:continue
        posterior=[n*law[j] for n,law in zip(ns,tilted)]
        total=sum(posterior);posterior=[p/total for p in posterior]
        counts=[min(n,int(m*p)) for n,p in zip(ns,posterior)]
        # Floating arithmetic proposes a mode only. Replay recomputes it
        # with rationals and does not trust these counts or this score.
        if sum(counts)>m:counts=[0]*len(ns)
        for _ in range(m-sum(counts)):
            i=max((i for i,p in enumerate(posterior) if p and counts[i]<ns[i]),
                  key=lambda i:posterior[i]/(counts[i]+1))
            counts[i]+=1
        value=lgamma(m+1)-sum(lgamma(c+1) for c in counts)
        value+=sum(c*log(p) for c,p in zip(counts,posterior) if c)
        logpeak=min(logpeak,value)
    return logpeak-logprior
