"""Binary64 shape-resolved nine-coordinate pilot for categorical conditioning.

No outward certificate. This deliberately omits the two mature-tail
coordinates. It retains the existing mass/density/fresh invariant and
averages actual input shapes, not entrywise shape maxima.
"""
from itertools import product
from math import exp,prod,factorial,comb
import numpy as np
from scipy.special import gammaln,xlogy
from occupancy_memory import Z,F,M,C,U,epoch_operators
from occupancy_multi_average import multiplicities
from fresh_collision import probability as fresh_collision_probability


def asarray(t):
    return np.array([[float(t[i,j]) for j in range(t.ncols())] for i in range(t.nrows())])


def raw_coarse(weights,spectrum,tilt,pair_bound):
    j=len(weights);weight=sum(weights);levels=sorted(spectrum);m=(1<<19)-1
    atom=min(1/((33-j)*max(comb(4,w) for w in weights)),pair_bound*32*31/((34-j)*(33-j)))
    f=exp(-tilt*max(0,48-weight));t=np.zeros((9,9))
    t[Z,Z]=exp(-tilt*weight)*atom;t[Z,M]=exp(-tilt*weight);t[Z,C]=t[Z,Z]
    fresh=f*min(atom,1/32)
    t[F,Z]=fresh/2+f/(2*m);t[F,M]=f/2;t[F,C]=fresh/2
    t[M,Z]=f/(2*m);t[M,M]=f/2;t[C,Z]=f/2;t[C,C]=f/2
    for i in (F,M):
        for k,v in enumerate(levels):t[i,U+k]=f*spectrum[v]/(2*m)
    for i,v in enumerate(levels):
        peak=exp(-tilt*max(0,v-weight))
        t[U+i,Z]=peak/(2*spectrum[v])+peak/(2*m)
        t[U+i,M]=peak/2;t[U+i,C]=peak/(2*spectrum[v])
        for k,w in enumerate(levels):t[U+i,U+k]=peak*spectrum[w]/(2*m)
    return t


def refine_shape(t,weights,spectrum,fresh,window,zeros,prepared,full,tilt):
    """Float counterparts of existing per-shape bounds, then two updates."""
    t=t.copy();j=len(weights);weight=sum(weights);levels=sorted(spectrum);m=(1<<19)-1
    _,distributions,fresh_moments=fresh
    def average(hist):return sum(c*exp(-tilt*w) for w,c in hist.items())/sum(hist.values())
    def put(i,k,value):t[i,k]=min(t[i,k],value)
    if j==0:
        moment=max(average(hist) for hist in distributions.values())
    elif j==1:
        b=weights[0]
        moment=max(average(hist) for (a,c),hist in fresh_moments.items() if c==b)
    else:
        moment=max(sum(c*exp(-tilt*max(0,v-weight)) for v,c in hist.items())/sum(hist.values())
                   for hist in distributions.values())
    if j:put(F,M,moment/2)
    for k,v in enumerate(levels):put(F,U+k,moment*spectrum[v]/(2*m))
    if j==1:
        mass,density,zero,den=window[weights[0]]
        mass/=den;density/=den;zero/=den
        put(M,M,mass/2);put(M,Z,mass/(2*m));put(C,Z,zero/2);put(C,C,density/2)
        for k,v in enumerate(levels):
            put(M,U+k,mass*spectrum[v]/(2*m));put(U+k,C,density/(2*spectrum[v]))
    if j>=2:
        moments={v:min(1.,exp(-tilt*v)*prod(((128-v)*exp(-tilt*w)+v*exp(tilt*w))/128 for w in weights))
                 for v in levels}
        arbitrary=max(moments.values())
        fresh_mean=max(sum(c*moments[v] for v,c in hist.items())/sum(hist.values()) for hist in distributions.values())
        put(F,M,fresh_mean/2);put(M,M,arbitrary/2);put(M,Z,arbitrary/(2*m))
        for k,v in enumerate(levels):
            put(F,U+k,fresh_mean*spectrum[v]/(2*m));put(M,U+k,arbitrary*spectrum[v]/(2*m))
            put(U+k,M,moments[v]/2)
            put(U+k,Z,exp(-tilt*abs(v-weight))/(2*spectrum[v])+moments[v]/(2*m))
            for l,w in enumerate(levels):put(U+k,U+l,moments[v]*spectrum[w]/(2*m))
    if j in (2,3):
        tables={1:{},2:{},**prepared[0][2]}
        cancellation=max(float(fresh_collision_probability(tables,a,weights)) for a in range(1,5))
        put(F,Z,exp(-tilt*max(0,48-weight))*(cancellation+1/m)/2)
        hist,den= zeros[j][tuple(sorted(weights))]
        put(C,Z,sum(n*exp(-tilt*w) for w,n in enumerate(hist))/den/2)
    # Same lazy/refresh decomposition as mixing_rounds.transform(r=2).
    old=t.copy();t[1:,:]*=.5;t[1:,U:U+5]=old[1:,U:U+5]*1.5
    for i in range(1,9):
        if j:
            if i==M:t[i,Z]=1.5*old[i,Z]
            elif i==F or U<=i<U+5:
                refresh=min(old[i,U+k]/spectrum[v] for k,v in enumerate(levels))
                t[i,Z]=.5*old[i,Z]+refresh
        elif U<=i<U+5:
            t[i,i]=1.5*old[i,i]-.5*exp(-tilt*levels[i-U])
    if j>=2 and tuple(sorted(weights)) in full:
        zero,peak,den,classes=full[tuple(sorted(weights))]
        p0=zero/den;nz=peak/den;all_peak=max(zero,peak)/den
        factor=exp(-tilt*weight);lower=exp(-tilt*max(0,48-weight))
        put(Z,Z,factor*p0);put(Z,M,factor*(1-p0));put(Z,C,factor*nz)
        put(F,C,.25*lower*all_peak)
        refresh=min(t[F,U+k]/spectrum[v] for k,v in enumerate(levels))
        put(F,Z,.25*lower*nz+refresh)
        put(C,Z,.25*sum(classes[v]*exp(-tilt*abs(v-weight)) for v in levels)/den)
        for i,v in enumerate(levels):
            refresh=min(t[U+i,U+k]/spectrum[w] for k,w in enumerate(levels))
            put(U+i,Z,.25*exp(-tilt*abs(v-weight))*classes[v]/(den*spectrum[v])+refresh)
            put(U+i,C,.25*exp(-tilt*abs(v-weight))*min(all_peak,1/spectrum[v]))
    assert (t>=0).all() and np.isfinite(t).all()
    return t


def build(prepared,fresh,window,zeros,full,tilt,maximum=32):
    tilt=float(tilt)
    raw=epoch_operators(prepared,str(tilt),detailed=True)
    spectrum=prepared[0][0][0]
    pair_bound=max(row[2]/row[0] for row in prepared[0][1][1].values())
    result=[]
    for j in range(maximum+1):
        shapes=[];matrices=[]
        for counts in multiplicities(j):
            weights=tuple(w for w,n in enumerate(counts,1) for _ in range(n))
            t=asarray(raw[0] if j==0 else raw[j][weights]) if j<=4 else raw_coarse(weights,spectrum,tilt,pair_bound)
            matrices.append(refine_shape(t,weights,spectrum,fresh,window,zeros,prepared,full,float(tilt)))
            shapes.append(counts)
        shapes=np.array(shapes,dtype=int)
        coefficient=gammaln(j+1)-gammaln(shapes+1).sum(axis=1)
        result.append((shapes,np.array(matrices),coefficient))
    return result


def mix(data,theta):
    theta=np.asarray(theta);assert theta.shape==(4,) and np.all(theta>=0) and abs(theta.sum()-1)<1e-12
    result=[]
    for counts,matrices,coefficient in data:
        probabilities=np.exp(coefficient+xlogy(counts,theta).sum(axis=1))
        assert abs(probabilities.sum()-1)<2e-12
        result.append(np.einsum('s,sij->ij',probabilities,matrices,optimize=False))
    return np.array(result)


def normalization_test(data):
    from fractions import Fraction as Q
    # Exact multinomial normalization through degree 32, and atomic limits.
    for counts,matrices,_ in data:
        j=int(counts[0].sum());theta=(Q(1,10),Q(2,10),Q(3,10),Q(4,10))
        assert sum(Q(factorial(j),prod(factorial(int(n)) for n in row))*prod(p**int(n) for p,n in zip(theta,row)) for row in counts)==1
    for i in range(4):
        theta=np.eye(4)[i];mixed=mix(data,theta)
        for j,(counts,matrices,_) in enumerate(data):
            expected=np.zeros(4,dtype=int);expected[i]=j
            index=np.flatnonzero((counts==expected).all(axis=1))[0]
            assert np.allclose(mixed[j],matrices[index],rtol=1e-13,atol=0)
    print('Categorical local mixtures: exact normalization and all atomic limits passed',flush=True)


def state_test(data,inputs,tilt):
    """Independent direct-state checks for empty and one-window inputs."""
    low,high,words,spectrum=inputs
    levels=sorted(spectrum);expansion=np.bitwise_count(low)+np.bitwise_count(high)
    states=sorted(set([1,17]+[int(np.flatnonzero(expansion==v)[0]) for v in levels]))
    m=(1<<19)-1;checks=0
    for b in range(5):
        j=int(b!=0);counts,matrices,_=data[j]
        shape=np.zeros(4,dtype=int)
        if j:shape[b-1]=1
        t=matrices[np.flatnonzero((counts==shape).all(axis=1))[0]]
        entries=words[b] if j else [(0,0,0)]
        for state in states:
            actual=np.zeros(9);density={}
            for syndrome,lo,hi in entries:
                w=(int(low[state])^lo).bit_count()+(int(high[state])^hi).bit_count()
                factor=exp(-tilt*w)/len(entries);target=state^syndrome
                if target:
                    actual[M]+=.25*factor
                    density[target]=density.get(target,0)+.25*factor
                else:actual[Z]+=.25*factor
                if syndrome:actual[Z]+=.75*factor/m
                for k,v in enumerate(levels):actual[U+k]+=.75*factor*spectrum[v]/m
            actual[C]=max(density.values(),default=0.)
            assert (actual<=t[M]+t[C]+1e-12).all(),(j,b,state,actual,t[M]+t[C])
            checks+=9
    print('Categorical direct-state empty/single-window checks:',checks,'passed',flush=True)
