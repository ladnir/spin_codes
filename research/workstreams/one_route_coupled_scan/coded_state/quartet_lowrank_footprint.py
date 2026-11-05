"""Positive130-state factorization of the forced-footprint upper kernel.

Exact toy checks; full-width diagnostics are FLOAT, not certificates.
Fixed quartet construction. One bounded/file-free Peach job at a time.
The loader may inject uniform_route_moment from the prior global helper.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb,log2,lgamma
import os,resource,shutil,sys,threading,time
import numpy as np

LD=np.longdouble
N=1<<21;W=943718;T=4096


def mulpoly(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def exact_checks():
    m=F(3,4);b=F(5,4);x=F(2,3);checks=0;path_checks=0;eligible_checks=0;renewal_checks=0;minority_checks=0
    for width in (4,8):
        groups=width//4;dim=groups+2;star=groups+1
        states=[(f,g) for g in range(0,width+1,4) for f in range(width-g+1)]
        index={s:i for i,s in enumerate(states)};zero=index[0,0]
        U=[[F(0)]*dim for _ in states]
        V=[[F(0)]*len(states) for _ in range(dim)]
        powers=[[F(1)]]
        for _ in range(groups):powers.append(mulpoly(powers[-1],[1,4,6,4]))
        for si,(f,g) in enumerate(states):
            if f+g==0:U[si][star]=1
            else:
                for k in range(groups+1):
                    if 4*k<f:continue
                    coeff=sum(F(comb(g,j)*comb(width-f-g,4*k-f-j))*b**j
                              for j in range(g+1) if 0<=4*k-f-j<=width-f-g)
                    U[si][k]=m**f*coeff/comb(width,4*k)
            for k in range(groups+1):
                if g==4*k and f<len(powers[groups-k]):
                    V[k][si]=comb(groups,k)*powers[groups-k][f]*x**(f+g)
            if g==0:V[star][si]=comb(width,f)*x**f
        C=[[sum(V[i][s]*U[s][j] for s in range(len(states))) for j in range(dim)] for i in range(dim)]
        transfer=[[sum(U[s][i]*V[i][d] for i in range(dim)) for d in range(len(states))] for s in range(len(states))]
        parts=[(15,)] if width==4 else []
        if width==8:
            for rest in combinations(range(1,8),3):
                first=1+sum(1<<i for i in rest);parts.append((first,255^first))
        for si,(f,g) in enumerate(states):
            forced=(1<<f)-1;optional=((1<<g)-1)<<f;truth=[F(0)]*len(states)
            for partition in parts:
                for covered in range(1<<width):
                    if any((forced&quad) and (covered&quad)!=quad for quad in partition):continue
                    if f+g==0:out=(covered.bit_count(),0);hits=0
                    else:
                        full=[quad for quad in partition if covered&quad==quad]
                        out=(covered.bit_count()-4*len(full),4*len(full))
                        hits=sum((optional&quad).bit_count() for quad in full)
                    truth[index[out]]+=m**f*b**hits*x**covered.bit_count()/len(parts)
            assert truth==transfer[si]
            checks+=len(states)
        terminal=[m**f*b**g for f,g in states];h=[F(int(i in (groups,star))) for i in range(dim)]
        assert [sum(U[s][i]*h[i] for i in range(dim)) for s in range(len(states))]==terminal
        r=[F(int(s==zero)) for s in range(len(states))];rr=[F(int(i==star)) for i in range(dim)]
        for _ in range(4):
            r=[sum(r[s]*transfer[s][d] for s in range(len(states))) for d in range(len(states))]
            rr=[sum(rr[i]*C[i][j] for i in range(dim)) for j in range(dim)]
            assert sum(r[i]*terminal[i] for i in range(len(states)))==sum(rr[i]*h[i] for i in range(dim))
            path_checks+=1
        float_matrix=matrix(width,m,b,x)
        for i in range(dim):
            for j in range(dim):
                truth=LD(C[i][j].numerator)/LD(C[i][j].denominator)
                assert abs(float_matrix[i,j]-truth)<=LD('1e-16')*max(1,abs(truth))
        p=F(1,3);rprob=F(2,5)
        Bh=[[V[k][s]*(p**states[s][0]*rprob**states[s][1] if k<star else p**states[s][0])
             for s in range(len(states))] for k in range(dim)]
        Ch=[[sum(Bh[i][s]*U[s][j] for s in range(len(states))) for j in range(dim)] for i in range(dim)]
        floating=matrix(width,m,b,x,p,rprob)
        for i in range(dim):
            for j in range(dim):
                truth=LD(Ch[i][j].numerator)/LD(Ch[i][j].denominator)
                assert abs(floating[i,j]-truth)<=LD('1e-16')*max(1,abs(truth))
                eligible_checks+=1
        minority=F(2,3)
        Bs=[[Bh[k][s]*minority**min(sum(states[s]),width-sum(states[s]))
             for s in range(len(states))] for k in range(dim)]
        Cs=[[sum(Bs[i][s]*U[s][j] for s in range(len(states))) for j in range(dim)] for i in range(dim)]
        weighted=matrix(width,m,b,x,p,rprob,float(log2(minority)))
        neutral=minority_matrix(width,LD(m.numerator)/m.denominator,LD(b.numerator)/b.denominator,
                                LD(x.numerator)/x.denominator,LD(p.numerator)/p.denominator,
                                LD(rprob.numerator)/rprob.denominator,LD(1))
        for i in range(dim):
            for j in range(dim):
                truth=LD(Cs[i][j].numerator)/LD(Cs[i][j].denominator)
                assert abs(weighted[i,j]-truth)<=LD('1e-15')*max(1,abs(truth))
                assert abs(neutral[i,j]-floating[i,j])<=LD('1e-15')*max(1,abs(floating[i,j]))
                minority_checks+=1
        positive,scale=positive_region(float_matrix)
        direct=np.linalg.matrix_power(float_matrix,16)
        direct[0,-1]-=1;direct[-1,-1]-=1
        assert np.max(abs(positive*np.exp2(scale)-direct))<=LD('1e-14')*max(1,direct.max())
        Q=[row.copy() for row in C];Q[0][star]-=1;Q[star][star]-=1
        paths={(0,star):F(1)};regions=6
        for _ in range(regions):
            after={}
            for (occupied,i),weight in paths.items():
                if i in (0,star):after[occupied,star]=after.get((occupied,star),F(0))+weight
                for j in range(dim):
                    key=(occupied+1,j);after[key]=after.get(key,F(0))+weight*Q[i][j]
            paths=after
        qfloat=np.array([[LD(v.numerator)/LD(v.denominator) for v in row] for row in Q])
        for minimum in (1,3,6):
            expected=sum(weight*h[i] for (occupied,i),weight in paths.items() if occupied>=minimum)
            truth=np.log2(LD(expected.numerator)/LD(expected.denominator))
            assert abs(occupied_tail(qfloat,0,regions,minimum)-truth)<LD('1e-14')
            renewal_checks+=1
    print(dict(exact_factorization_entries=checks,exact_composed_paths=path_checks,
               eligible_rank_substitution_comparisons=eligible_checks,
               rank_one_renewal_comparisons=renewal_checks,
               minority_factorization_comparisons=minority_checks,
               full_code_certificate=False),flush=True)


def matrix(width,m,b,x,forced_probability=1.,optional_probability=1.,minority_log=0.):
    """Build B_x A directly by univariate coefficient extraction."""
    m,b,x=map(lambda z:LD(z.numerator)/LD(z.denominator) if isinstance(z,F) else LD(z),(m,b,x))
    if minority_log:
        return minority_matrix(width,m,b,x,LD(forced_probability),LD(optional_probability),
                               np.exp2(LD(minority_log)))
    forced_x=x*LD(forced_probability);optional_x=x*LD(optional_probability)
    groups=width//4;dim=groups+2;star=groups+1;s=m*forced_x
    local=np.array([1,4*(1+s),6*(1+s)**2,4*(1+s)**3,1+4*s+6*s*s+4*s**3],dtype=LD)
    baseline=np.array([1,4,6,4,1],dtype=LD)
    powers=[np.array([1],dtype=LD)];diff=np.array([0],dtype=LD);base=powers[0]
    for _ in range(groups):
        # Positive recurrence for Q^n-(1+t)^(4n); no large subtraction.
        diff=np.convolve(diff,local)+np.convolve(base,local-baseline)
        base=np.convolve(base,baseline)
        powers.append(np.convolve(powers[-1],local))
    denom=np.array([LD(comb(width,4*j)) for j in range(groups+1)])
    result=np.zeros((dim,dim),dtype=LD)
    result[0,:groups+1]=diff[::4]/denom
    for i in range(1,groups+1):
        factor=np.array([LD(comb(4*i,j))*optional_x**(4*i)*b**j for j in range(4*i+1)])
        coeff=np.convolve(factor,powers[groups-i])
        result[i,:groups+1]=LD(comb(groups,i))*coeff[::4]/denom
    result[0,star]=1
    result[star,:groups+1]=np.expm1(np.arange(groups+1,dtype=LD)*4*np.log1p(s))
    result[star,star]=1
    return result


def minority_matrix(width,m,b,x,forced_probability,optional_probability,minority):
    """Insert the actual footprint minority in B, using positive binomial transforms."""
    groups=width//4;dim=groups+2;star=groups+1
    forced=m*x*forced_probability;optional=x*optional_probability
    factorial=np.ones(width+1,dtype=LD)
    for j in range(1,width+1):factorial[j]=factorial[j-1]*j
    inverse=1/factorial
    p3=[np.ones(1,dtype=LD)]
    for _ in range(groups):p3.append(np.convolve(p3[-1],np.array([1,4,6,4],dtype=LD)))
    denom=np.array([LD(comb(width,4*j)) for j in range(groups+1)])
    result=np.zeros((dim,dim),dtype=LD)
    for i in range(groups+1):
        n=width-4*i;raw=p3[groups-i];ff=np.arange(len(raw))
        weights=raw*forced**ff*minority**np.minimum(ff+4*i,width-ff-4*i)
        weights/=np.array([LD(comb(n,int(f))) for f in ff])
        if i==0:weights[0]=0 # Physical empty state uses the separate star column.
        trans=np.convolve(weights*inverse[:len(raw)],inverse[:n+1])[:n+1]
        trans*=factorial[n]/factorial[n-np.arange(n+1)]
        factor=np.array([LD(comb(4*i,j))*b**j for j in range(4*i+1)])
        result[i,:groups+1]=LD(comb(groups,i))*optional**(4*i)*np.convolve(factor,trans)[::4]/denom
    ff=np.arange(width+1)
    weights=forced**ff*minority**np.minimum(ff,width-ff);weights[0]=0
    trans=np.convolve(weights*inverse,inverse)[:width+1]*factorial
    result[star,:groups+1]=trans[::4]
    result[0,star]=result[star,star]=1
    return result


def fresh_birth_matrix(width,m,b,p,r,coverage=1.):
    """Coverage-size upper kernel retaining every forced birth's eligible-row cost.

    Drops forced-future containment, unlike matrix(). The physical basis is
    independent of the quota tilt p. This is a separate rigorous relaxation.
    """
    m,b,p,r,coverage=map(LD,(m,b,p,r,coverage))
    groups=width//4;a=p*m/b;delta=(a/r)**4
    local=np.array([delta,4*b,6*b*b,4*b**3,b**4],dtype=LD)
    powers=[np.ones(1,dtype=LD)];p3=[np.ones(1,dtype=LD)]
    for _ in range(groups):
        powers.append(np.convolve(powers[-1],local))
        p3.append(np.convolve(p3[-1],np.array([1,4,6,4],dtype=LD)))
    denom=np.array([LD(comb(width,u)) for u in range(width+1)])
    U=np.empty((width+1,groups+1),dtype=LD)
    V=np.zeros((groups+1,width+1),dtype=LD)
    for k in range(groups+1):
        rest=width-4*k
        U[:,k]=np.convolve(powers[k],np.array([LD(comb(rest,u)) for u in range(rest+1)]))/denom
        ff=np.arange(len(p3[groups-k]))
        V[k,4*k:4*k+len(ff)]=LD(comb(groups,k))*(r*coverage)**(4*k)*p3[groups-k]*(a*coverage)**ff
    return U@V


def fresh_birth_checks():
    m=F(3,4);b=F(5,4);p=F(1,3);r=F(2,5);x=F(2,3);checks=0
    for width in (4,8):
        parts=[(15,)] if width==4 else []
        if width==8:
            for rest in combinations(range(1,8),3):
                first=1+sum(1<<i for i in rest);parts.append((first,255^first))
        floating=fresh_birth_matrix(width,float(m),float(b),float(p),float(r),float(x))
        for u in range(width+1):
            old=(1<<u)-1;truth=[F(0)]*(width+1)
            for partition in parts:
                for covered in range(1<<width):
                    v=covered.bit_count();forced=0;allowed=0
                    for quad in partition:
                        hits=(covered&quad).bit_count()
                        if hits==4:
                            incoming=(old&quad).bit_count();allowed+=incoming
                            if incoming==0:forced+=4
                        else:forced+=hits
                    truth[v]+=x**v*p**forced*r**(v-forced)*b**allowed*(m/b)**forced/len(parts)
            for v in range(width+1):
                exact=LD(truth[v].numerator)/LD(truth[v].denominator)
                assert abs(floating[u,v]-exact)<LD('1e-14')*max(1,abs(exact))
                checks+=1
        Q,scale=physical_positive_region(floating)
        direct=np.linalg.matrix_power(floating,16);direct[:,0]-=1
        assert np.max(abs(Q*np.exp2(scale)-direct))<LD('1e-13')*max(1,direct.max())
    print(dict(fresh_birth_partition_entries=checks,full_code_certificate=False),flush=True)


def physical_positive_region(C):
    """Sixteen updates minus Z=ones*e0^T, through positive arithmetic."""
    Q=C.copy();Q[:,0]=0;scale=LD(0)
    for _ in range(4):
        norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
        common=max(scale,2*scale);new=(Q@Q)*np.exp2(2*scale-common)
        linear=np.exp2(scale-common)
        new+=linear*Q[0,:][None,:]
        new[:,0]+=linear*Q.sum(axis=1)
        Q=new;scale=common
    norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    return Q,scale


def physical_occupied_tail(Q,scale,b,regions=256,minimum=129):
    row=np.zeros(len(Q),dtype=LD);row[0]=1;rs=LD(0)
    terminal=LD(b)**np.arange(len(Q));aa=np.empty(regions+1,dtype=LD);bb=aa.copy();ln2=np.log(LD(2))
    for k in range(regions+1):
        a=row.sum();v=row@terminal
        aa[k]=np.log2(a)+rs if a else -np.inf
        bb[k]=np.log2(v)+rs if v else -np.inf
        if k<regions:
            row=row@Q;norm=row.sum();assert norm>0;row/=norm;rs+=np.log2(norm)
    poly=bb.copy();logs=[]
    for zero in range(regions-minimum+1):
        occupied=regions-zero;logs.append(poly[occupied]+scale*occupied)
        after=np.full(regions+1,LD('-inf'))
        for j in range(regions+1):
            after[j:]=np.logaddexp(after[j:]*ln2,(poly[:regions+1-j]+aa[j])*ln2)/ln2
        poly=after
    largest=max(logs)
    return largest+np.log2(sum(np.exp2(v-largest) for v in logs))


def fresh_birth_value(h,theta,logz,tilts=5):
    rows=8192;a=LD(2)**(32*theta);actual=LD(h)/rows
    m=a*LD(1+LD(2)**-32)**theta*LD(2)**logz;b=(1+m**(1/theta))**theta
    choose=(lgamma(rows+1)-lgamma(h+1)-lgamma(rows-h+1))/float(np.log(LD(2)))
    low=actual/(a-(a-1)*actual);logs=None
    for fraction in np.linspace(0.,1.,tilts):
        odds=(low/(1-low))**(1-fraction)*(actual/(1-actual))**fraction;p=odds/(1+odds)
        logpi=choose+h*np.log2(p)+(rows-h)*np.log2(1-p);r=p+(1-p)/a
        C=fresh_birth_matrix(512,m,b,p,r);Q,scale=physical_positive_region(C)
        with np.errstate(divide='ignore'):candidate=np.log2(Q)+scale-logpi
        logs=candidate if logs is None else np.minimum(logs,candidate)
    scale=logs.max();Q=np.exp2(logs-scale);norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    value=physical_occupied_tail(Q,scale,b)
    return float(theta-32*theta-W*logz+choose-32*theta*128*h+value)


def fresh_birth_diagnostic():
    for h,theta,logz in ((4096,.0625,-1.),(2048,.0625,-.5),(512,.125,-.25)):
        start=time.monotonic();value=fresh_birth_value(h,theta,logz)
        print(dict(fresh_birth_rows=h,theta=theta,logz=logz,tilts=5,log2_diagnostic=value,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def coarse_fresh_prepare(width,m,b,limits):
    """Common A for a guaranteed fresh-quartet count, plus sparse B rows.

    An old footprint u<=upper meets at most upper fully covered new quartets.
    Select max(0,k-upper) genuinely fresh ones and retain all their forced
    positions. This strengthens the130-state kernel without changing setup.
    """
    m,b=map(LD,(m,b));groups=width//4
    limits=tuple(sorted(v for v in limits if 0<v<groups))+(width,)
    latent=[(k,k) for k in range(groups+1)]
    latent += [(k,upper) for upper in limits[:-1] for k in range(upper+1,groups+1)]
    lookup={pair:i for i,pair in enumerate(latent)};star=len(latent)
    dispatch=np.array([[lookup[k,min(k,upper)] for k in range(groups+1)] for upper in limits])
    choose=[np.array([LD(comb(n,j)) for j in range(n+1)]) for n in range(width+1)]
    denominator=choose[width][::4]
    bpolys={o:choose[o]*b**np.arange(o+1) for o in range(0,width+1,4)}
    p3=[np.ones(1,dtype=LD)]
    for _ in range(groups):p3.append(np.convolve(p3[-1],np.array([1,4,6,4],dtype=LD)))
    cache={};rows=[]
    for k,c in latent+[(-1,-1)]:
        if k<0:fs=np.arange(width+1);o=0;coeff=choose[width]
        else:
            fs=np.arange(len(p3[groups-k]))+4*(k-c);o=4*c
            coeff=choose[groups][k]*p3[groups-k]
        vectors=np.zeros((len(fs),groups+1),dtype=LD)
        bins=np.searchsorted(limits,fs+o)
        for i,f in enumerate(fs):
            if f+o==0:continue
            key=(int(f),o)
            if key not in cache:
                poly=np.convolve(bpolys[o],choose[width-f-o])
                dest=np.arange(groups+1)*4-f;valid=dest>=0
                vec=np.zeros(groups+1,dtype=LD)
                vec[valid]=m**f*poly[dest[valid]]/denominator[valid]
                cache[key]=vec
            vectors[i]=cache[key]
        rows.append((fs,o,coeff,vectors,bins))
    terminal=np.array([int(k==groups) for k,c in latent]+[1],dtype=LD)
    return dict(width=width,latent=latent,star=star,dispatch=dispatch,rows=rows,terminal=terminal,
                physical_states=len(cache))


def coarse_fresh_matrix(prepared,p,r,coverage=1.):
    p,r,coverage=map(LD,(p,r,coverage));star=prepared['star']
    C=np.zeros((star+1,star+1),dtype=LD)
    for i,(fs,o,coeff,vectors,bins) in enumerate(prepared['rows']):
        weights=coeff*(p*coverage)**fs*(r*coverage)**o
        for group,columns in enumerate(prepared['dispatch']):
            selected=(bins==group)&(fs+o>0)
            if selected.any():C[i,columns]+=weights[selected]@vectors[selected]
        if o==0 and fs[0]==0:C[i,star]=weights[0]
    return C


def coarse_fresh_checks():
    checks=0;m=F(3,4);b=F(5,4);p=F(1,3);r=F(2,5)
    for width,limits in ((4,()),(8,(1,))):
        prepared=coarse_fresh_prepare(width,float(m),float(b),limits);groups=width//4
        latent=prepared['latent'];star=prepared['star'];dim=star+1;lookup={pair:i for i,pair in enumerate(latent)}
        states=[(f,o) for o in range(0,width+1,4) for f in range(width-o+1)];index={s:i for i,s in enumerate(states)}
        A=[[F(0)]*dim for _ in states];B=[[F(0)]*len(states) for _ in range(dim)]
        parts=[(15,)] if width==4 else []
        if width==8:
            for rest in combinations(range(1,8),3):
                first=1+sum(1<<i for i in rest);parts.append((first,255^first))
        powers=[[F(1)]]
        for _ in range(groups):powers.append(mulpoly(powers[-1],[1,4,6,4]))
        for si,(f,o) in enumerate(states):
            if f+o==0:A[si][star]=1
            else:
                upper=next(v for v in tuple(limits)+(width,) if f+o<=v)
                for k in range(groups+1):
                    if 4*k<f:continue
                    coeff=sum(F(comb(o,j)*comb(width-f-o,4*k-f-j))*b**j
                              for j in range(o+1) if 0<=4*k-f-j<=width-f-o)
                    A[si][lookup[k,min(k,upper)]]=m**f*coeff/comb(width,4*k)
        for i,(k,c) in enumerate(latent):
            for f0,coeff in enumerate(powers[groups-k]):
                f=f0+4*(k-c);o=4*c
                B[i][index[f,o]]=comb(groups,k)*coeff*p**f*r**o
        for f in range(width+1):B[star][index[f,0]]=comb(width,f)*p**f
        for si,(f,o) in enumerate(states):
            forced=(1<<f)-1;old=(1<<(f+o))-1;optional=old^forced
            truth=[F(0)]*len(states)
            upper=next(v for v in tuple(limits)+(width,) if f+o<=v)
            for partition in parts:
                for covered in range(1<<width):
                    if any((forced&quad) and (covered&quad)!=quad for quad in partition):continue
                    if f+o==0:out=(covered.bit_count(),0);hits=0
                    else:
                        full=[quad for quad in partition if covered&quad==quad];k=len(full)
                        j=max(0,k-upper)
                        assert sum(not(old&quad) for quad in full)>=j
                        out=(covered.bit_count()-4*k+4*j,4*(k-j))
                        hits=sum((optional&quad).bit_count() for quad in full)
                    truth[index[out]]+=m**f*b**hits*p**out[0]*r**out[1]/len(parts)
            product=[sum(A[si][i]*B[i][d] for i in range(dim)) for d in range(len(states))]
            assert truth==product;checks+=len(states)
        C=coarse_fresh_matrix(prepared,float(p),float(r))
        for i in range(dim):
            for j in range(dim):
                exact=sum(B[i][s]*A[s][j] for s in range(len(states)))
                value=LD(exact.numerator)/LD(exact.denominator)
                assert abs(C[i,j]-value)<LD('1e-14')*max(1,abs(value))
        terminal=prepared['terminal']
        for si,(f,o) in enumerate(states):
            assert sum(A[si][i]*int(terminal[i]) for i in range(dim))==m**f*b**o
    print(dict(coarse_fresh_partition_entries=checks,terminal_identity=True,full_code_certificate=False),flush=True)


def coarse_fresh_value(h,theta,logz,tilts=5,limits=(16,64,127),logoutside=0.):
    rows=8192;a=LD(2)**(32*theta);actual=LD(h)/rows
    m=a*LD(1+LD(2)**-32)**theta*LD(2)**logz;b=(1+m**(1/theta))**theta
    prepared=coarse_fresh_prepare(512,m,b,limits)
    choose=(lgamma(rows+1)-lgamma(h+1)-lgamma(rows-h+1))/float(np.log(LD(2)))
    low=actual/(a-(a-1)*actual);logs=None
    for fraction in np.linspace(0.,1.,tilts):
        odds=(low/(1-low))**(1-fraction)*(actual/(1-actual))**fraction;p=odds/(1+odds)
        logpi=choose+h*np.log2(p)+(rows-h)*np.log2(1-p)
        r=p+(1-p)*LD(2)**logoutside/a
        C=coarse_fresh_matrix(prepared,p,r);Q,scale=positive_region(C)
        with np.errstate(divide='ignore'):candidate=np.log2(Q)+scale-logpi
        logs=candidate if logs is None else np.minimum(logs,candidate)
    scale=logs.max();Q=np.exp2(logs-scale);norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    value=occupied_tail(Q,scale,terminal=prepared['terminal'])
    return float(theta-32*theta-W*logz+choose-32*theta*128*h+value
                 -(1009255-128*h)*logoutside)


def coarse_fresh_diagnostic():
    for h,theta,logz in ((4096,.0625,-1.),(2048,.0625,-.5),(512,.125,-.25)):
        start=time.monotonic();value=coarse_fresh_value(h,theta,logz)
        print(dict(coarse_fresh_rows=h,theta=theta,logz=logz,tilts=5,dimension=307,log2_diagnostic=value,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def complete_force_matrix(width,m,b,p,r,coverage=1.):
    """An optional full quartet must meet old coverage; other quartets are forced.

    Allows forced full quartets to contain old optional coverage, whose state
    is then constrained to zero. Keeps every genuine forced propagation.
    """
    m,b,p,r,coverage=map(LD,(m,b,p,r,coverage));groups=width//4;s=m*p*coverage
    constant=np.expm1(4*np.log1p(s))
    local=np.array([constant,4*(1+s)**3*b,6*(1+s)**2*b*b,4*(1+s)*b**3,b**4],dtype=LD)
    rows=np.arange(groups+1);prefactor=np.array([LD(comb(groups,i))/LD(comb(width,4*i)) for i in rows])*(r*coverage)**(4*rows)
    C=np.empty((groups+1,groups+1),dtype=LD);power=np.ones(1,dtype=LD)
    for k in range(groups+1):
        rest=width-4*k
        poly=np.convolve(power,np.array([LD(comb(rest,j)) for j in range(rest+1)]))
        C[:,k]=prefactor*poly[::4]
        if k<groups:power=np.convolve(power,local)
    terminal=np.array([LD(comb(groups,i)) for i in rows])*(r*coverage*b)**(4*rows)*(1+s)**(width-4*rows)
    terminal[0]=np.expm1(width*np.log1p(s)) # Positive terminal part, subtract e0 exactly.
    return C,terminal


def complete_force_region(C,terminal_positive):
    P=C.copy();P[0,0]=0;Q=P.copy();scale=LD(0)
    for _ in range(4):
        norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
        common=max(scale,2*scale);new=(Q@Q)*np.exp2(2*scale-common)
        linear=np.exp2(scale-common);new[0,:]+=linear*Q[0,:];new[:,0]+=linear*Q[:,0]
        Q=new;scale=common
    norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    # C^15 (e0+terminal_positive)-e0, without subtracting an empty history.
    end=terminal_positive.copy();norm=end.max();end/=norm;end_scale=np.log2(norm)
    for _ in range(15):
        common=max(LD(0),end_scale)
        end=(C@end)*np.exp2(end_scale-common)+P[:,0]*np.exp2(-common)
        norm=end.max();end/=norm;end_scale=common+np.log2(norm)
    return Q,scale,end,end_scale


def complete_force_tail(Q,scale,end,end_scale,regions=256,minimum=129,profile=False):
    prefix=regions-1;row=np.zeros(len(Q),dtype=LD);row[0]=1;rs=LD(0);ln2=np.log(LD(2))
    aa=np.empty(prefix+1,dtype=LD);bb=aa.copy()
    for k in range(prefix+1):
        a=row[0];v=row@end
        aa[k]=np.log2(a)+rs if a else -np.inf
        bb[k]=np.log2(v)+rs if v else -np.inf
        if k<prefix:
            row=row@Q;norm=row.sum();assert norm>0;row/=norm;rs+=np.log2(norm)
    logs=[];occupancies=[]
    for last_positive in (False,True):
        poly=bb.copy() if last_positive else aa.copy()
        lower=minimum-int(last_positive)
        for zero in range(prefix-lower+1):
            occupied=prefix-zero
            logs.append(poly[occupied]+scale*occupied+(end_scale if last_positive else 0))
            occupancies.append(occupied+int(last_positive))
            after=np.full(prefix+1,LD('-inf'))
            for j in range(prefix+1):
                after[j:]=np.logaddexp(after[j:]*ln2,(poly[:prefix+1-j]+aa[j])*ln2)/ln2
            poly=after
    largest=max(logs)
    weights=np.exp2(np.array(logs)-largest);value=largest+np.log2(weights.sum())
    if not profile:return value
    forward=np.zeros(len(Q),dtype=LD);forward[0]=1
    for _ in range(prefix//2):
        forward=forward@Q;forward/=forward.sum()
    backward=end.copy();backward/=backward.max()
    for _ in range(prefix-prefix//2):
        backward=Q@backward;backward/=backward.max()
    marginal=forward*backward;marginal/=marginal.sum()
    top=np.argsort(marginal)[-6:][::-1]
    full_log=bb[prefix]+scale*prefix+end_scale
    return value,dict(mean_occupied_regions=float(np.dot(weights,occupancies)/weights.sum()),
                      all_regions_fraction=float(np.exp2(full_log-value)),
                      all_positive_midpoint_mean_optional=float(np.dot(marginal,np.arange(len(Q))*4)),
                      all_positive_midpoint_top=[(int(i),float(marginal[i])) for i in top])


def complete_force_checks():
    m=F(3,4);b=F(5,4);p=F(1,3);r=F(2,5);checks=0;composed=0
    for width in (4,8):
        groups=width//4;states=[(f,o) for o in range(0,width+1,4) for f in range(width-o+1)]
        A=[[F(0)]*(groups+1) for _ in states];B=[[F(0)]*len(states) for _ in range(groups+1)]
        parts=[(15,)] if width==4 else []
        if width==8:
            for rest in combinations(range(1,8),3):
                first=1+sum(1<<i for i in rest);parts.append((first,255^first))
        for si,(f,o) in enumerate(states):
            forced=(1<<f)-1;covered=(1<<(f+o))-1;optional=covered^forced
            for k in range(groups+1):
                for partition in parts:
                    for selected in combinations(partition,k):
                        union=0
                        for quad in selected:union|=quad
                        if forced&~union or any(not(covered&quad) for quad in selected):continue
                        A[si][k]+=m**f*b**((optional&union).bit_count())/(len(parts)*comb(groups,k))
                if o==4*k:B[k][si]=comb(groups,k)*comb(width-o,f)*p**f*r**o
        exact=[[sum(B[i][s]*A[s][j] for s in range(len(states))) for j in range(groups+1)] for i in range(groups+1)]
        C,gpos=complete_force_matrix(width,float(m),float(b),float(p),float(r))
        for i in range(groups+1):
            for j in range(groups+1):
                value=LD(exact[i][j].numerator)/LD(exact[i][j].denominator)
                assert abs(C[i,j]-value)<LD('1e-14')*max(1,abs(value));checks+=1
        terminal=[m**f*b**o for f,o in states]
        g=[sum(B[i][s]*terminal[s] for s in range(len(states))) for i in range(groups+1)]
        physical=[[sum(A[s][i]*B[i][d] for i in range(groups+1)) for d in range(len(states))] for s in range(len(states))]
        if width==8:
            singleton=states.index((1,0))
            # The earlier relaxation admitted1->8->1. Fresh forced positions forbid it.
            assert sum(physical[singleton][s]*physical[s][singleton]
                       for s,(f,o) in enumerate(states) if f+o==width)==0
        row=[F(int(s==0)) for s in range(len(states))];small=[F(int(i==0)) for i in range(groups+1)]
        for t in range(1,5):
            row=[sum(row[s]*physical[s][d] for s in range(len(states))) for d in range(len(states))]
            assert sum(row[s]*terminal[s] for s in range(len(states)))==sum(small[i]*g[i] for i in range(groups+1))
            small=[sum(small[i]*exact[i][j] for i in range(groups+1)) for j in range(groups+1)];composed+=1
        Q,scale,end,end_scale=complete_force_region(C,gpos)
        direct=np.linalg.matrix_power(C,16);direct[0,0]-=1
        assert np.max(abs(Q*np.exp2(scale)-direct))<LD('1e-13')*max(1,direct.max())
        gfloat=gpos.copy();gfloat[0]+=1
        direct_end=np.linalg.matrix_power(C,15)@gfloat;direct_end[0]-=1
        assert np.max(abs(end*np.exp2(end_scale)-direct_end))<LD('1e-13')*max(1,direct_end.max())
        # Exhaust the occupied-region count for three regions, including the final block.
        P=Q*np.exp2(scale);d=end*np.exp2(end_scale)
        paths={(0,0):LD(1)}
        for _ in range(2):
            after={}
            for (occupied,i),weight in paths.items():
                if i==0:after[occupied,0]=after.get((occupied,0),0)+weight
                for j in range(groups+1):
                    key=(occupied+1,j);after[key]=after.get(key,0)+weight*P[i,j]
            paths=after
        for minimum in (1,2,3):
            truth=sum(weight*((int(i==0) if occupied>=minimum else 0)+(d[i] if occupied+1>=minimum else 0))
                      for (occupied,i),weight in paths.items())
            result=complete_force_tail(Q,scale,end,end_scale,regions=3,minimum=minimum)
            assert abs(result-np.log2(truth))<LD('1e-12')
    print(dict(complete_force_factorization_entries=checks,exact_composed_paths=composed,
               terminal_and_occupied_region_checks=True,full_code_certificate=False),flush=True)


def complete_force_value(h,theta,logz,tilts=9,logoutside=0.,logcoverage=0.,logrankgap=0.,
                         logavailable=0.,profile=False,loghcoverage=0.,wide_tilts=False,logforced=0.,
                         loghlower=0.,birth_run_bonus=False):
    rows=8192;a=LD(2)**(32*theta);actual=LD(h)/rows
    m=a*LD(1+LD(2)**-32)**theta*LD(2)**logz;b=(1+m**(1/theta))**theta
    run_log=min(LD(0),np.log2(b)-theta*np.log2(LD(2)**32+1)) if birth_run_bonus else LD(0)
    m*=LD(2)**run_log
    assert logrankgap>=0
    m*=LD(2)**(logrankgap+logavailable);b*=LD(2)**(logrankgap+logavailable)
    choose=(lgamma(rows+1)-lgamma(h+1)-lgamma(rows-h+1))/float(np.log(LD(2)))
    low=actual/(a-(a-1)*actual);logs=None;endlogs=None
    fractions=list(np.linspace(0.,1.,tilts))
    candidates=[]
    for fraction in fractions:
        odds=(low/(1-low))**(1-fraction)*(actual/(1-actual))**fraction;p=odds/(1+odds)
        candidates.append(p)
    if wide_tilts:
        for shift in (-16,-12,-8,-6,-4,-2,2,4):
            odds=actual/(1-actual)*LD(2)**shift;candidates.append(odds/(1+odds))
    for p in candidates:
        logpi=choose+h*np.log2(p)+(rows-h)*np.log2(1-p)
        forced=p*LD(2)**(loghcoverage+loghlower)
        r=forced+(1-p)*LD(2)**(logoutside-logrankgap)/a
        C,gpos=complete_force_matrix(512,m,b,forced*LD(2)**logforced,r,LD(2)**logcoverage)
        Q,scale,end,end_scale=complete_force_region(C,gpos)
        if wide_tilts and (not np.all(Q>0) or not np.all(end>0)):
            raise ArithmeticError('Wide-tilt underflow: do not use zero entries as upper bounds')
        with np.errstate(divide='ignore'):
            candidate=np.log2(Q)+scale-logpi;candidate_end=np.log2(end)+end_scale-logpi
        logs=candidate if logs is None else np.minimum(logs,candidate)
        endlogs=candidate_end if endlogs is None else np.minimum(endlogs,candidate_end)
    scale=logs.max();Q=np.exp2(logs-scale);norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    end_scale=endlogs.max();end=np.exp2(endlogs-end_scale)
    result=complete_force_tail(Q,scale,end,end_scale,profile=profile)
    value=result[0] if profile else result
    bound=float(theta-32*theta-W*logz+choose-32*theta*128*h+value
                -(1009255-128*h)*logoutside-1736703*logcoverage+(65536-128*h)*logrankgap
                -256*h*loghcoverage-129*h*loghlower-512*run_log)
    return (bound,result[1]) if profile else bound


def complete_force_profile():
    """Coarse sensitivity NEGATIVE CONTROL: phase switches prevent joint means."""
    step=1/1024
    for h,theta,logz,gap in ((4096,1/16,-1.,.5),(2048,1/16,-.5,0.),(512,1/8,-.25,0.)):
        start=time.monotonic();base,detail=complete_force_value(h,theta,logz,logrankgap=gap,profile=True)
        common=dict(logrankgap=gap)
        coverage=complete_force_value(h,theta,logz,logcoverage=-step,**common)
        outside=complete_force_value(h,theta,logz,logoutside=-step,**common)
        available=complete_force_value(h,theta,logz,logavailable=-step,**common)
        # Left sensitivities of the upper envelope, not statistics of a sampled code.
        Fmean=1736703+(base-coverage)/step
        Rmean=1009255+(base-outside)/step
        Lmean=(base-available)/step
        print(dict(pressure_rows=h,theta=theta,logz=logz,logrankgap=gap,log2_diagnostic=base,
                   left_sensitivity_F=Fmean,left_sensitivity_R=Rmean,left_sensitivity_L_eff=Lmean,
                   left_sensitivity_R_minus_L=Rmean-Lmean,profile=detail,
                   coarse_sensitivities_are_not_joint_means=True,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def complete_force_capacity():
    for h,theta,logz,gap in ((512,1/8,-.25,0.),(2048,1/16,-.5,0.),(4096,1/16,-1.,.5)):
        start=time.monotonic();values=[]
        for marker in (0.,-.25,-.5,-1.,-2.):
            values.append((marker,complete_force_value(h,theta,logz,logrankgap=gap,loghcoverage=marker)))
        print(dict(capacity_rows=h,theta=theta,logz=logz,logrankgap=gap,markers_and_log2_diagnostics=values,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def complete_force_wide():
    for h,theta,logz,gap in ((512,1/8,-.25,0.),(2048,1/16,-.5,0.),(4096,1/16,-1.,.5)):
        start=time.monotonic();values=[]
        for marker,outside in ((0.,0.),(-1/1024,0.),(-1/64,0.),(-1/16,0.),(0.,-1/16)):
            value=complete_force_value(h,theta,logz,logrankgap=gap,loghcoverage=marker,
                                       logoutside=outside,wide_tilts=True)
            values.append((marker,outside,value))
        print(dict(wide_quota_rows=h,theta=theta,logz=logz,logrankgap=gap,
                   markers_and_log2_diagnostics=values,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def complete_force_sensitivity():
    step=2**-20
    for h,theta,logz,gap in ((512,1/8,-.25,0.),(2048,1/16,-.5,0.),(4096,1/16,-1.,.5)):
        start=time.monotonic();values={}
        common=dict(logrankgap=gap)
        base=complete_force_value(h,theta,logz,**common)
        for name,offset in (('logcoverage',1736703),('logoutside',1009255),
                            ('loghcoverage',256*h),('logavailable',0),('logforced',0)):
            minus=complete_force_value(h,theta,logz,**common,**{name:-step})
            plus=complete_force_value(h,theta,logz,**common,**{name:step})
            values[name]=dict(left=offset+(base-minus)/step,right=offset+(plus-base)/step,
                             central=offset+(plus-minus)/(2*step))
        print(dict(sensitivity_rows=h,step=step,log2_diagnostic=base,values=values,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def complete_force_joint():
    samples=((512,1/8,-.25,0.),(2048,1/16,-.5,0.),(4096,1/16,-1.,.5))
    for h,theta,logz,gap in samples:
        start=time.monotonic();values=[]
        trials=([dict(loghlower=x) for x in (0.,1/64,1/16,1/4,1.)] if h==512 else
                [dict(output_shift=x) for x in (0.,1/64,1/16,1/8,1/4)] if h==2048 else
                [dict(logoutside=y,loghcoverage=x) for x,y in ((0.,0.),(-1/16,-1/16),
                    (-1/4,-1/16),(-1/16,-1/4),(-1/4,-1/4))])
        for original in trials:
            trial=original.copy();output=trial.pop('output_shift',0.)
            value=complete_force_value(h,theta,logz+output,logrankgap=gap,**trial)
            values.append((original,value))
        print(dict(joint_rows=h,theta=theta,logz=logz,logrankgap=gap,
                   parameters_and_log2_diagnostics=values,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def availability_run_checks():
    closed_checks=0;open_checks=0
    for q in (2,4,16,1<<32):
        for x in (F(1),F(1,2),F(1,4),F(1,16),F(1,q)):
            lam=1+(q+1)*x;c=(q+1)*x;e=min(F(1),lam/(q+1))
            # Begin a run with its required birth; H advances an allowed state.
            row=[F(0),q*x]
            for length in range(1,25):
                closed=row[0]+x*row[1];opened=row[0]+row[1]
                assert closed<=e*c*lam**(length-1)
                assert opened<=c*lam**(length-1)
                closed_checks+=1;open_checks+=1
                row=[row[0]+x*row[1],q*x*row[0]+(1+q*x)*row[1]]
    print(dict(exact_forced_start_closed_runs=closed_checks,exact_forced_start_terminal_runs=open_checks,
               full_code_certificate=False),flush=True)


def complete_force_run_bonus():
    availability_run_checks()
    for h,theta,logz,gap in ((2048,1/16,-.4375,0.),(4096,1/16,-1.,.5),(512,1/8,-.25,0.)):
        start=time.monotonic();values=[]
        for output_shift in (0.,-1/16,-1/4):
            value=complete_force_value(h,theta,logz+output_shift,logrankgap=gap,birth_run_bonus=True)
            values.append((output_shift,value))
        print(dict(run_bonus_rows=h,theta=theta,logz=logz,logrankgap=gap,
                   shifts_and_log2_diagnostics=values,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def complete_force_run_search():
    for h in (2048,4096):
        start=time.monotonic();best=(float('inf'),None)
        for theta in (1/32,1/16,3/32,1/8,3/16,1/4):
            rowbest=(float('inf'),None)
            for fraction in (.125,.1875,.25,.3125,.375,.5):
                logz=-32*theta*fraction
                value=complete_force_value(h,theta,logz,tilts=3,birth_run_bonus=True)
                if value<rowbest[0]:rowbest=(value,logz)
                if value<best[0]:best=(value,(theta,logz))
            print(dict(run_search_rows=h,theta=theta,best_for_theta=rowbest,
                       seconds=round(time.monotonic()-start,3),certificate=False),flush=True)
        print(dict(run_search_rows=h,best=best,certificate=False),flush=True)


def singleton_basis(width,limits):
    groups=width//4;limits=tuple(sorted(set(v for v in limits if 0<v<width)))+(width,)
    latent=[(k,0) for k in range(groups+1)]
    for upper in limits[:-1]:
        for k in range(upper//2+1,min(upper,groups)+1):
            pair=(k,2*k-upper)
            if pair not in latent:latent.append(pair)
    lookup={pair:i for i,pair in enumerate(latent)}
    pascal=np.zeros((width+1,width+1),dtype=LD)
    for n in range(width+1):pascal[n,:n+1]=[LD(comb(n,j)) for j in range(n+1)]
    return dict(width=width,groups=groups,limits=limits,latent=latent,lookup=lookup,pascal=pascal)


def singleton_matrix(basis,b,p,r,m=None):
    """Guaranteed singleton propagation, retaining forced state AND output.

    Each selected singleton's outside-H successors are forced. Select on
    canonical terms, then enlarge A using the coverage-bin upper bound.
    """
    b,p,r=map(LD,(b,p,r));m=b if m is None else LD(m)
    width=basis['width'];groups=basis['groups']
    latent=basis['latent'];limits=basis['limits'];lookup=basis['lookup'];pascal=basis['pascal']
    dimension=len(latent);C=np.zeros((dimension,dimension),dtype=LD)
    birth_weights=pascal*(m*p)**np.arange(width+1)[None,:]
    denominator=pascal[width,::4]
    bpowers=[pascal[n,:n+1]*b**np.arange(n+1) for n in range(width+1)]
    terminal=np.zeros(dimension,dtype=LD)
    for row,(i,j) in enumerate(latent):
        n=width-4*i;old_optional=4*(i-j);selected=4*j
        selected_slope=p*b+(r-p)*m
        selected_poly=pascal[selected,:selected+1]*p**(selected-np.arange(selected+1))*selected_slope**np.arange(selected+1)
        factor=pascal[groups,i]*r**old_optional*np.convolve(bpowers[old_optional],selected_poly)
        lower=1
        for upper in limits:
            lo=max(0,lower-4*i);hi=min(n,upper-4*i)
            if lo<=hi:
                summed=birth_weights[:n+1,lo:hi+1].sum(axis=1)
                poly=pascal[n,:n+1]*summed
                coefficients=np.convolve(factor,poly)[::4]/denominator
                for k in range(min(groups,upper)+1):
                    C[row,lookup[k,max(0,2*k-upper)]]+=coefficients[k]
            lower=upper+1
        if i==0:C[row,lookup[0,0]]+=1
        terminal[row]=pascal[groups,i]*(r*b)**old_optional*selected_slope**selected*(1+p*m)**n
    terminal[0]=np.expm1(width*np.log1p(p*m))
    return C,terminal


def singleton_checks():
    m=F(3,4);b=F(5,4);p=F(1,3);outside=F(1,5);r=p+outside;checks=0
    for width,limits in ((4,(1,)),(8,(2,3))):
        basis=singleton_basis(width,limits);groups=width//4;latent=basis['latent'];dim=len(latent)
        states=[(f,o) for o in range(width+1) for f in range(width-o+1)];index={s:i for i,s in enumerate(states)}
        A=[[F(0)]*dim for _ in states];B=[[F(0)]*len(states) for _ in range(dim)]
        for si,(f,o) in enumerate(states):
            if f+o==0:A[si][0]=1;continue
            upper=next(v for v in basis['limits'] if f+o<=v)
            for k in range(min(groups,upper)+1):
                if 4*k<f:continue
                coeff=sum(F(comb(o,h)*comb(width-f-o,4*k-f-h))*b**h
                          for h in range(o+1) if 0<=4*k-f-h<=width-f-o)
                A[si][basis['lookup'][k,max(0,2*k-upper)]]=m**f*coeff/comb(width,4*k)
        for row,(i,j) in enumerate(latent):
            for f0 in range(width-4*i+1):
                for h in range(4*j+1):
                    f=f0+4*j-h;o=4*(i-j)+h
                    B[row][index[f,o]]+=comb(groups,i)*comb(width-4*i,f0)*p**f0*r**(4*(i-j))*comb(4*j,h)*p**h*outside**(4*j-h)
        C,gpos=singleton_matrix(basis,float(b),float(p),float(r),float(m))
        for i in range(dim):
            for j in range(dim):
                exact=sum(B[i][s]*A[s][j] for s in range(len(states)))
                value=LD(exact.numerator)/LD(exact.denominator)
                assert abs(C[i,j]-value)<LD('1e-13')*max(1,abs(value));checks+=1
            exact=sum(B[i][s]*m**states[s][0]*b**states[s][1] for s in range(len(states)))
            actual=gpos[i]+int(i==0)
            assert abs(actual-LD(exact.numerator)/LD(exact.denominator))<LD('1e-13')*max(1,abs(actual))
    print(dict(singleton_factorization_entries=checks,terminal_identity=True,
               full_code_certificate=False),flush=True)


def singleton_value(h,theta,logz,tilts=9,limits=(64,128,192,255),loghlower=0.,
                    loghcoverage=0.,logoutside=0.,logrankgap=0.):
    rows=8192;a=LD(2)**(32*theta);actual=LD(h)/rows
    m=a*LD(1+LD(2)**-32)**theta*LD(2)**logz;b=(1+m**(1/theta))**theta
    b*=LD(2)**logrankgap;m*=LD(2)**logrankgap;basis=singleton_basis(512,limits)
    choose=(lgamma(rows+1)-lgamma(h+1)-lgamma(rows-h+1))/float(np.log(LD(2)))
    low=actual/(a-(a-1)*actual);logs=None;endlogs=None
    for fraction in np.linspace(0.,1.,tilts):
        odds=(low/(1-low))**(1-fraction)*(actual/(1-actual))**fraction;p=odds/(1+odds)
        logpi=choose+h*np.log2(p)+(rows-h)*np.log2(1-p)
        eligible=p*LD(2)**(loghlower+loghcoverage)
        r=eligible+(1-p)*LD(2)**(logoutside-logrankgap)/a
        C,gpos=singleton_matrix(basis,b,eligible,r,m)
        Q,scale,end,end_scale=complete_force_region(C,gpos)
        assert np.all(Q>0) and np.all(end>0)
        candidate=np.log2(Q)+scale-logpi;candidate_end=np.log2(end)+end_scale-logpi
        logs=candidate if logs is None else np.minimum(logs,candidate)
        endlogs=candidate_end if endlogs is None else np.minimum(endlogs,candidate_end)
    scale=logs.max();Q=np.exp2(logs-scale);norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    end_scale=endlogs.max();end=np.exp2(endlogs-end_scale)
    value=complete_force_tail(Q,scale,end,end_scale)
    return float(theta-32*theta-W*logz+choose-32*theta*128*h+value
                 -129*h*loghlower-256*h*loghcoverage-(1009255-128*h)*logoutside
                 +(65536-128*h)*logrankgap)


def singleton_diagnostic():
    for h,theta,logz,gap in ((2048,1/16,-.4375,0.),(4096,1/16,-1.,.5),(512,1/8,-.25,0.)):
        start=time.monotonic();value=singleton_value(h,theta,logz,logrankgap=gap)
        print(dict(singleton_rows=h,theta=theta,logz=logz,logrankgap=gap,
                   dimension=258,log2_diagnostic=value,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def rankgap_certificate():
    gap=65537;q=1<<32
    # (1+2/q)^N <= 1/(1-2N/q), by a positive geometric majorant.
    correction=F(1024,989)/(1-F(2*N,q))
    assert N-32*gap-32==-64 and correction<2
    print(dict(global_rank_minus_available_slots=gap,exact_probability_upper=f'{correction} * 2^-64',
               failure_bits_strictly_greater_than=63,full_code_certificate=False),flush=True)


def complete_force_rankgap():
    rankgap_certificate()
    for h,theta,logz in ((4096,.0625,-1.),(2048,.0625,-.5),(512,.125,-.25)):
        start=time.monotonic();values=[]
        for marker in (0.,1/16,1/4,1/2,1.):
            values.append((marker,complete_force_value(h,theta,logz,logrankgap=marker)))
        print(dict(rankgap_rows=h,theta=theta,logz=logz,markers_and_log2_diagnostics=values,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def complete_force_diagnostic():
    for h,theta,logz in ((4096,.0625,-1.),(2048,.0625,-.5),(512,.125,-.25)):
        start=time.monotonic();value=complete_force_value(h,theta,logz)
        print(dict(complete_force_rows=h,theta=theta,logz=logz,tilts=9,dimension=129,log2_diagnostic=value,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def complete_force_markers():
    h=4096;theta=1/16;logz=-1.
    for outside,coverage,output in ((0.,0.,0.),(0.,-1/256,0.),(-1/256,0.,0.),
                                    (0.,0.,-1/256),(0.,-1/16,0.),(-1/16,0.,0.),
                                    (0.,-1/4,0.),(-1/4,0.,0.)):
        start=time.monotonic();value=complete_force_value(h,theta,logz+output,logoutside=outside,logcoverage=coverage)
        print(dict(complete_force_rows=h,theta=theta,logz=logz+output,logoutside=outside,logcoverage=coverage,
                   log2_diagnostic=value,seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def power_log(C,steps):
    norm=C.sum(axis=1).max();base=C/norm;scale=np.log2(norm)
    row=np.zeros(len(C),dtype=LD);row[-1]=1;row_scale=LD(0)
    while steps:
        if steps&1:
            row=row@base;norm=row.sum();assert norm>0
            row/=norm;row_scale+=scale+np.log2(norm)
        steps>>=1
        if steps:
            base=base@base;norm=base.sum(axis=1).max();assert norm>0
            base/=norm;scale=2*scale+np.log2(norm)
    return row_scale+np.log2(row[-1]+row[-2])


def positive_region(C):
    """C^16-Z, Z=(e0+e*)e*^T, with positive arithmetic."""
    Q=C.copy();Q[0,-1]=0;Q[-1,-1]=0;scale=LD(0)
    for _ in range(4):
        norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
        common=max(scale,2*scale);new=(Q@Q)*np.exp2(2*scale-common)
        linear=np.exp2(scale-common)
        new[0,:]+=linear*Q[-1,:];new[-1,:]+=linear*Q[-1,:]
        new[:,-1]+=linear*(Q[:,0]+Q[:,-1])
        Q=new;scale=common
    norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    return Q,scale


def occupied_tail(Q,scale,regions=256,minimum=129,terminal=None):
    row=np.zeros(len(Q),dtype=LD);row[-1]=1;rs=LD(0)
    aa=np.empty(regions+1,dtype=LD);bb=aa.copy();ln2=np.log(LD(2))
    for k in range(regions+1):
        a=row[0]+row[-1];b=row[-2]+row[-1] if terminal is None else row@terminal
        aa[k]=np.log2(a)+rs if a else -np.inf
        bb[k]=np.log2(b)+rs if b else -np.inf
        if k<regions:
            row=row@Q;norm=row.sum();assert norm>0;row/=norm;rs+=np.log2(norm)
    poly=bb.copy();logs=[]
    for zero in range(regions-minimum+1):
        occupied=regions-zero;logs.append(poly[occupied]+scale*occupied)
        after=np.full(regions+1,LD('-inf'))
        for j in range(regions+1):
            after[j:]=np.logaddexp(after[j:]*ln2,(poly[:regions+1-j]+aa[j])*ln2)/ln2
        poly=after
    largest=max(logs)
    return largest+np.log2(sum(np.exp2(v-largest) for v in logs))


def evaluate(theta,logz,logx,occupied=True):
    m=np.exp2(LD(32)*theta+theta*np.log2(1+LD(2)**-32)+logz)
    b=np.exp2(theta*np.log2(1+np.exp2(np.log2(m)/theta)))
    C=matrix(512,m,b,np.exp2(LD(logx)))
    if occupied:
        Q,scale=positive_region(C);value=occupied_tail(Q,scale)
    else:value=power_log(C,T)
    return float(value)


def diagnostic():
    for fraction,theta,logz,logx in ((.125,.2,0.,-6.4),(.25,.125,0.,-4.),
        (.45,.125,0.,-4.),(.5,.125,-1.,-3.),(.625,.125,-2.,-2.),
        (.75,.125,-3.,-1.),(.8,.2,-6.,1.)):
        start=time.monotonic();covered=round(fraction*N)
        count=evaluate(theta,logz,logx)
        route=uniform_route_moment(covered,theta)
        bound=theta-32*theta-W*logz-covered*logx+route+count
        print(dict(coverage=fraction,theta=theta,logz=logz,logx=logx,
                   log2_diagnostic=bound,log2_count=count,log2_route=route,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def eligible_value(h,theta,logz,logeligible=0.,logoutside=0.,tilts=3,
                   remainder=None,logcoverage=0.):
    rows=8192;a=LD(2)**(32*theta);actual=LD(h)/rows
    m=a*LD(1+LD(2)**-32)**theta*LD(2)**logz
    b=(1+m**(1/theta))**theta
    extra=0.;minority_log=0.;outside_cap=128*(rows-h)
    if remainder is not None:
        eta1,eta2=remainder
        assert eta1>=0 and eta2>=0 and logoutside<=0 and logcoverage<=0 and h<=7884
        total=eta1+eta2;constant=7217-4096*h
        extra=eta1*constant+eta2*(constant+11*W)
        # CH >= 129h is an independent marker from the remainder's CH coefficient.
        extra+=129*h*32*total
        logeligible+=32*total;minority_log=-total
        m*=LD(2)**(-11*eta2);b*=LD(2)**(-11*eta2)
        outside_cap=1009255-128*h
    choose=(lgamma(rows+1)-lgamma(h+1)-lgamma(rows-h+1))/float(np.log(LD(2)))
    low=actual/(a-(a-1)*actual);logs=None
    for fraction in (np.linspace(0.,1.,tilts) if h<rows else (0.,)):
        if h==rows:p=LD(1);logpi=0.
        else:
            odds=(low/(1-low))**(1-fraction)*(actual/(1-actual))**fraction
            p=odds/(1+odds)
            logpi=choose+h*np.log2(p)+(rows-h)*np.log2(1-p)
        route=p*LD(2)**logeligible+(1-p)*LD(2)**logoutside/a
        forced=p*LD(2)**logeligible
        C=matrix(512,m,b,LD(2)**logcoverage,forced,route,minority_log)
        Q,scale=positive_region(C)
        with np.errstate(divide='ignore'):candidate=np.log2(Q)+scale-logpi
        logs=candidate if logs is None else np.minimum(logs,candidate)
    scale=logs.max();Q=np.exp2(logs-scale)
    norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    value=occupied_tail(Q,scale)
    return float(theta-32*theta-W*logz+choose-32*theta*128*h+value
                 -129*h*logeligible-outside_cap*logoutside+extra
                 -1736703*logcoverage)


def remainder_diagnostic():
    for h,theta,logz in ((512,.125,-.25),(2048,.0625,-.5),(4096,.0625,-1.)):
        start=time.monotonic();values=[]
        for eta1,eta2 in ((0.,0.),(1/256,0.),(1/64,0.),(0.,1/256),(0.,1/64)):
            value=eligible_value(h,theta,logz,tilts=9,remainder=(eta1,eta2))
            values.append((eta1,eta2,value))
        print(dict(remainder_rows=h,theta=theta,logz=logz,tilts=9,
                   markers_and_log2_diagnostics=values,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def eligible_diagnostic():
    rows=8192
    for h,theta in ((1,.125),(8,.125),(64,.125),(512,.125),(2048,.125),(4096,.125),(6144,.125),(8192,.125)):
        start=time.monotonic();actual=LD(h)/rows
        logz=-32*theta*float(actual)
        bound=eligible_value(h,theta,logz)
        print(dict(eligible_rows=h,theta=theta,logz=logz,log2_diagnostic=float(bound),
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def eligible_search():
    for h in (512,2048,4096):
        best=(float('inf'),None);start=time.monotonic()
        for theta in (1/16,1/8,.2,1/4):
            for fraction in (0.,.25,.5,.75,1.):
                logz=-32*theta*fraction;value=eligible_value(h,theta,logz)
                if value<best[0]:best=(value,(theta,logz,0.,0.))
        # Recheck a small degree-marker neighborhood at the best field tilt.
        theta,logz,_,_=best[1]
        for x in (.25,1.):
            for y in (0.,-.25,-1.):
                value=eligible_value(h,theta,logz,x,y)
                if value<best[0]:best=(value,(theta,logz,x,y))
        print(dict(searched_eligible_rows=h,log2_diagnostic=best[0],parameters=best[1],
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def quota_diagnostic():
    for h,theta,logz in ((512,.125,-.25),(2048,.0625,-.5),(4096,.0625,-1.)):
        start=time.monotonic();base=eligible_value(h,theta,logz,tilts=17)
        shift=1/64
        output=eligible_value(h,theta,logz-shift,tilts=17)
        eligible=eligible_value(h,theta,logz,shift,0.,tilts=17)
        outside=eligible_value(h,theta,logz,0.,-shift,tilts=17)
        print(dict(quota_rows=h,theta=theta,logz=logz,tilts=17,
                   base=base,lower_output_marker=output,eligible_marker=eligible,outside_marker=outside,
                   seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


def main():
    if shutil.disk_usage('/tmp').free<20*(1<<30):raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(2*(1<<30),2*(1<<30)))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    if '--checks' in sys.argv:exact_checks()
    elif '--fresh-checks' in sys.argv:fresh_birth_checks()
    elif '--fresh' in sys.argv:fresh_birth_diagnostic()
    elif '--coarse-checks' in sys.argv:coarse_fresh_checks()
    elif '--coarse' in sys.argv:coarse_fresh_diagnostic()
    elif '--complete-checks' in sys.argv:complete_force_checks()
    elif '--complete' in sys.argv:complete_force_diagnostic()
    elif '--complete-markers' in sys.argv:complete_force_markers()
    elif '--rankgap' in sys.argv:complete_force_rankgap()
    elif '--profile' in sys.argv:complete_force_profile()
    elif '--capacity' in sys.argv:complete_force_capacity()
    elif '--wide' in sys.argv:complete_force_wide()
    elif '--sensitivity' in sys.argv:complete_force_sensitivity()
    elif '--joint' in sys.argv:complete_force_joint()
    elif '--run-checks' in sys.argv:availability_run_checks()
    elif '--run-bonus' in sys.argv:complete_force_run_bonus()
    elif '--run-search' in sys.argv:complete_force_run_search()
    elif '--singleton-checks' in sys.argv:singleton_checks()
    elif '--singleton' in sys.argv:singleton_diagnostic()
    elif '--all-checks' in sys.argv:
        exact_checks();fresh_birth_checks();coarse_fresh_checks();complete_force_checks();rankgap_certificate()
        singleton_checks();availability_run_checks()
    elif '--eligible' in sys.argv:eligible_diagnostic()
    elif '--search' in sys.argv:eligible_search()
    elif '--quota' in sys.argv:quota_diagnostic()
    elif '--remainder' in sys.argv:remainder_diagnostic()
    else:diagnostic()
    timer.cancel()


if __name__=='__main__':main()
