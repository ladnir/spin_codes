"""Actual two-update positive operator on four coordinate paths.

Unlike label lookahead alone, this keeps pair correlations across two
successive matching bits. Construction unchanged; bounded FLOAT proposal.
"""
import time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from layered_general_components import (tables,STATES,M,b,N,envelopes,mixed_dimension)
from check_layered_component_count import log_upper_units

CONTROL_MARKERS=[F(n,d) for n,d in (
    (8081456898423955,1208925819614629174706176),
    (2831158716032241,18014398509481984),
    (6093940992491657,1073741824),(4981563812423823,140737488355328),
    (6503579565136291,1073741824),(2054538321893257,1208925819614629174706176),
    (8086853783136951,1125899906842624),(7651793841660657,140737488355328),
    (7190657962946819,1125899906842624),(820636957697015,17592186044416),
    (3595328981408141,562949953421312),(3282547830741853,70368744177664),
    (3631502026833111,562949953421312),(1491645815884459,35184372088832))]


def positive_kernel(profile,markers,m):
    terms,_,sources=tables(profile,2);K=[[F(0)]*9 for _ in range(9)]
    for src,dst,gain,z,r,f,e,B,W in terms:
        zeros=STATES[src//3][0]+STATES[src%3][0]+B-W
        value=F(2)**((32-m)*gain-m*zeros)
        for marker,count in zip(markers[:6],(z,r,f,e,B,W)):value*=marker**count
        K[src][dst]+=value
    return K,sources


def exact_rows(profile,markers,m):
    K,sources=positive_kernel(profile,markers,m)
    h=[[F(1),markers[6+2*j],markers[7+2*j]] for j in range(4)];out=[]
    for label in range(4):
        end=h[label].copy()
        if profile[1]==profile[2]==0 and not label&1:end[0]=F(0)
        v=[sum(K[i][j]*end[j//3]*end[j%3] for j in range(9)) for i in range(9)]
        middle=[[v[3*(i//3)+j//3]*v[3*(i%3)+j%3] for j in range(9)] for i in range(9)]
        left=[[sum(K[i][k]*middle[k][j] for k in range(9)) for j in range(9)] for i in range(9)]
        result=[[sum(left[i][k]*K[j][k] for k in range(9)) for j in range(9)] for i in range(9)]
        out.append([result[i][j]/h[0][i//3]/h[1][i%3]/h[2][j//3]/h[3][j%3]
                    for i in sources for j in sources])
    return out


def forced_outer_dimension(profile):
    R,z,_,_,_,B,_,_=profile
    forced=max(0,2*z-B-b)
    # Each RS row has at most R remaining positions. For g<=R,
    # max(0,g-128) <= (R-128)*g/R. No random-route assertion is used.
    assert R>128
    return ((R-128)*(R*M-forced))//R,forced


def exact_replay(profile,markers=CONTROL_MARKERS,m=28,power=512,retain_forced_zeros=False):
    from layered_general_components import rational_operator
    R,z,r,f,e,B,dim,_=profile;P=(R*M-b)//2;V=P-z+f;D=V+r-2*e;KR=R*M-N//2;W=400*N//1000
    if retain_forced_zeros:KR=min(KR,forced_outer_dimension(profile)[0])
    h,_,ends=rational_operator(profile,markers,m,2);blocks=exact_rows(profile,markers,m)
    factors=[]
    for c in range(2):
        first=h[2*c][1]*markers[5]
        if z or r:first+=1
        factors.append((first*(1+h[2*c+1][1]*markers[5]),b//4))
    factors.extend((max(row),(P-b//2)//16) for row in blocks)
    factors.extend((max(row.values()),b//4) for row in ends)
    factors.extend((1/marker,count) for marker,count in zip(markers[:6],(z,r,f,e,B,W)))
    value=power*(1+m*(V+b+dim)+(32-m)*KR-32*(D+e))
    value+=sum(count*log_upper_units(v,power) for v,count in factors)
    return -((-value)//power)


def exact_tensor_test():
    # Independent dense81x81 composition with integer local weights.
    # theta=0, all markers1, and h=(1,2,4). Large integers remain exact.
    from itertools import product
    profile=(192,1,1,1,0,0,0,0);markers=[F(1)]*6+[F(2),F(4)]*4
    K,_=positive_kernel(profile,markers,0);states=list(product(range(3),repeat=4))
    initial=np.array([[K[3*a[0]+a[1]][3*c[0]+c[1]]*K[3*a[2]+a[3]][3*c[2]+c[3]] for c in states] for a in states],dtype=object)
    following=np.array([[K[3*a[0]+a[2]][3*c[0]+c[2]]*K[3*a[1]+a[3]][3*c[1]+c[3]] for c in states] for a in states],dtype=object)
    end=np.array([F(2)**sum(c) for c in states],dtype=object)
    direct=initial@(following@end);row=exact_rows(profile,markers,0)[0]
    assert all(direct[i]/end[i]==row[i] for i in range(81))
    return dict(exact_dense_rows=81,successive_matchings_checked=True)


def logsum2(a,axis=None):return logsumexp(np.asarray(a)*np.log(2),axis=axis)/np.log(2)


def logmat(A,B):return logsum2(A[:,:,None]+B[None,:,:],axis=1)


def rows(x,profile):
    theta=x[-1];h=np.zeros((4,3));h[:,1:]=x[6:14].reshape(4,2)
    terms,last,sources=tables(profile,2)
    K=np.full((9,9),-1e100)
    for src,dst,gain,z,r,f,e,B,W in terms:
        zero=STATES[src//3][0]+STATES[src%3][0]+B-W
        value=32*((1-theta)*gain-theta*zero)-np.dot(x[:6],(z,r,f,e,B,W))
        K[src,dst]=np.logaddexp(K[src,dst]*np.log(2),value*np.log(2))/np.log(2)
    den=(h[0,:,None]+h[1,None,:]).reshape(9,1)+(h[2,:,None]+h[3,None,:]).reshape(1,9)
    out=[]
    for label in range(4):
        end=h[label].copy()
        if profile[1]==profile[2]==0 and not label&1:end[0]=-1e100
        v=logsum2(K+(end[:,None]+end[None,:]).reshape(1,9),axis=1).reshape(3,3)
        # After the first layer, the second matching is (0,2),(1,3).
        middle=(v[:,None,:,None]+v[None,:,None,:]).reshape(9,9)
        value=logmat(logmat(K,middle),K.T)-den
        out.append(value[np.ix_(sources,sources)].ravel())
    ends=[{s:[] for s in sources} for _ in range(2)]
    for src,gain,z,r,f,e in last:
        for c in range(2):
            ends[c][src].append(32*(1-theta)*gain-np.dot(x[:4],(z,r,f,e))-
                                h[2*c,src//3]-h[2*c+1,src%3])
    return np.array(out),np.array([[logsum2(v[s]) for s in sources] for v in ends])


def objective(x,profile):
    R,z,r,f,e,B,dim,_=profile;P=(R*M-b)//2;V=P-z+f;D=V+r-2*e;KR=R*M-N//2
    h=x[6:14].reshape(4,2);w=x[5];theta=x[-1];initial=0
    for c in range(2):
        first=h[2*c,0]-w
        if z or r:first=logsum2([0,first])
        initial+=first+logsum2([0,h[2*c+1,0]-w])
    return (initial*(b//4)+(P-b//2)//16*sum(x[14:18])+b//4*sum(x[18:20])+
            np.dot(x[:6],(z,r,f,e,B,400*N//1000))+
            32*(theta*(V+b+dim)+(1-theta)*KR-D-e)+1)/(R*M)


def propose(profile):
    started=time.monotonic();calls=0;latest=None
    # Near the one-update optimum for the retained substantial-00 profile.
    x=np.array([27.,2.5,-22.5,-5.2,-22.,29.,3.,6.,3.,6.,3.,6.,3.,6.,0.,0.,0.,0.,0.,0.,.875])
    rr,ee=rows(x,profile);x[14:18]=rr.max(axis=1);x[18:20]=ee.max(axis=1)
    def fn(x):
        nonlocal calls,latest
        calls+=1
        if time.monotonic()-started>30:raise TimeoutError('30s two-layer proposal cap')
        latest=x.copy()
        return objective(x,profile)
    def constraints(x):
        rr,ee=rows(x,profile)
        return np.concatenate(((x[14:18,None]-rr).ravel(),(x[18:20,None]-ee).ravel()))
    m=28
    bounds=[(-256,256)]*5+[(0,128)]+[(-128,128)]*8+[(-2048,2048)]*4+[(-512,512)]*2+[(m/32,m/32)]
    try:
        result=minimize(fn,x,method='SLSQP',constraints=[dict(type='ineq',fun=constraints)],bounds=bounds,
                        options=dict(maxiter=120,ftol=1e-10))
        x=result.x;stop=str(result.message)
    except TimeoutError as error:
        x=latest;stop=str(error)
    # A capped/infeasible optimizer iterate is only a marker proposal. Replace
    # every proposed norm by the actual maximum row before reporting a bound.
    rr,ee=rows(x,profile);x[14:18]=rr.max(axis=1);x[18:20]=ee.max(axis=1)
    markers=[F(2.0**(-v if i<6 else v)) for i,v in enumerate(x[:14])]
    return dict(m=m,markers=[(v.numerator,v.denominator) for v in markers],
                floating_log2_upper=objective(x,profile)*profile[0]*M,calls=calls,stop=stop,
                elapsed_seconds=time.monotonic()-started,whole_code_certificate=False)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    R,z,r,f,e,B=192,49136,24568,589632,368520,49136
    dim=mixed_dimension(R,(R*M-b)//2,z,r,f,e,envelopes(R))
    profile=(R,z,r,f,e,B,dim,0)
    if args.verify:
        started=time.monotonic();print(exact_tensor_test())
        print(dict(exact_log2_upper=exact_replay(profile),
                   with_forced_input_zeros=exact_replay(profile,retain_forced_zeros=True),
                   outer_dimension_and_zero_floor=forced_outer_dimension(profile),
                   class_closed=False,whole_code_certificate=False))
        assert time.monotonic()-started<20
    else:print(propose(profile))
