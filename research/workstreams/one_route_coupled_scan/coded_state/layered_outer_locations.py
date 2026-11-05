"""Joint fixed-input locations, temporal components, and active RS rows.

Research bound for the unchanged layered construction. Float proposals are
diagnostics; --exact replays fixed rational inputs. No encoder benchmarking.
"""
import math
import time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, gammaln
from layered_general_components import M, b, N, STATES, coordinate_step, envelopes, mixed_dimension
from layered_two_layer_components import CONTROL_MARKERS
from check_layered_component_count import log_upper_units
from layered_interpolation import exact_tau

LOCATION_MARKERS=[F(n,d) for n,d in (
    (3516193804967843,144115188075855872),(2582595331189461,288230376151711744),
    (3836966332952327,8796093022208),(8434367173395669,549755813888),
    (1571897855966257,137438953472),(2856970697196541,75557863725914323419136),
    (589583120677623,281474976710656),(2829102902432735,8796093022208),
    (589583120677173,281474976710656),(707275725607825,2199023255552),
    (1179166241354685,562949953421312),(2829102902431767,8796093022208),
    (4716664965417853,2251799813685248),(2829102902431349,8796093022208))]


def location_terms():
    out=[];last=[]
    for i,(si,vi) in enumerate(STATES):
        for j,(sj,vj) in enumerate(STATES):
            choices=((-1,-1,1,0),(0,1,1,1)) if si and sj else (
                ((1,-1,0,0),) if si else (((-1,-1,0,0),) if sj else ((0,0,0,0),)))
            z=int(not si and not sj);r=int(not si and sj)
            for g,h,f,e in choices:
                last.append((3*i+j,vi+vj+(g==-1)+(h==-1),z,r,f,e))
                for di,gi,bi,wi in coordinate_step(i,g):
                    for dj,gj,bj,wj in coordinate_step(j,h):
                        def classify(src,dst,gamma,w):
                            # Tags classify ACTUAL input values, not the auxiliary kernel.
                            if gamma<0:return (0,0,0)
                            s=STATES[src][0];t=STATES[dst][0]
                            if gamma==0:return (1-t,t,0)
                            # gamma=1: X=Qnext+Q. At Y=0, delta != 1.
                            return (int(not s and not t),int(s!=t or (s and t and not w)),int(s and t and w))
                        ji,ui,ci=classify(i,di,g,wi);jj,uj,cj=classify(j,dj,h,wj)
                        out.append((3*i+j,3*di+dj,gi+gj,z,r,f,e,bi+bj,wi+wj,ji+jj,ui+uj,ci+cj))
    return out,last


TERMS,LAST=location_terms()


def hardest():
    R,z,r,f,e,B=192,49136,24568,589632,368520,49136
    d=mixed_dimension(R,(R*M-b)//2,z,r,f,e,envelopes(R))
    return (R,z,r,f,e,B,d,0)


def positive_rows(markers,m,alpha,two=False,split=True):
    tau=exact_tau(m);q=F(1<<32)
    h=[[F(1),markers[6+2*j],markers[7+2*j]] for j in range(4)]
    zero=1-alpha+alpha*tau/q
    K=[[F(0)]*9 for _ in range(9)]
    for src,dst,gain,z,r,f,e,B,W,J,U,C in TERMS:
        zeros=STATES[src//3][0]+STATES[src%3][0]+B-W
        val=q**gain/tau**(gain+zeros)*zero**J*alpha**U
        if split:val*=(alpha+zero/tau)**C
        for marker,count in zip(markers[:6],(z,r,f,e,B,W)):val*=marker**count
        K[src][dst]+=val
    rows=[]
    for label in range(4):
        c=label&1
        rows.append([sum(K[src][dst]*h[label][dst//3]*h[label][dst%3] for dst in range(9))
                     /h[2*c][src//3]/h[2*c+1][src%3] for src in range(9)])
    if two:
        rows=[]
        for label in range(4):
            v=[sum(K[i][j]*h[label][j//3]*h[label][j%3] for j in range(9)) for i in range(9)]
            middle=[[v[3*(i//3)+j//3]*v[3*(i%3)+j%3] for j in range(9)] for i in range(9)]
            left=[[sum(K[i][k]*middle[k][j] for k in range(9)) for j in range(9)] for i in range(9)]
            result=[[sum(left[i][k]*K[j][k] for k in range(9)) for j in range(9)] for i in range(9)]
            rows.append([result[i][j]/h[0][i//3]/h[1][i%3]/h[2][j//3]/h[3][j%3]
                         for i in range(9) for j in range(9)])
    ends=[]
    for c in range(2):
        vals=[F(0)]*9
        for src,gain,z,r,f,e in LAST:
            val=(q/tau)**gain
            for marker,count in zip(markers[:4],(z,r,f,e)):val*=marker**count
            vals[src]+=val/h[2*c][src//3]/h[2*c+1][src%3]
        ends.append(vals)
    return h,rows,ends,zero


def exact(profile,s,markers=CONTROL_MARKERS,m=28,power=512,two=False,split=True,colour=None):
    # Historical integer-exponent replay. General rational tau weights use
    # the directed coefficient replay in cover_layered_location_regions.
    assert isinstance(m,int)
    R,z,r,f,e,B,dbar,_=profile;P=(R*M-b)//2;V=P-z+f;D=V+r-2*e
    assert 1<=s<=M and 0<=m<=32 and 0<markers[5]<=1
    alpha=F(s,M) if colour is None else colour
    assert 0<alpha<=1 and (alpha<1 or s==M)
    h,rows,ends,zero=positive_rows(markers,m,alpha,two,split)
    # d_S <= dbar and <= dim(C_S); deliberately retain the old route envelope.
    ds=min(dbar,(R-128)*s);factors=[]
    for c in range(2):
        factors.append(((zero+alpha*h[2*c][1]*markers[5])*
                        (zero+alpha*h[2*c+1][1]*markers[5]),b//4))
    factors.extend((max(row),(P-b//2)//(16 if two else 4)) for row in rows)
    factors.extend((max(row),b//4) for row in ends)
    factors.extend((1/marker,count) for marker,count in zip(markers[:6],(z,r,f,e,B,400*N//1000)))
    # Any Bernoulli colour probability gives the SAME uniform s-subset law
    # after conditioning. Keep its exact normalization, not a heuristic tilt.
    choose=F(math.comb(M,s))
    factors.extend(((choose,1),(1/choose,R),(1/alpha,R*s)))
    if s<M:factors.append((1/(1-alpha),R*(M-s)))
    value=power*(1+m*(V+b+ds)+(32-m)*(R-128)*s-32*(D+e))
    value+=sum(count*log_upper_units(v,power) for v,count in factors)
    return -((-value)//power)


def logsum2(a,axis=None):return logsumexp(np.asarray(a)*np.log(2),axis=axis)/np.log(2)


def rows(x,s,m,two=False,split=True):
    m=float(m)
    alpha=1/(1+2.**(-x[20])) if len(x)>20 else s/M
    lj=math.log2(1-alpha+alpha*2**(-(32-m)));lu=math.log2(alpha)
    h=np.zeros((4,3));h[:,1:]=x[6:14].reshape(4,2)
    K=np.full((9,9),-1e100)
    lc=math.log2(alpha+2**(lj-m)) if split else 0.
    for src,dst,gain,z,r,f,e,B,W,J,U,C in TERMS:
        zeros=STATES[src//3][0]+STATES[src%3][0]+B-W
        val=(32-m)*gain-m*zeros-np.dot(x[:6],(z,r,f,e,B,W))+J*lj+U*lu+C*lc
        K[src,dst]=np.logaddexp(K[src,dst]*np.log(2),val*np.log(2))/np.log(2)
    rr=[]
    for label in range(4):
        c=label&1
        rr.append(logsum2(K+(h[label,:,None]+h[label,None,:]).reshape(1,9),axis=1)
                  -(h[2*c,:,None]+h[2*c+1,None,:]).reshape(9))
    if two:
        from layered_two_layer_components import logmat
        rr=[]
        den=(h[0,:,None]+h[1,None,:]).reshape(9,1)+(h[2,:,None]+h[3,None,:]).reshape(1,9)
        for label in range(4):
            v=logsum2(K+(h[label,:,None]+h[label,None,:]).reshape(1,9),axis=1).reshape(3,3)
            middle=(v[:,None,:,None]+v[None,:,None,:]).reshape(9,9)
            rr.append((logmat(logmat(K,middle),K.T)-den).ravel())
    ee=[]
    for c in range(2):
        vals=[[] for _ in range(9)]
        for src,gain,z,r,f,e in LAST:
            vals[src].append((32-m)*gain-np.dot(x[:4],(z,r,f,e))-
                             h[2*c,src//3]-h[2*c+1,src%3])
        ee.append([logsum2(v) for v in vals])
    return np.array(rr),np.array(ee)


def objective(x,profile,s,m,two=False):
    m=float(m)
    R,z,r,f,e,B,dbar,_=profile;P=(R*M-b)//2;V=P-z+f;D=V+r-2*e
    alpha=1/(1+2.**(-x[20])) if len(x)>20 else s/M
    lj=math.log2(1-alpha+alpha*2**(-(32-m)));lu=math.log2(alpha)
    choose=(gammaln(M+1)-gammaln(s+1)-gammaln(M-s+1))/np.log(2)
    conditioning=-choose-s*lu-((M-s)*math.log2(1-alpha) if s<M else 0)
    h=x[6:14].reshape(4,2)
    initial=sum(logsum2([lj,lu+h[c,0]-x[5]]) for c in range(4))*(b//4)
    return (initial+(P-b//2)//(16 if two else 4)*sum(x[14:18])+b//4*sum(x[18:20])+
            np.dot(x[:6],(z,r,f,e,B,400*N//1000))+
            1+m*(V+b+min(dbar,(R-128)*s))+(32-m)*(R-128)*s-32*(D+e)+
            choose+R*conditioning)/(R*M)


def start(markers=CONTROL_MARKERS):
    return np.array([(-1 if i<6 else 1)*math.log2(v) for i,v in enumerate(markers)]+[0.]*6)


def evaluate(profile,s,m,x,two=False):
    x=x.copy();rr,ee=rows(x,s,m,two);x[14:18]=rr.max(axis=1);x[18:20]=ee.max(axis=1)
    return objective(x,profile,s,m,two)*profile[0]*M


def propose(profile,s,m,seconds=25,two=False):
    began=time.monotonic();latest=start(LOCATION_MARKERS if m<=24 else CONTROL_MARKERS);calls=0
    if s<M:latest=np.append(latest,math.log2(s/(M-s)))
    rr,ee=rows(latest,s,m,two);latest[14:18]=rr.max(axis=1);latest[18:20]=ee.max(axis=1)
    def fn(x):
        nonlocal calls,latest
        if time.monotonic()-began>seconds:raise TimeoutError('bounded location proposal')
        calls+=1;latest=x.copy()
        return objective(x,profile,s,m,two)
    def cons(x):
        rr,ee=rows(x,s,m,two)
        return np.concatenate(((x[14:18,None]-rr).ravel(),(x[18:20,None]-ee).ravel()))
    try:
        res=minimize(fn,latest,method='SLSQP',constraints=[dict(type='ineq',fun=cons)],
                     bounds=[(-256,256)]*5+[(0,128)]+[(-128,128)]*8+[(-2048,2048)]*6+([(-16,24)] if s<M else []),
                     options=dict(maxiter=140,ftol=1e-10))
        latest=res.x;stop=str(res.message)
    except TimeoutError as error:stop=str(error)
    markers=[F(2.**(-v if i<6 else v)) for i,v in enumerate(latest[:14])]
    colour=F(1/(1+2.**(-latest[20]))) if s<M else F(1)
    return markers,dict(s=s,m=m,two=two,float_upper=evaluate(profile,s,m,latest,two),
                        exact_upper=exact(profile,s,markers,m,two=two,colour=colour),
                        colour=(colour.numerator,colour.denominator),calls=calls,
                        seconds=time.monotonic()-began,stop=stop)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--s',type=int);p.add_argument('--m',type=int,default=28)
    p.add_argument('--location-screen',action='store_true')
    p.add_argument('--two',action='store_true')
    p.add_argument('--json',action='store_true')
    p.add_argument('--offsets',nargs=5,type=int,default=(0,0,0,0,0),metavar=('Z','R','F','E','B'))
    p.add_argument('--exact',action='store_true');args=p.parse_args();profile=hardest()
    if any(args.offsets):
        counts=tuple(v+d for v,d in zip(profile[1:6],args.offsets));R=profile[0]
        dim=mixed_dimension(R,(R*M-b)//2,*counts[:4],envelopes(R))
        profile=(R,)+counts+(dim,0)
    if args.location_screen:
        x=start(LOCATION_MARKERS)
        print(dict(selected_s=[(s,round(evaluate(profile,s,23,x),2))
                   for s in [1]+list(range(256,M+1,256))],whole_code_certificate=False))
    elif args.s:
        if args.exact:print(dict(s=args.s,m=args.m,two=args.two,exact_upper=exact(profile,args.s,m=args.m,two=args.two)))
        else:
            markers,result=propose(profile,args.s,args.m,two=args.two)
            if args.json:
                import json
                # Preserve integers through clients whose JSON numbers are doubles.
                print(json.dumps(dict(result=result,markers=[(str(v.numerator),str(v.denominator)) for v in markers])))
            else:
                print(result);print([(v.numerator,v.denominator) for v in markers])
    else:
        print(dict(profile=profile,selected_s=[(s,round(evaluate(profile,s,args.m,start()),2))
                   for s in (1,512,1024,2048,2394,3072,4096,5120,6144,7168,7680,8192)],
                   whole_code_certificate=False))
