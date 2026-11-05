"""Bounded all-pair-type component counts. Research only; no code changes.

The first replay reproduces the previous selected screen's forty failures.
Selected profiles and proposed numerical bounds are not whole-code coverage.
"""
import time
import argparse
from itertools import product
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from layered_mixed_gate import (M,b,N,envelopes,mixed_dimension,no_empty_point,
                               occupation_markers,occupation_floor)
from check_layered_component_count import MARKERS,rational_rows,log_upper_units
from layered_component_transfer import STATES,coordinate_step,logsum2


def all_terms():
    terms=[];last=[]
    for i,(s,v) in enumerate(STATES):
        for j,(ss,vv) in enumerate(STATES):
            if s and ss:choices=((-1,-1,1,0),(0,1,1,1))
            elif s:choices=((1,-1,0,0),)
            elif ss:choices=((-1,-1,0,0),)
            else:choices=((0,0,0,0),)
            z=int(not s and not ss);r=int(not s and ss)
            for g,gg,full,exc in choices:
                last.append((3*i+j,v+vv+(g==-1)+(gg==-1),z,r,full,exc))
                for di,k,bb,ww in coordinate_step(i,g):
                    for dj,kk,bbb,www in coordinate_step(j,gg):
                        terms.append((3*i+j,3*di+dj,k+kk,z,r,full,exc,bb+bbb,ww+www))
    return terms,last


ALL_TERMS,ALL_LAST=all_terms()


def tables(profile,depth):
    _,z,r,_,_,_,_,_=profile
    terms=[v for v in ALL_TERMS if (z or not v[3]) and (r or not v[4])]
    last=[v for v in ALL_LAST if (z or not v[2]) and (r or not v[3])]
    sources=sorted({v[0] for v in last})
    return terms,last,sources


def general_rows(x,profile,depth):
    size=1<<depth;contexts=size//2;theta=x[-1]
    h=np.zeros((size,3));h[:,1:]=x[6:6+2*size].reshape(size,2)
    terms,last,sources=tables(profile,depth)
    rows=[{s:[] for s in sources} for _ in range(size)]
    ends=[{s:[] for s in sources} for _ in range(contexts)]
    forced=profile[1]==profile[2]==0
    for src,dst,gain,z,r,f,e,B,W in terms:
        zeros=STATES[src//3][0]+STATES[src%3][0]+B-W
        base=32*((1-theta)*gain-theta*zeros)-np.dot(x[:6],(z,r,f,e,B,W))
        for label in range(size):
            if forced and not label&1 and (dst//3==0 or dst%3==0):continue
            context=label&(contexts-1)
            rows[label][src].append(base+h[label,dst//3]+h[label,dst%3]-h[2*context,src//3]-h[2*context+1,src%3])
    for src,gain,z,r,f,e in last:
        for context in range(contexts):
            ends[context][src].append(32*(1-theta)*gain-np.dot(x[:4],(z,r,f,e))-
                                      h[2*context,src//3]-h[2*context+1,src%3])
    return (np.array([[logsum2(row[s]) for s in sources] for row in rows]),
            np.array([[logsum2(row[s]) for s in sources] for row in ends]))


def general_objective(x,profile,depth):
    R,z,r,f,e,B,dim,_=profile;P=(R*M-b)//2;V=P-z+f;D=V+r-2*e;KR=max(0,R*M-N//2)
    size=1<<depth;contexts=size//2;W=400*N//1000;theta=x[-1]
    h=x[6:6+2*size].reshape(size,2);w=x[5];initial=0
    for c in range(contexts):
        first=h[2*c,0]-w
        if z or r:first=logsum2([0,first])
        initial+=first+logsum2([0,h[2*c+1,0]-w])
    value=(initial*(b//size)+(P-b//2)//size*sum(x[6+2*size:6+3*size])+
           b//size*sum(x[6+3*size:-1])+np.dot(x[:6],(z,r,f,e,B,W))+
           32*(theta*(V+b+dim)+(1-theta)*KR-D-e)+1)
    return value/(R*M)


def propose(profile,depth=2):
    started=time.monotonic();size=1<<depth;contexts=size//2;calls=0
    def objective(x):
        nonlocal calls
        calls+=1
        if time.monotonic()-started>45:raise TimeoutError('45s general-profile proposal cap')
        return general_objective(x,profile,depth)
    def constraints(x):
        rows,last=general_rows(x,profile,depth)
        return np.concatenate(((x[6+2*size:6+3*size,None]-rows).ravel(),
                               (x[6+3*size:-1,None]-last).ravel()))
    hs=[v for _ in range(contexts) for v in (2.68,6.50,1.96,5.72)]
    x=np.array([20.,20.,-19.84,-6.72,-21.04,28.35]+hs+[0.]*(size+contexts)+[.84375])
    rows,last=general_rows(x,profile,depth);x[6+2*size:6+3*size]=rows.max(axis=1);x[6+3*size:-1]=last.max(axis=1)
    bounds=[(-256,256)]*5+[(0,128)]+[(-128,128)]*(2*size)+[(-512,512)]*(size+contexts)+[(0,1)]
    result=minimize(objective,x,method='SLSQP',constraints=[dict(type='ineq',fun=constraints)],
                    bounds=bounds,options=dict(maxiter=180,ftol=1e-10))
    m=max(0,min(32,round(32*result.x[-1])));bounds[-1]=(m/32,m/32)
    result=minimize(objective,result.x,method='SLSQP',constraints=[dict(type='ineq',fun=constraints)],
                    bounds=bounds,options=dict(maxiter=120,ftol=1e-10))
    x=result.x;rows,last=general_rows(x,profile,depth);x[6+2*size:6+3*size]=rows.max(axis=1);x[6+3*size:-1]=last.max(axis=1)
    markers=[F(2.0**(-v if i<6 else v)) for i,v in enumerate(x[:6+2*size])]
    return dict(depth=depth,m=m,markers=[(v.numerator,v.denominator) for v in markers],
                floating_log2_upper=objective(x)*profile[0]*M,calls=calls,
                elapsed_seconds=time.monotonic()-started)


def rational_operator(profile,markers,m,depth=2):
    size=1<<depth;contexts=size//2;terms,last,sources=tables(profile,depth)
    h=[[F(1),markers[6+2*c],markers[7+2*c]] for c in range(size)]
    rows=[{s:F(0) for s in sources} for _ in range(size)]
    ends=[{s:F(0) for s in sources} for _ in range(contexts)]
    forced=profile[1]==profile[2]==0
    for src,dst,gain,z,r,f,e,B,W in terms:
        zeros=STATES[src//3][0]+STATES[src%3][0]+B-W
        base=F(2)**((32-m)*gain-m*zeros)
        for marker,count in zip(markers[:6],(z,r,f,e,B,W)):base*=marker**count
        for label in range(size):
            if forced and not label&1 and (dst//3==0 or dst%3==0):continue
            c=label&(contexts-1)
            rows[label][src]+=base*h[label][dst//3]*h[label][dst%3]/h[2*c][src//3]/h[2*c+1][src%3]
    for src,gain,z,r,f,e in last:
        base=F(2)**((32-m)*gain)
        for marker,count in zip(markers[:4],(z,r,f,e)):base*=marker**count
        for c in range(contexts):ends[c][src]+=base/h[2*c][src//3]/h[2*c+1][src%3]
    return h,rows,ends


def general_exact(profile,markers,m,depth=2,power=512):
    size=1<<depth;contexts=size//2
    h,rows,ends=rational_operator(profile,markers,m,depth)
    forced=profile[1]==profile[2]==0
    R,z,r,f,e,B,dim,_=profile;P=(R*M-b)//2;V=P-z+f;D=V+r-2*e;KR=max(0,R*M-N//2);W=400*N//1000
    factors=[]
    for c in range(contexts):
        first=h[2*c][1]*markers[5]
        if not forced:first+=1
        factors.append((first*(1+h[2*c+1][1]*markers[5]),b//size))
    factors.extend((max(row.values()),(P-b//2)//size) for row in rows)
    factors.extend((max(row.values()),b//size) for row in ends)
    factors.extend((1/marker,count) for marker,count in zip(markers[:6],(z,r,f,e,B,W)))
    upper=power*(1+m*(V+b+dim)+(32-m)*KR-32*(D+e))
    upper+=sum(count*log_upper_units(value,power) for value,count in factors)
    return -((-upper)//power)


def retained_failures():
    started=time.monotonic();failed=[];considered=0;impossible=0
    for R in (129,160,168,176,192,208,224,240,256):
        P=(R*M-b)//2;env=envelopes(R);markers=occupation_markers(R)
        for f in sorted({P*j//16 for j in range(17)}):
            for r in sorted({(P-f)*j//8 for j in range(9)}):
                for e in sorted({f*j//8 for j in range(9)}):
                    W=400*N//1000;dim=mixed_dimension(R,P,0,r,f,e,env)
                    for B in sorted({min(P-f,W//2)*j//4 for j in range(5)}):
                        if min(r,P-f-r)<occupation_floor(R,f,B,markers):
                            impossible+=1;continue
                        old=no_empty_point(R,f,r,e,env,births=B,dimension=dim)
                        considered+=1
                        if old['best']>-100:failed.append((R,0,r,f,e,B,dim,old['best']))
                    if time.monotonic()-started>40:raise TimeoutError('40s retained-screen cap')
    assert len(failed)==40 and considered==45967 and impossible==9194
    return failed


def fixed_markers_replay(profile,power=512):
    R,z,r,f,e,B,dim,_=profile
    if z or r:return None
    P=(R*M-b)//2;V=P+f;D=P+f-2*e;KR=R*M-N//2;W=400*N//1000
    h,rows,last=rational_rows(MARKERS,2)
    factors=[(h[2*c][1]*MARKERS[3]*(1+h[2*c+1][1]*MARKERS[3]),b//4) for c in range(2)]
    factors.extend((max(row),(P-b//2)//4) for row in rows)
    factors.extend((max(row),b//4) for row in last)
    factors.extend((1/m,count) for m,count in zip(MARKERS[:4],(f,e,B,W)))
    upper=power*(1+27*(V+b+dim)+5*KR-32*(D+e))
    upper+=sum(count*log_upper_units(value,power) for value,count in factors)
    return -((-upper)//power)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',type=int,default=-1)
    parser.add_argument('--bank',action='store_true')
    parser.add_argument('--depth',type=int,default=2,choices=(1,2,3))
    parser.add_argument('--empty-count',type=int,default=0)
    parser.add_argument('--empty-screen',action='store_true')
    args=parser.parse_args()
    failed=retained_failures()
    if args.empty_screen:
        from layered_component_markers import BANK,depth_of
        _,m,markers=BANK[1];R,_,r,f,e,B,_,_=failed[35];P=(R*M-b)//2;env=envelopes(R)
        values=[]
        for z in (0,1,2048,4096,8192,16384,32768,49136,98304):
            dim=mixed_dimension(R,P,z,r,f,e,env);p=(R,z,r,f,e,B,dim,0)
            values.append((z,dim,general_exact(p,markers,m,depth_of(markers))))
        print(dict(selected_empty_profiles=values,whole_code_certificate=False))
        return
    if args.bank:
        from layered_component_markers import BANK,depth_of
        values=[]
        for i,p in enumerate(failed):
            best=min((general_exact(p,markers,m,depth=depth_of(markers)),name) for name,m,markers in BANK)
            values.append((i,p,best))
        print(dict(profiles=40,passed=sum(v[2][0]<-100 for v in values),
                   worst=max(values,key=lambda v:v[2][0]),
                   remaining=[v for v in values if v[2][0]>=-100],whole_code_certificate=False))
        return
    if args.case>=0:
        profile=failed[args.case]
        if args.empty_count:
            R,_,r,f,e,B,_,old=profile;z=args.empty_count;P=(R*M-b)//2
            assert 0<z<=P-r-f
            dim=mixed_dimension(R,P,z,r,f,e,envelopes(R))
            profile=(R,z,r,f,e,B,dim,old)
        proposal=propose(profile,args.depth)
        exact=general_exact(profile,[F(*v) for v in proposal['markers']],proposal['m'],depth=args.depth)
        print(dict(case=args.case,profile=profile,proposal=proposal,exact_log2_upper=exact,whole_code_certificate=False))
        return
    values=[(p,fixed_markers_replay(p)) for p in failed]
    print(dict(retained_failures=len(failed),r_zero=sum(p[2]==0 for p in failed),
               fixed_markers_passed=sum(v is not None and v<-100 for _,v in values),
               profiles=values,whole_code_certificate=False))


if __name__=='__main__':main()
