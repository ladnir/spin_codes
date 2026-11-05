"""Bounded joint component/exception diagnostic for the retained r=0 class.

Three states: inactive, active-pinned, active-live. Pair-tag agreement is
relaxed, but both endpoint tag totals and the actual nine-bit wiring remain.
Floating optimization only proposes dual weights; integer max-plus verifies
the resulting uniform kernel cap. No whole-code certificate or native timing.
"""
import time
import argparse
import numpy as np
from scipy.optimize import minimize

STATES=((0,0),(1,0),(1,1))
R=192;M=8192;b=512;T=R*M//b-1
f=589632;e=368520;B=49136;W=838860
COUNTS=np.array([f,f,e,e,B,W],dtype=np.int64)


def branches(role,s):
    # Tags: full-head, full-tail, exception-head, exception-tail.
    # role=0 is the first coordinate of the current pair.
    if role==0:
        if not s:return ()
        return ((-1,(1,0,0,0)),(0,(1,0,1,0)),(1,(0,0,0,0)))
    if s:return ((-1,(0,1,0,0)),(1,(0,1,0,1)))
    return ((-1,(0,0,0,0)),)


def local_terms(role):
    out=[]
    for src,(s,live) in enumerate(STATES):
        for gamma,tags in branches(role,s):
            for sn in (0,1):
                for zero in ((True,False) if s and sn else (s==sn,)):
                    linked=zero or gamma==1
                    pinned=not sn or gamma==0 or (zero and gamma==1)
                    gain=0 if linked else live
                    vn=(live and not pinned) if linked else (not pinned)
                    dst=STATES.index((sn,int(vn)))
                    out.append((src,dst,gain,tags+(int(not s and sn),int(not zero))))
    return out


TERMS=(local_terms(0),local_terms(1))


def coordinate_step(src,gamma):
    s,live=STATES[src]
    for sn in (0,1):
        for zero in ((True,False) if s and sn else (s==sn,)):
            linked=zero or gamma==1
            pinned=not sn or gamma==0 or (zero and gamma==1)
            vn=(live and not pinned) if linked else (not pinned)
            yield STATES.index((sn,int(vn))),0 if linked else live,int(not s and sn),int(not zero)


def paired_terms():
    terms=[];last=[]
    for i,(s,v) in enumerate(STATES):
        for j,(ss,vv) in enumerate(STATES):
            if not s:continue  # hard class has no00 or01 pairs
            choices=((-1,-1,1,0),(0,1,1,1)) if ss else ((1,-1,0,0),)
            for g,gg,full,exc in choices:
                last.append((3*i+j,v+vv+(g==-1)+(gg==-1),full,exc))
                for di,k,bb,ww in coordinate_step(i,g):
                    for dj,kk,bbb,www in coordinate_step(j,gg):
                        terms.append((3*i+j,3*di+dj,k+kk,full,exc,bb+bbb,ww+www))
    return terms,last


PAIR_TERMS,LAST_TERMS=paired_terms()


def logsum2(values):
    return np.logaddexp.reduce(np.asarray(values)*np.log(2))/np.log(2)


def pair_rows(x):
    # x=(full price, exception price, birth price, weight price, h1, h2,
    #    log2 internal norm, log2 terminal norm).
    a,c,d,w,h1,h2=x[:6];h=(0,h1,h2);theta=x[8] if len(x)>8 else 0
    rows={src:[] for src in range(3,9)};last={src:[] for src in range(3,9)}
    for src,dst,gain,full,exc,birth,weight in PAIR_TERMS:
        # Internal zero rank is union(current,next) minus output weight.
        active=STATES[src//3][0]+STATES[src%3][0]
        zeros=active+birth-weight
        rows[src].append(32*((1-theta)*gain-theta*zeros)-a*full-c*exc-d*birth-w*weight+
                         h[dst//3]+h[dst%3]-h[src//3]-h[src%3])
    for src,gain,full,exc in LAST_TERMS:
        last[src].append(32*(1-theta)*gain-a*full-c*exc-h[src//3]-h[src%3])
    return np.array([logsum2(rows[s]) for s in rows]),np.array([logsum2(last[s]) for s in last])


def paired_proposal():
    started=time.monotonic();calls=0;P=(R*M-b)//2
    def objective(x):
        nonlocal calls
        calls+=1
        if time.monotonic()-started>30:raise TimeoutError('30s paired proposal cap')
        a,c,d,w,h1,h2,lam,term,theta=x
        return (b*logsum2([0,h1-w])+(P-b//2)*lam+b//2*term+
                a*f+c*e+d*B+w*W+
                32*(theta*(1375808+b+163923)+(1-theta)*524288-638768-368520)+1)/(R*M)
    def constraints(x):
        rows,last=pair_rows(x)
        return np.concatenate((x[6]-rows,x[7]-last))
    x=np.array([60.,-60.,32.,2.,0.,32.,0.,0.,0.])
    rows,last=pair_rows(x);x[6]=max(rows);x[7]=max(last)
    result=minimize(objective,x,method='SLSQP',constraints=[dict(type='ineq',fun=constraints)],
                    bounds=[(-256,256)]*3+[(0,128)]+[(-128,128)]*2+[(-512,512)]*2+[(0,1)],
                    options=dict(maxiter=200,ftol=1e-10))
    x=result.x;rows,last=pair_rows(x);x[6]=max(rows);x[7]=max(last)
    print(dict(kind='paired_positive_potential',calls=calls,stop=str(result.message),
               weights=x.tolist(),floating_log2_upper=objective(x)*(R*M),
               elapsed_seconds=time.monotonic()-started,whole_code_certificate=False))


def oriented_rows(x):
    a,c,d,w,h10,h20,h11,h21=x[:8];theta=x[11]
    h=((0,h10,h20),(0,h11,h21))
    rows=[{src:[] for src in range(3,9)} for _ in (0,1)]
    last={src:[] for src in range(3,9)}
    for src,dst,gain,full,exc,birth,weight in PAIR_TERMS:
        active=STATES[src//3][0]+STATES[src%3][0]
        zeros=active+birth-weight
        for role in (0,1):
            if role==0 and (dst//3==0 or dst%3==0):continue
            score=32*((1-theta)*gain-theta*zeros)-a*full-c*exc-d*birth-w*weight
            score+=h[role][dst//3]+h[role][dst%3]-h[0][src//3]-h[1][src%3]
            rows[role][src].append(score)
    for src,gain,full,exc in LAST_TERMS:
        last[src].append(32*(1-theta)*gain-a*full-c*exc-h[0][src//3]-h[1][src%3])
    return (np.array([[logsum2(rows[r][s]) for s in range(3,9)] for r in (0,1)]),
            np.array([logsum2(last[s]) for s in range(3,9)]))


def oriented_proposal():
    started=time.monotonic();P=(R*M-b)//2;calls=0
    def objective(x):
        nonlocal calls
        calls+=1
        if time.monotonic()-started>30:raise TimeoutError('30s oriented proposal cap')
        a,c,d,w,h10,h20,h11,h21,l0,l1,term,theta=x
        initial=b//2*(h10-w+logsum2([0,h11-w]))
        return (initial+(P-b//2)//2*(l0+l1)+b//2*term+a*f+c*e+d*B+w*W+
                32*(theta*(1375808+b+163923)+(1-theta)*524288-638768-368520)+1)/(R*M)
    def constraints(x):
        rows,last=oriented_rows(x)
        return np.concatenate((x[8]-rows[0],x[9]-rows[1],x[10]-last))
    x=np.array([-19.54,-6.87,-19.84,28.23,2.85,6.80,2.85,6.80,0.,0.,0.,.835])
    rows,last=oriented_rows(x);x[8:10]=rows.max(axis=1);x[10]=max(last)
    result=minimize(objective,x,method='SLSQP',constraints=[dict(type='ineq',fun=constraints)],
                    bounds=[(-256,256)]*3+[(0,128)]+[(-128,128)]*4+[(-512,512)]*3+[(0,1)],
                    options=dict(maxiter=250,ftol=1e-10))
    x=result.x;rows,last=oriented_rows(x);x[8:10]=rows.max(axis=1);x[10]=max(last)
    print(dict(kind='oriented_paired_positive_potential',calls=calls,stop=str(result.message),
               weights=x.tolist(),floating_log2_upper=objective(x)*(R*M),
               elapsed_seconds=time.monotonic()-started,whole_code_certificate=False))


def lookahead_rows(x,depth):
    size=1<<depth;contexts=size//2
    a,c,d,w=x[:4];theta=x[-1]
    h=np.zeros((size,3));h[:,1:]=x[4:4+2*size].reshape(size,2)
    rows=[[[] for _ in range(6)] for _ in range(size)]
    last=[[[] for _ in range(6)] for _ in range(contexts)]
    for src,dst,gain,full,exc,birth,weight in PAIR_TERMS:
        active=STATES[src//3][0]+STATES[src%3][0]
        zeros=active+birth-weight
        base=32*((1-theta)*gain-theta*zeros)-a*full-c*exc-d*birth-w*weight
        for z in range(size):
            if not z&1 and (dst//3==0 or dst%3==0):continue
            context=z&(contexts-1)
            score=base+h[z,dst//3]+h[z,dst%3]-h[2*context,src//3]-h[2*context+1,src%3]
            rows[z][src-3].append(score)
    for src,gain,full,exc in LAST_TERMS:
        for context in range(contexts):
            last[context][src-3].append(32*(1-theta)*gain-a*full-c*exc-
                                     h[2*context,src//3]-h[2*context+1,src%3])
    return np.array([[logsum2(v) for v in row] for row in rows]),np.array([[logsum2(v) for v in row] for row in last])


def lookahead_proposal(depth):
    assert 1<=depth<=3
    started=time.monotonic();P=(R*M-b)//2;calls=0;size=1<<depth;contexts=size//2
    def objective(x):
        nonlocal calls
        calls+=1
        if time.monotonic()-started>45:raise TimeoutError('45s lookahead proposal cap')
        a,c,d,w=x[:4];theta=x[-1];h=x[4:4+2*size].reshape(size,2)
        initial=sum(h[2*j,0]-w+logsum2([0,h[2*j+1,0]-w]) for j in range(contexts))*(b//size)
        value=(initial+(P-b//2)//size*sum(x[4+2*size:4+3*size])+
               b//size*sum(x[4+3*size:-1])+a*f+c*e+d*B+w*W+
               32*(theta*(1375808+b+163923)+(1-theta)*524288-638768-368520)+1)
        return value/(R*M)
    def constraints(x):
        rows,last=lookahead_rows(x,depth)
        return np.concatenate(((x[4+2*size:4+3*size,None]-rows).ravel(),
                               (x[4+3*size:-1,None]-last).ravel()))
    hs=[v for _ in range(contexts) for v in (2.68,6.50,1.96,5.72)]
    x=np.array([-19.84,-6.72,-21.04,28.35]+hs+[0.]*(size+contexts)+[.8415])
    rows,last=lookahead_rows(x,depth);x[4+2*size:4+3*size]=rows.max(axis=1);x[4+3*size:-1]=last.max(axis=1)
    result=minimize(objective,x,method='SLSQP',constraints=[dict(type='ineq',fun=constraints)],
                    bounds=[(-256,256)]*3+[(0,128)]+[(-128,128)]*(2*size)+[(-512,512)]*(size+contexts)+[(27/32,27/32)],
                    options=dict(maxiter=200,ftol=1e-10))
    x=result.x;rows,last=lookahead_rows(x,depth);x[4+2*size:4+3*size]=rows.max(axis=1);x[4+3*size:-1]=last.max(axis=1)
    print(dict(kind='lookahead_paired_positive_potential',depth=depth,calls=calls,stop=str(result.message),
               weights=x.tolist(),floating_log2_upper=objective(x)*(R*M),
               elapsed_seconds=time.monotonic()-started,whole_code_certificate=False))


def product(A,C,temperature=None):
    vals=[A[:,:,k,None]+C[:,None,k,:] for k in range(3)]
    if temperature is None:return np.maximum.reduce(vals)
    return temperature*np.logaddexp.reduce(np.array(vals)/temperature,axis=0)


def evaluate(weights,temperature=None,unit=1,theta=0):
    """Return dual upper value in kernel-dimension units (integer if unit>1)."""
    integer=np.issubdtype(np.asarray(weights).dtype,np.integer) and temperature is None and theta==0
    neg=-(10**14) if integer else -1e100
    mats=[]
    for role in (0,1):
        mat=np.full((3,3),neg,dtype=np.int64 if integer else float)
        for src,dst,gain,counts in TERMS[role]:
            zero=int(STATES[src][0] and counts[5]==0)
            value=unit*((1-theta)*gain-theta*zero)-np.dot(weights,counts)
            if temperature is None:mat[src,dst]=max(mat[src,dst],value)
            else:mat[src,dst]=temperature*np.logaddexp(mat[src,dst]/temperature,value/temperature)
        mats.append(mat)
    identity=np.full((b,3,3),neg,dtype=mats[0].dtype)
    for j in range(3):identity[:,j,j]=0
    labels=np.arange(b)
    phases=[np.stack(mats)[(labels>>(t%9))&1] for t in range(1,10)]
    cycle=identity.copy()
    for mat in phases:cycle=product(cycle,mat,temperature)
    periods,tail=divmod(T-1,9)
    value=identity.copy();power=cycle
    while periods:
        if periods&1:value=product(value,power,temperature)
        periods//=2
        if periods:power=product(power,power,temperature)
    for mat in phases[:tail]:value=product(value,mat,temperature)
    initial=np.array([0,-weights[5],neg],dtype=mats[0].dtype)
    # Q1 is pinned by X0=Q1. Last output constraints are omitted; each free
    # terminal input introduces one additional component. Fixed inputs do not.
    terminal=np.full((2,3),neg,dtype=mats[0].dtype)
    for role in (0,1):
        for src,(s,live) in enumerate(STATES):
            for gamma,tags in branches(role,s):
                score=unit*(1-theta)*(live+(gamma==-1))-np.dot(weights[:4],tags)
                if temperature is None:terminal[role,src]=max(terminal[role,src],score)
                else:terminal[role,src]=temperature*np.logaddexp(terminal[role,src]/temperature,score/temperature)
    end=terminal[(labels>>(T%9))&1]
    values=initial[None,:,None]+value+end[:,None,:]
    if temperature is None:total=np.max(values,axis=(1,2)).sum()
    else:total=temperature*np.logaddexp.reduce((values/temperature).reshape(b,-1),axis=1).sum()
    return total+np.dot(weights,COUNTS)


def main(blend=False):
    started=time.monotonic();calls=0;best=[float('inf'),None]
    def objective(x):
        nonlocal calls
        if time.monotonic()-started>40:raise TimeoutError('40s proposal cap')
        calls+=1;theta=x[6] if blend else 0
        v=evaluate(x[:6],temperature=1/(32*np.log(2)),theta=theta)
        if blend:v+=theta*(1375808+b+163923-524288)
        if v<best[0]:best[:]=[v,x.copy()]
        return v/(R*M)
    try:
        start=[.9,1.,-.91,-.93,1.06,.07]+([0.] if blend else [])
        result=minimize(objective,np.array(start),method='L-BFGS-B',
                        bounds=[(-4,4)]*5+[(0,4)]+([(0,1)] if blend else []),
                        options=dict(maxfun=2200,maxiter=250,eps=1e-6,ftol=1e-13,gtol=1e-8))
        reason=result.message
    except TimeoutError as error:reason=str(error)
    x=best[1];denom=4096;integer=np.rint(x*denom).astype(np.int64)
    cap=None if blend else int(evaluate(integer,unit=denom))//denom
    log_bound=32*(best[0]+524288-638768-368520)+1
    print(dict(proposal_calls=calls,stop=str(reason),weights=x.tolist(),
               integer_weights=integer.tolist(),denominator=denom,
               verified_relaxed_kernel_upper=cap,target_kernel_upper=407155,
               floating_joint_log2_upper=log_bound,
               elapsed_seconds=time.monotonic()-started,whole_code_certificate=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--paired',action='store_true');parser.add_argument('--blend',action='store_true');parser.add_argument('--oriented',action='store_true')
    parser.add_argument('--lookahead',type=int,default=0)
    args=parser.parse_args()
    if args.lookahead:lookahead_proposal(args.lookahead)
    elif args.oriented:oriented_proposal()
    elif args.paired:paired_proposal()
    else:main(args.blend)
