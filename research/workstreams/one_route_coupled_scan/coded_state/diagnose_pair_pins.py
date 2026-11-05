"""Rejected-uniform SL2 research model: exact pin identity and bounded probes.

Not the old two-shear ensemble. No encoder or frozen verifier is changed.
Only compact stdout; no output files. One numerical process at a time.
"""
import argparse
import itertools
import math
import time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, gammaln
from diagnose_uniform_pairs import field

M, b, N = 8192, 512, 1 << 21

# Frozen proposals. Decimal literals round-trip the original binary64 values.
# They select positive rational proof inputs; the replay does not trust floats.
PIN_BANK = (
    (28,(7273803141626959,9525602955312207),
     (-4.221337159826098,-0.39052653993327585,3.496463508483695,5.394776627977944e-14,
      2.580572043186524,2.580572043186525,2.5805720431865065,2.5805720431865384,1.6916313045010383)),
    (0,(5991779214278069,24006177723760053),
     (-14.452094575928886,-13.178453446578766,-10.817202869757294,1.4901180776885245e-08,
      2.5837376166538277,2.5837376163509913,2.583737616394343,2.5837376160914993,-1.5880941124948043)),
    (0,(7163984599985015,7163984602082167),
     (-61.89291067500755,-0.09629909146985126,4.597357075306449,3.015864135877516e-13,
      2.5843064376543072,2.5843064351812224,2.5843064361131964,2.584306434605766,31.66968365935721)),
    (0,(5023537322463411,14030736577204403),
     (-14.31987050550407,-12.595284294905357,-10.236274915556006,1.4901185898573537e-08,
      2.5817327170515108,2.581732717255463,2.5817327170515045,2.5817327170054174,-0.8423749813544846)),
)
SECOND_SLICE = (-1.8094818746982524,2.314019437883824,0.8043434283053125,
                7.650916116406272e-15,2.045191535769004,2.045191535769156,
                2.0451915357690904,2.045191535768939,1.8681069398912338)


def terms():
    out = []
    for src in range(4):
        a, c = divmod(src, 2)
        for dst in range(4):
            d, e = divmod(dst, 2)
            for w0 in ((0, 1) if a and d else (int(a != d),)):
                for w1 in ((0, 1) if c and e else (int(c != e),)):
                    birth = (not a and d) + (not c and e)
                    weight = w0 + w1
                    zeros = a + c + birth - weight
                    empty = int(src == 0)
                    # Each birth out of an empty pair pins a NEW component.
                    pins = empty * (d + e)
                    fixed_zero = empty * (2 - d - e)
                    out.append((src, dst, empty, int(src == 3), birth,
                                weight, zeros, pins, fixed_zero))
    return np.array(out, dtype=np.int64)


TERMS = terms()


def log2sum(x, axis=None):
    return logsumexp(np.asarray(x) * math.log(2), axis=axis) / math.log(2)


def logs(x, s, m, no_singleton=False, two=False):
    # Log inverse markers z,f,B,W; four label potentials; log2 colour odds.
    prices, odds = x[:4], x[8]
    lp = -np.logaddexp(0., -odds * math.log(2)) / math.log(2)
    lnotp = -np.logaddexp(0., odds * math.log(2)) / math.log(2)
    if s == M:
        lp, lnotp = 0., -np.inf
    lbeta = log2sum([lnotp, lp - (32 - m)])
    src, dst, z, f, birth, weight, zeros, pins, J = TERMS.T
    values = (-TERMS[:, 2:6] @ prices - 32 * zeros
              + pins * (lp - (32 - m)) + J * lbeta)
    K = np.array([[log2sum(values[(src == i) & (dst == j)])
                   for j in range(4)] for i in range(4)])
    if no_singleton:K[1:3,:]=-np.inf
    internal, terminal = [], []
    for label in range(4):
        c = label & 1
        pot = x[4 + label]
        dest = np.array([0., pot, pot, 2 * pot])
        inp = np.array([0., x[4 + 2*c + 1], x[4 + 2*c],
                        x[4 + 2*c] + x[4 + 2*c + 1]])
        internal.append(log2sum(K + dest, axis=1) - inp)
    if two:
        def mat(A,B):return log2sum(A[:,:,None]+B[None,:,:],axis=1)
        den=(np.array([0.,x[5],x[4],x[4]+x[5]])[:,None]
             +np.array([0.,x[7],x[6],x[6]+x[7]])[None,:])
        internal=[]
        for label in range(4):
            pot=x[4+label];dest=np.array([0.,pot,pot,2*pot])
            v=log2sum(K+dest,axis=1).reshape(2,2)
            middle=(v[:,None,:,None]+v[None,:,None,:]).reshape(4,4)
            internal.append((mat(mat(K,middle),K.T)-den).ravel())
    for c in range(2):
        inp = np.array([0., x[4+2*c+1], x[4+2*c], x[4+2*c]+x[4+2*c+1]])
        terminal.append(np.array([-prices[0], 0., 0., -prices[1]]) - inp)
        if no_singleton:terminal[-1][1:3]=-np.inf
    initial = sum(log2sum([lbeta, lp-(32-m)-prices[3]+v]) for v in x[4:8])
    choose = (gammaln(M+1)-gammaln(s+1)-gammaln(M-s+1))/math.log(2)
    conditioning = -choose - s*lp - ((M-s)*lnotp if s<M else 0)
    return initial, np.array(internal), np.array(terminal), choose, conditioning


def value(x, profile, s, m, W, two=False):
    R, z, f, B, dbar = profile
    P = (R*M-b)//2
    V, D = P-z+f, 2*(P-z)
    initial, rows, ends, choose, cond = logs(x, s, m,z+f==P,two)
    return (b//4*initial + (P-b//2)//(16 if two else 4)*rows.max(axis=1).sum()
            + b//4*ends.max(axis=1).sum() + np.dot(x[:4], (z,f,B,W))
            + 1 + 32*(V+b-D) + m*min(dbar,(R-128)*s)
            + (32-m)*(R-128)*s + choose + R*cond)


def exact(x, profile, s, m, W, power=1<<24, colour=None, coefficients=False,two=False):
    """Directed integer replay; m must be integer, others rounded to dyadics."""
    from rational_log_bounds import log2_upper_scaled as lup
    assert isinstance(m, int) and 0<=m<=32
    markers = [F(2.**(-float(v))) for v in x[:4]]
    assert all(v>0 for v in markers) and markers[3]<=1
    h = [F(2.**float(v)) for v in x[4:8]]
    odds = F(2.**float(x[8])); p = odds/(1+odds) if s<M else F(1)
    if colour is not None:p=colour
    beta = 1-p+p/F(1 << (32-m))
    pin = p/F(1 << (32-m)); q = F(1 << 32)
    K = [[F(0) for _ in range(4)] for _ in range(4)]
    for src,dst,z,f,B,ww,Z,H,J in TERMS.tolist():
        v = q**(-Z)*pin**H*beta**J
        for mark,count in zip(markers,(z,f,B,ww)): v *= mark**count
        K[src][dst] += v
    R,z,f,B,dbar = profile; P=(R*M-b)//2; V=P-z+f; D=2*(P-z)
    if z+f==P:K[1]=[F(0)]*4;K[2]=[F(0)]*4
    factors = [(beta+pin*markers[3]*v,b//4) for v in h]
    for label in range(4):
        c=label&1; dest=(F(1),h[label],h[label],h[label]**2)
        inp=(F(1),h[2*c+1],h[2*c],h[2*c]*h[2*c+1])
        if not two:
            rows=[sum(K[i][j]*dest[j] for j in range(4))/inp[i] for i in range(4)]
        else:
            v=[sum(K[i][j]*dest[j] for j in range(4)) for i in range(4)]
            mid=[[v[2*(i//2)+j//2]*v[2*(i%2)+j%2] for j in range(4)] for i in range(4)]
            left=[[sum(K[i][k]*mid[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
            hp=(F(1),h[1],h[0],h[0]*h[1]);hq=(F(1),h[3],h[2],h[2]*h[3])
            rows=[sum(left[i][k]*K[j][k] for k in range(4))/hp[i]/hq[j]
                  for i in range(4) for j in range(4)]
        factors.append((max(rows),(P-b//2)//(16 if two else 4)))
    for c in range(2):
        inp=(F(1),h[2*c+1],h[2*c],h[2*c]*h[2*c+1])
        end=(markers[0],F(int(z+f!=P)),F(int(z+f!=P)),markers[1])
        factors.append((max(v/hh for v,hh in zip(end,inp)),b//4))
    factors.extend((1/v,c) for v,c in zip(markers,(z,f,B,W)))
    base=(1+32*(V+b-D))*power+sum(n*lup(v,power) for v,n in factors)
    if p<1:
        base+=R*M*lup(1/(1-p),power)
        slope=(32-m)*(R-128)*power+R*lup((1-p)/p,power)
        if coefficients:return base,slope,m*power
        units=base+slope*s+m*power*min(dbar,(R-128)*s)+(R-1)*lup(F(1,math.comb(M,s)),power)
    else:
        assert s==M and not coefficients
        units=base+(m*min(dbar,(R-128)*s)+(32-m)*(R-128)*s)*power
    return -((-units)//power)


def propose(profile,s,m,W,seconds=8,two=False):
    began=time.monotonic(); best=[float('inf'),None]; calls=0
    def objective(x):
        nonlocal calls
        if time.monotonic()-began>seconds: raise TimeoutError('bounded pair-pin probe')
        calls+=1
        v=value(x,profile,s,m,W,two)
        if v<best[0]: best[:]=[v,x.copy()]
        return v/(profile[0]*M)
    # Smooth epigraph optimization avoids nonsmooth max-row proposal failures.
    x=np.array([0.,0.,0.,16.,0.,0.,0.,0.,math.log2(s/(M-s)) if s<M else 0.])
    def full_objective(y):
        objective(y[:9])
        R,z,f,B,dbar=profile;P=(R*M-b)//2;V=P-z+f;D=2*(P-z)
        initial,_,_,choose,cond=logs(y[:9],s,m,z+f==P,two)
        return (b//4*initial+(P-b//2)//(16 if two else 4)*sum(y[9:13])+b//4*sum(y[13:15])
                +np.dot(y[:4],(z,f,B,W))+1+32*(V+b-D)+m*min(dbar,(R-128)*s)
                +(32-m)*(R-128)*s+choose+R*cond)/(R*M)
    def cons(y):
        P=(profile[0]*M-b)//2
        _,rows,ends,_,_=logs(y[:9],s,m,profile[1]+profile[2]==P,two)
        vals=np.concatenate(((y[9:13,None]-rows).ravel(),(y[13:15,None]-ends).ravel()))
        return vals[np.isfinite(vals)]
    P=(profile[0]*M-b)//2
    _,rows,ends,_,_=logs(x,s,m,profile[1]+profile[2]==P,two)
    y=np.r_[x,rows.max(axis=1),ends.max(axis=1)]
    try:
        result=minimize(full_objective,y,method='SLSQP',constraints=[dict(type='ineq',fun=cons)],
            bounds=[(-96,96)]*3+[(0,64)]+[(-64,64)]*4+[(-40,40)]+[(-512,512)]*6,
            options=dict(maxiter=120,ftol=1e-9))
        objective(result.x[:9]); reason=str(result.message)
    except TimeoutError as err: reason=str(err)
    replay=exact(best[1],profile,s,m,W,two=two)
    return dict(s=s,m=m,two=two,float_bound=best[0],directed_slice_bound=replay,
                weights=best[1].tolist(),calls=calls,seconds=time.monotonic()-began,stop=reason)


def rank(matrix, table):
    if not matrix:return 0
    a=[list(row) for row in matrix]; pivot=0
    for col in range(len(a[0])):
        j=next((j for j in range(pivot,len(a)) if a[j][col]),None)
        if j is None:continue
        a[pivot],a[j]=a[j],a[pivot]
        inv=next(v for v in range(1,len(table)) if table[a[pivot][col]][v]==1)
        a[pivot]=[table[inv][v] for v in a[pivot]]
        for j in range(len(a)):
            if j!=pivot and a[j][col]:
                c=a[j][col];a[j]=[v^table[c][w] for v,w in zip(a[j],a[pivot])]
        pivot+=1
        if pivot==len(a):break
    return pivot


def test():
    start=time.monotonic();t=field(4,7);checked=0
    # Exhaust all three-vertex support/zero-output histories on four coordinates.
    # Zero constraints form paths; fixed inputs occur initially and after 00 pairs.
    width,vertices=4,3;nv=width*vertices
    for masks in itertools.product(range(1<<width),repeat=vertices):
        supports=[[(v>>j)&1 for j in range(width)] for v in masks]
        possible=[(i,j) for i in range(vertices-1) for j in range(width)
                  if supports[i][j] and supports[i+1][j]]
        for zz in range(1<<len(possible)):
            zeros=[v for k,v in enumerate(possible) if zz>>k&1]
            equations=[]
            for i in range(vertices):
                for j in range(width):
                    if not supports[i][j]:
                        row=[0]*nv;row[i*width+j]=1;equations.append(row)
            for i,j in zeros:
                row=[0]*nv;row[i*width+j]=2;row[(i+1)*width+j]=1;equations.append(row)
            pins=[]
            for j in range(width):
                row=[0]*nv;row[j]=1;pins.append(row)
            H=sum(supports[0])
            for i in range(vertices-1):
                bit=1<<(i%2)
                for j in range(width):
                    if not supports[i][j] and not supports[i][j^bit]:
                        row=[0]*nv;row[(i+1)*width+j]=1;pins.append(row)
                        H+=supports[i+1][j]
            actual=nv-rank(equations,t)
            kernel=nv-rank(equations+pins,t)
            assert actual==sum(map(sum,supports))-len(zeros)
            assert actual-kernel==H,(masks,zeros,actual,kernel,H)
            checked+=1
            if time.monotonic()-start>25:raise TimeoutError('25s exhaustive pin test')
    print(dict(exact_rank_identities=checked,seconds=time.monotonic()-start),flush=True)


def geometry():
    """Exact toy check of the no-birth/no-singleton structural lemma."""
    receipts=[]
    for d in (1,2,3,4):
        width=1<<d;full=(1<<width)-1
        closed=[]
        for phase in range(d):
            bit=1<<phase;pairs=[(j,j^bit) for j in range(width) if not j&bit]
            closed.append([sum((1<<j)|(1<<k) for n,(j,k) in enumerate(pairs) if mask>>n&1)
                           for mask in range(1<<len(pairs))])
        starts=transitions=0
        for start in closed[0]:
            if start==full:continue
            reachable={start};starts+=1
            for phase in range(1,d):
                nxt=set()
                for src in reachable:
                    for dst in closed[phase]:
                        if dst&~src==0:nxt.add(dst);transitions+=1
                reachable=nxt
            assert reachable=={0},(d,start,reachable)
        receipts.append(dict(width=width,proper_starts=starts,transitions=transitions))
    print(dict(no_birth_no_singleton_exact_checks=receipts,
               K20_statement='A proper support cannot stay nonempty across nine distinct pairing phases.',
               whole_code_certificate=False),flush=True)


def tensor_test():
    """Independent literal four-coordinate sum for the two-layer contraction."""
    checks=0
    for exclude in (False,True):
        K=[[F(1+4*i+j,17) if not (exclude and i in (1,2)) else F(0)
            for j in range(4)] for i in range(4)]
        for h in (F(1,3),F(2)):
            v=[sum(K[i][j]*h**(j.bit_count()) for j in range(4)) for i in range(4)]
            mid=[[v[2*(i//2)+j//2]*v[2*(i%2)+j%2] for j in range(4)] for i in range(4)]
            for start in range(16):
                i,j=divmod(start,4)
                contracted=sum(K[i][a]*mid[a][c]*K[j][c] for a in range(4) for c in range(4))
                direct=F(0)
                for middle in range(16):
                    a,c=divmod(middle,4);ac=2*(a//2)+c//2;bd=2*(a%2)+c%2
                    for end in range(16):
                        e,f=divmod(end,4);ef=2*(e//2)+f//2;gh=2*(e%2)+f%2
                        direct+=K[i][a]*K[j][c]*K[ac][ef]*K[bd][gh]*h**end.bit_count()
                assert direct==contracted
                checks+=1
    print(dict(exact_two_layer_contractions=checks),flush=True)


def cover(profile,W,seconds=90):
    """Directed enumeration of all active-row counts for ONE aggregate profile."""
    from diagnose_layered_families import source_binomials
    from cover_layered_location_regions import POWER
    began=time.monotonic();svec=np.arange(1,M+1,dtype=np.int64)
    R,_,_,_,dbar=profile;bn=source_binomials()
    best=np.full(M,2**62,dtype=np.int64);bank=[]
    seeds=[(4096,28),(1024,0),(6144,0)]
    for i in range(14):
        if time.monotonic()-began>seconds:break
        if i<len(seeds):s,m=seeds[i]
        else:s,m=int(best.argmax())+1,(0,16,28,32)[(i-len(seeds))%4]
        proposal=propose(profile,s,m,W,seconds=min(6,max(1,seconds-(time.monotonic()-began))))
        x=np.array(proposal['weights']);odds=F(2.**float(x[8]));p=odds/(1+odds)
        base,slope,rate=exact(x,profile,s,m,W,colour=p,coefficients=True)
        vals=base+slope*svec+rate*np.minimum(dbar,(R-128)*svec)+(R-1)*bn
        best=np.minimum(best,vals)
        bank.append(dict(m=m,colour=(p.numerator,p.denominator),weights=x.tolist()))
        worst=-((-int(best.max()))//POWER)
        print(dict(proposal=i,s=s,m=m,per_s_upper=worst,worst_s=int(best.argmax())+1,
                   failed_source_sizes=int((best>=-100*POWER).sum())),flush=True)
        if worst<-100:break
    print(dict(profile=profile,W=W,per_profile_good_route_upper=-((-int(best.max()))//POWER)+13,
               all_source_sizes_checked=M,all_source_sizes_pass=bool(np.all(best<-100*POWER)),
               bank=bank,seconds=time.monotonic()-began,whole_code_certificate=False),flush=True)


def replay(distance=400):
    from rational_log_bounds import log2_upper_scaled as lup
    from cover_layered_location_regions import region_env
    began=time.monotonic();R,z,f,B=192,393088,294816,98272
    P=(R*M-b)//2;D=2*(P-z);_,E,u=region_env(R,'factorial')
    dbar=min(D//2,(R-128)*M,max(0,int((np.minimum(D,2*E)-128*u).max())))
    assert dbar==64000
    profile=R,z,f,B,dbar;W=distance*N//1000;svec=np.arange(1,M+1,dtype=np.int64)
    out=[]
    for power in (1<<24,1<<28):
        bn=np.array([lup(F(1,math.comb(M,s)),power) for s in range(1,M+1)],dtype=np.int64)
        best=np.full(M,2**62,dtype=np.int64)
        for m,colour,x in PIN_BANK:
            base,slope,rate=exact(x,profile,1,m,W,power=power,colour=F(*colour),coefficients=True)
            # Python integers independently check the vector arithmetic limits.
            assert abs(base)+abs(slope)*M+abs(rate)*dbar+(R-1)*int(max(abs(bn)))<2**62
            vals=base+slope*svec+rate*np.minimum(dbar,(R-128)*svec)+(R-1)*bn
            best=np.minimum(best,vals)
        bound=-((-int(best.max()))//power)+13
        assert bound<-100
        out.append(dict(log_scale=power,good_route_class_upper=bound,worst_s=int(best.argmax())+1))
        if time.monotonic()-began>45:raise TimeoutError('45s selected-profile replay cap')
    print(dict(profile=profile,distance_permille=distance,all_source_counts=M,replays=out,
               seconds=time.monotonic()-began,whole_code_certificate=False),flush=True)


def replay_second():
    from cover_layered_location_regions import region_env
    R,z,f,B=160,327552,327552,81888;P=(R*M-b)//2;D=2*(P-z)
    _,E,u=region_env(R,'factorial')
    dbar=min(D//2,(R-128)*M,max(0,int((np.minimum(D,2*E)-128*u).max())))
    assert dbar==0
    bounds=[exact(SECOND_SLICE,(R,z,f,B,dbar),4096,27,400*N//1000,power=p,two=True)
            for p in (1<<24,1<<28)]
    assert all(v<-100 for v in bounds)
    print(dict(second_profile_source_slice=4096,directed_bounds=bounds,
               whole_profile_closed=False,whole_code_certificate=False),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--test',action='store_true')
    p.add_argument('--geometry',action='store_true')
    p.add_argument('--replay',action='store_true')
    p.add_argument('--two',action='store_true')
    p.add_argument('--tensor-test',action='store_true')
    p.add_argument('--replay-second',action='store_true')
    p.add_argument('--case',type=int,default=0);p.add_argument('--s',type=int,default=M)
    p.add_argument('--m',type=int,default=0);p.add_argument('--seconds',type=float,default=8)
    p.add_argument('--distance',type=int,default=400);p.add_argument('--cover',action='store_true');args=p.parse_args()
    assert 1<=args.s<=M and 0<=args.m<=32 and 0<args.seconds<=30
    if args.test:return test()
    if args.tensor_test:return tensor_test()
    if args.geometry:return geometry()
    if args.replay:return replay(args.distance)
    if args.replay_second:return replay_second()
    from cover_layered_location_regions import region_env
    R,z,f,B=((192,393088,294816,98272),(160,327552,327552,81888),
             (256,524160,524160,0),(200,49136,589632,49136))[args.case]
    P=(R*M-b)//2;D=2*(P-z);_,E,u=region_env(R,'factorial')
    dbar=min(D//2,(R-128)*M,max(0,int((np.minimum(D,2*E)-128*u).max())))
    profile=(R,z,f,B,dbar)
    if args.cover:return cover(profile,args.distance*N//1000)
    result=propose(profile,args.s,args.m,args.distance*N//1000,args.seconds,two=args.two)
    print(dict(profile=profile,distance_permille=args.distance,result=result,
               whole_profile_closed=False,whole_code_certificate=False),flush=True)


if __name__=='__main__':main()
