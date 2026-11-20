"""Q1 bridge with an explicit deterministic activation state.

States: zero, arbitrary nonzero, and live with density <= kappa times
uniform nonzero. See README.md for the weighted-measure invariant.
No import or write in the source worktree. Outputs are write-once.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import bridge as base


def matrices(t, s, spectrum, z, number):
    den = (1 << s)-1
    kappa = number(den)/(den-1)
    m0 = sum((number(c)*z**w/den for w,c in spectrum.items()), number(0))
    m1 = sum((number(c)*(w*z**(w-1)+(t-w)*z**(w+1))/(t*den)
              for w,c in spectrum.items()), number(0))
    d = min(spectrum)
    one, zero = number(1), number(0)
    return ((one,zero,zero, zero,zero,z**d, zero,zero,kappa*m0),
            (zero,z,zero, z**(d-1)/den,zero,z**(d-1), kappa*m1/den,zero,kappa*m1))


def logmul(a,b):
    # Three fixed intermediary states, no dynamic matrix kernel dispatch.
    return np.logaddexp(np.logaddexp(a[:,0,None]+b[0,None,:],
                                    a[:,1,None]+b[1,None,:]),
                                    a[:,2,None]+b[2,None,:])


def log_regions(t,s,spectrum,lam):
    den=(1<<s)-1
    m0=base.lse([math.log(c/den)-lam*w for w,c in spectrum.items()])
    m1=base.lse([math.log(c/den*w/t)-lam*(w-1) for w,c in spectrum.items()]
                +[math.log(c/den*(t-w)/t)-lam*(w+1) for w,c in spectrum.items() if w<t])
    k=math.log1p(1/(den-1))
    d=min(spectrum)
    zero=np.full((3,3),-np.inf)
    zero[0,0]=0.; zero[1,2]=-lam*d; zero[2,2]=k+m0
    active=np.full((3,3),-np.inf)
    active[0,1]=-lam
    active[1,0]=-lam*(d-1)-math.log(den); active[1,2]=-lam*(d-1)
    active[2,0]=k+m1-math.log(den); active[2,2]=k+m1
    rz=np.full((3,3),-np.inf); np.fill_diagonal(rz,0.)
    ra=np.full((3,3),-np.inf)
    for _ in range(base.ROWS//t):
        ra=np.logaddexp(logmul(ra,zero),logmul(rz,active))
        rz=logmul(rz,zero)
    return rz,ra-math.log(base.ROWS//t)


def log_vector_mul(v,m):
    return np.logaddexp(np.logaddexp(v[:,0,None]+m[0,None,:],
                                    v[:,1,None]+m[1,None,:]),
                                    v[:,2,None]+m[2,None,:])


def log_coefficients(rz,ra,length=base.LENGTH):
    current=np.full((length+1,3),-np.inf); current[0,0]=0.
    for n in range(length):
        updated=np.full_like(current,-np.inf)
        updated[:n+1]=log_vector_mul(current[:n+1],rz)
        updated[1:n+2]=np.logaddexp(updated[1:n+2],log_vector_mul(current[:n+1],ra))
        current=updated
    return np.logaddexp.reduce(current,axis=1)-np.array([math.log(math.comb(length,w)) for w in range(length+1)])


def screen(name):
    t,s,spectrum=base.load_map(name)
    best=np.full(base.LENGTH+1,np.inf); witness=np.zeros(base.LENGTH+1,dtype=int)
    for tenth in range(-120,1):
        lam=math.exp(tenth/10)
        value=log_coefficients(*log_regions(t,s,spectrum,lam))+base.CUTOFF*lam
        value=math.log(base.ROWS)+np.minimum(0.,value)
        witness[value<best]=tenth; best=np.minimum(best,value)
    co={w:F.from_float(math.exp(max(-700.,float(best[w])))) for w in base.WEIGHTS}
    upper,factor,rest=base.bch_bound(co)
    payload=dict(status='THREE_STATE_BINARY64_Q1_SCREEN_ONLY',configuration=name,
        coefficient_rows={str(w):dict(log_coefficient=float(best[w]),witness_tenth=int(witness[w])) for w in base.WEIGHTS},
        diagnostic_margin_bits=math.log2(upper.denominator)-math.log2(upper.numerator),
        dual_domination_factor=float(factor),all_occupations_certified=False,
        source_sha256={'activation_bridge.py':base.sha(Path(__file__)),
                       'bridge.py':base.sha(Path(base.__file__)),
                       'inputs/manifest.json':base.sha(base.HERE/'inputs/manifest.json')})
    base.write_new(base.HERE/'generated'/f'{name}_activation_screen.json',payload)
    print(name,'activation-aware Q1 screen margin',payload['diagnostic_margin_bits'],flush=True)


def positive_mul(a,b):
    return tuple(a[3*i]*b[j]+a[3*i+1]*b[3+j]+a[3*i+2]*b[6+j]
                 for i in range(3) for j in range(3))


def positive_regions(t,s,spectrum,z,number,epochs):
    zero,active=matrices(t,s,spectrum,z,number)
    rz=tuple(number(int(i==j)) for i in range(3) for j in range(3))
    ra=(number(0),)*9
    for _ in range(epochs):
        ra=tuple(x+y for x,y in zip(positive_mul(ra,zero),positive_mul(rz,active)))
        rz=positive_mul(rz,zero)
    return rz,tuple(x/epochs for x in ra)


def positive_vector_mul(v,m):
    x,y,z=v
    return (x*m[0]+y*m[3]+z*m[6], x*m[1]+y*m[4]+z*m[7], x*m[2]+y*m[5]+z*m[8])


def positive_coefficients(rz,ra,number,length):
    current=[(number(1),number(0),number(0))]
    for n in range(length):
        updated=[]
        for w in range(n+2):
            a=b=c=number(0)
            if w<=n:
                x,y,z=positive_vector_mul(current[w],rz)
                a+=x*(n+1-w)/(n+1); b+=y*(n+1-w)/(n+1); c+=z*(n+1-w)/(n+1)
            if w:
                x,y,z=positive_vector_mul(current[w-1],ra)
                a+=x*w/(n+1); b+=y*w/(n+1); c+=z*w/(n+1)
            updated.append((a,b,c))
        current=updated
    return [sum(v,number(0)) for v in current]


def certify(name):
    from flint import arb,ctx
    import sys
    sys.path.insert(0,str(base.BCH/'code'))
    from audit_bch_q1_full_arb import rational
    t,s,spectrum=base.load_map(name)
    path=base.HERE/'generated'/f'{name}_activation_screen.json'
    screen_data=base.read(path)
    for filename,digest in screen_data['source_sha256'].items():
        assert base.sha(base.HERE/filename)==digest
    groups={}
    for w in base.WEIGHTS:
        j=screen_data['coefficient_rows'][str(w)]['witness_tenth']
        assert isinstance(j,int) and -120<=j<=0
        groups.setdefault(j,[]).append(w)
    ctx.prec=256; co={}
    for j,weights in sorted(groups.items()):
        lam=(arb(j)/10).exp(); z=(-lam).exp()
        values=positive_coefficients(*positive_regions(t,s,spectrum,z,arb,base.ROWS//t),arb,base.LENGTH)
        correction=base.ROWS*(base.CUTOFF*lam).exp()
        for w in weights:
            co[w]=min(F(base.ROWS),rational((values[w]*correction).upper()))
            assert co[w]>0
        print(name,'Arb witness',j,'weights',weights,flush=True)
    upper,factor,rest=base.bch_bound(co)
    sources={**screen_data['source_sha256'],str(path.relative_to(base.HERE)):base.sha(path)}
    sources.update({f'BCH_joint_{f}':base.sha(base.BCH/'generated/shift_rank_oa29_joint'/f'{f}.json') for f in ('audit','objective')})
    payload=dict(status='ACTIVATION_AWARE_OUTWARD_Q1_BOUND',configuration=name,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_length=256,outer_rows=8192,
                        distance_cutoff=base.CUTOFF,step_bits=t,state_bits=s),
        coefficient_upper={str(w):base.encode(v) for w,v in co.items()},Q1_upper=base.encode(upper),
        dual_domination_factor=base.encode(factor),remaining_shell_upper=base.encode(rest),
        margin_bits_diagnostic=math.log2(upper.denominator)-math.log2(upper.numerator),
        Q1_below_2_to_minus_40=upper<F(1,1<<40),all_occupations_certified=False,
        assumption_scope='Three-state weighted-measure envelope; existing exact BCH outer certificate; fresh independent nonzero field multipliers each epoch',
        source_sha256=sources)
    base.write_new(base.HERE/'generated'/f'{name}_activation_q1_outward.json',payload)
    print(name,'OUTWARD Q1 margin',payload['margin_bits_diagnostic'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=('screen','certify'))
    parser.add_argument('--configuration',required=True,choices=base.CONFIGS)
    args=parser.parse_args()
    (screen if args.mode=='screen' else certify)(args.configuration)
