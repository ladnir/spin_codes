"""Read-only independent positive-Arb replay of all three Q1 bounds.

Run the existing BCH verifier separately to reconstruct its outer certificate.
This script uses unnormalized support counts and binary region powering;
the producing script uses normalized support counts and sequential regions.
"""
import math
import sys
from fractions import Fraction as F

from flint import arb,ctx
import bridge as base

sys.path.insert(0,str(base.BCH/'code'))
from audit_bch_q1_full_arb import rational


def coefficients(t,s,spectrum,tenth):
    lam=(arb(tenth)/10).exp(); z=(-lam).exp()
    den=(1<<s)-1; kappa=arb(den)/(den-1); d=min(spectrum)
    m0=sum((arb(c)*z**w for w,c in spectrum.items()),arb(0))/den
    m1=sum((arb(c)*(w*z**(w-1)+(t-w)*z**(w+1)) for w,c in spectrum.items()),arb(0))/(t*den)
    zero=[[arb(1),arb(0),arb(0)],[arb(0),arb(0),z**d],[arb(0),arb(0),kappa*m0]]
    active=[[arb(0),z,arb(0)],[z**(d-1)/den,arb(0),z**(d-1)],[kappa*m1/den,arb(0),kappa*m1]]
    def mul(a,b):
        return [[sum((a[i][k]*b[k][j] for k in range(3)),arb(0)) for j in range(3)] for i in range(3)]
    def add(a,b):
        return [[a[i][j]+b[i][j] for j in range(3)] for i in range(3)]
    def polymul(a,b):
        return [mul(a[0],b[0]),add(mul(a[0],b[1]),mul(a[1],b[0]))]
    identity=[[arb(int(i==j)) for j in range(3)] for i in range(3)]
    empty=[[arb(0) for _ in range(3)] for _ in range(3)]
    exponent=base.ROWS//t; result=[identity,empty]; power=[zero,active]
    while exponent:
        if exponent&1:
            result=polymul(result,power)
        exponent>>=1
        if exponent:
            power=polymul(power,power)
    rz=result[0]; ra=[[v/(base.ROWS//t) for v in row] for row in result[1]]
    current=[[arb(1),arb(0),arb(0)]]
    for n in range(256):
        updated=[]
        for w in range(n+2):
            v=[arb(0)]*3
            if w<=n:
                v=[sum((current[w][k]*rz[k][j] for k in range(3)),arb(0)) for j in range(3)]
            if w:
                v=[v[j]+sum((current[w-1][k]*ra[k][j] for k in range(3)),arb(0)) for j in range(3)]
            updated.append(v)
        current=updated
    correction=base.ROWS*(base.CUTOFF*lam).exp()
    return [sum(v,arb(0))*correction/math.comb(256,w) for w,v in enumerate(current)]


def verify(name):
    t,s,spectrum=base.load_map(name)
    cert=base.read(base.HERE/'generated'/f'{name}_activation_q1_outward.json')
    assert cert['configuration']==name
    assert cert['parameters']==dict(message_bits=1<<20,output_bits=1<<21,outer_length=256,
                                   outer_rows=8192,distance_cutoff=209716,step_bits=t,state_bits=s)
    for filename,digest in cert['source_sha256'].items():
        path=base.BCH/'generated/shift_rank_oa29_joint'/f'{filename[10:]}.json' if filename.startswith('BCH_joint_') else base.HERE/filename
        assert base.sha(path)==digest,filename
    screen=base.read(base.HERE/'generated'/f'{name}_activation_screen.json')
    co={int(w):base.decode(v) for w,v in cert['coefficient_upper'].items()}
    assert set(co)==set(base.WEIGHTS)
    groups={}
    for w in co:
        j=screen['coefficient_rows'][str(w)]['witness_tenth']
        assert isinstance(j,int) and -120<=j<=0
        groups.setdefault(j,[]).append(w)
    ctx.prec=512
    for tenth,weights in groups.items():
        values=coefficients(t,s,spectrum,tenth)
        for w in weights:
            assert 0<min(F(base.ROWS),rational(values[w].upper()))<=co[w],(name,w)
    upper,factor,rest=base.bch_bound(co)
    assert (upper,factor,rest)==tuple(base.decode(cert[k]) for k in ('Q1_upper','dual_domination_factor','remaining_shell_upper'))
    assert upper<F(1,1<<(48 if t==256 else 49))
    assert cert['Q1_below_2_to_minus_40'] is True and cert['all_occupations_certified'] is False
    print(name,'92 independent Arb coefficients and exact BCH aggregation passed',flush=True)


if __name__=='__main__':
    for name in base.CONFIGS:
        verify(name)
