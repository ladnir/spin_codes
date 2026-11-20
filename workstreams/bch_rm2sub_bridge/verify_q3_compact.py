"""Read-only independent 512-bit replay of every Q3 conditioning box."""
import itertools
import math
from fractions import Fraction as F

from flint import arb,ctx
import bridge as base
import occupation_three as q3


def independent_regions(t,s,spectrum,z,length=8192):
    den=(1<<s)-1;kappa=arb(den)/(den-1);d=min(spectrum)
    epoch=[]
    for j in range(4):
        moment=arb(0)
        for w,count in spectrum.items():
            for v in range(max(0,j-(t-w)),min(w,j)+1):
                moment+=arb(count*math.comb(w,v)*math.comb(t-w,j-v))*z**(w+j-2*v)/(den*math.comb(t,j))
        if j==0:
            matrix=[[arb(1),arb(0),arb(0)],[arb(0),arb(0),z**d],[arb(0),arb(0),kappa*moment]]
        else:
            matrix=[[arb(0),z**j,arb(0)],[z**(d-j)/den,arb(0),z**(d-j)],
                    [kappa*moment/den,arb(0),kappa*moment]]
        epoch.append([[v*math.comb(t,j) for v in row] for row in matrix])
    def mul(a,b):
        return [[sum((a[i][k]*b[k][j] for k in range(3)),arb(0)) for j in range(3)] for i in range(3)]
    def polymul(a,b):
        out=[]
        for degree in range(4):
            terms=[mul(a[j],b[degree-j]) for j in range(degree+1)]
            out.append([[sum((m[i][k] for m in terms),arb(0)) for k in range(3)] for i in range(3)])
        return out
    identity=[[arb(int(i==j)) for j in range(3)] for i in range(3)]
    empty=[[arb(0) for _ in range(3)] for _ in range(3)]
    result=[identity,empty,empty,empty];power=epoch;exponent=length//t
    while exponent:
        if exponent&1:result=polymul(result,power)
        exponent>>=1
        if exponent:power=polymul(power,power)
    return [[[v/math.comb(length,j) for v in row] for row in result[j]] for j in range(4)]


def independent_moment(region,ps):
    matrix=[[arb(0) for _ in range(3)] for _ in range(3)]
    for bits in itertools.product((0,1),repeat=3):
        probability=math.prod(p if bit else 1-p for p,bit in zip(ps,bits))
        factor=arb(probability.numerator)/probability.denominator
        for i in range(3):
            for j in range(3):matrix[i][j]+=factor*region[sum(bits)][i][j]
    state=[arb(1),arb(0),arb(0)]
    for _ in range(256):
        state=[sum((state[k]*matrix[k][j] for k in range(3)),arb(0)) for j in range(3)]
    return sum(state,arb(0))


def verify():
    import sys
    sys.path.insert(0,str(base.BCH/'code'))
    from audit_bch_q1_full_arb import rational
    path=base.HERE/'generated'/'t128_s15_q3_compact_outward.json';receipt=base.read(path)
    for filename,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/filename)==digest
    for filename,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/filename)==digest
    assert receipt['parameters']==dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,
                                       step_bits=128,state_bits=15,occupation=3,distance_cutoff=209716)
    t,s,spectrum=base.load_map('t128_s15');caps=q3.q2.deterministic_caps();ctx.prec=512
    bands=receipt['bands'];assert sorted(w for band in bands for w in band)==list(base.WEIGHTS)
    expected=list(itertools.combinations_with_replacement(range(len(bands)),3))
    assert [tuple(row['bands']) for row in receipt['boxes']]==expected
    cache={};total=F(0)
    for row in receipt['boxes']:
        box=row['bands'];ps=[base.decode(v) for v in row['p']];tenth=row['witness_tenth']
        assert len(ps)==3 and -120<=tenth<=0
        assert all(0<p<=1 and (p<1 or bands[i]==[256]) for i,p in zip(box,ps))
        if tenth not in cache:
            lam=(arb(tenth)/10).exp()
            cache[tenth]=(independent_regions(t,s,spectrum,(-lam).exp()),(209716*lam).exp())
        region,correction=cache[tenth]
        recomputed=independent_moment(region,ps)*correction
        moment=base.decode(row['moment_upper'])
        assert 0<rational(recomputed.upper())<=moment
        costs=[base.decode(v) for v in row['conditioning_upper']]
        for band,p,cost in zip((bands[i] for i in box),ps,costs):
            exact=sum((F(caps[w],math.comb(256,w))/(p**w*(1-p)**(256-w)) for w in band),F(0))
            assert exact<=cost<=exact*(1+F(1,1<<190))
        multiplicity=len(set(itertools.permutations(box)))
        term=math.comb(8192,3)*multiplicity*math.prod(costs)*moment
        assert term==base.decode(row['contribution_upper']);total+=term
    assert total==base.decode(receipt['Q3_upper'])<F(1,1<<126)
    first=base.read(base.HERE/'generated'/'t128_s15_activation_q1_outward.json')
    second=base.read(base.HERE/'generated'/'t128_s15_q2_j-75.json')
    partial=total+base.decode(first['Q1_upper'])+base.decode(second['Q2_upper'])
    assert partial==base.decode(receipt['Q123_upper'])<F(1,1<<49)
    assert receipt['all_occupations_certified'] is False
    print('All 455 boxes independently replayed at 512 bits; exact conditioning and Q123 aggregation passed',flush=True)


if __name__=='__main__':verify()
