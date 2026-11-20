"""Outward OA29 bounds on pairs of BCH tail words with small difference."""
import argparse
import math
from pathlib import Path
from fractions import Fraction as F
from flint import arb,ctx
import bridge as base
import christoffel_caps as christoffel
import screen_exponential_modes as cap_source


def run(verify=False):
    output=base.HERE/'generated/bch_tail_overlap_outward.json';old=base.read(output) if verify else None
    if not verify:assert not output.exists()
    if old:
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        caps=christoffel.deterministic_caps();sources=[base.HERE/name for name in old['cap_receipts']]
        for path in sources:
            receipt=base.read(path)
            for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
            for name,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
            assert receipt['rational_primal_dual_checks_passed']
            w=receipt['weight'];caps[w]=caps[256-w]=min(caps[w],receipt['cap'])
    else:caps,sources=cap_source.latest_caps()
    christoffel.build(verify=True);ctx.prec=512 if verify else 256
    tail_path=base.HERE/'generated/joint_tail_lower_80/lower.json';tail=base.read(tail_path)
    for name,digest in tail['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    for name,digest in tail['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    assert tail['status']=='EXACT_BCH_TAIL_LOWER' and tail['rational_primal_dual_checks_passed']
    a0=tail['lower'];upper=sum(caps[w] for w in range(38,81,2));rows=[];close=F(1,a0)
    from audit_bch_q1_full_arb import rational
    def basis(n):
        return [[arb(christoffel.kraw(n,a,i)**2)/math.comb(n,a) for i in range(n+1)] for a in range(15)]
    for c in range(38,81,2):
        first=basis(c);second=basis(256-c)
        cumulative=[[sum((second[b][j] for b in range(15-a)),arb(0)) for j in range(257-c)] for a in range(15)]
        total=0
        for i in range(c+1):
            for j in range(257-c):
                a=i+j;b=c-i+j
                if not (38<=a<=80 and 38<=b<=80 and a%2==0 and b%2==0):continue
                kernel=sum((first[d][i]*cumulative[d][j] for d in range(15)),arb(0)).lower()
                assert kernel>0
                cap=math.floor(rational((arb(2)**128/kernel).upper()))
                total+=min(cap,math.comb(c,i)*math.comb(256-c,j),caps[a],caps[b])
        total=min(total,upper);close+=F(caps[c]*total,a0*a0)
        rows.append(dict(difference_weight=c,intersection_upper=total))
    assert close<F(1,128)
    if old:
        assert len(rows)==len(old['rows'])
        for row,saved in zip(rows,old['rows']):
            assert row['difference_weight']==saved['difference_weight'] and row['intersection_upper']<=saved['intersection_upper']
        assert close<=base.decode(old['close_pair_probability_upper'])
        print('512-bit BCH tail overlap replay passed',flush=True);return
    local=[Path(__file__),Path(christoffel.__file__),Path(cap_source.__file__),tail_path,
        base.HERE/'generated/christoffel_oa29_caps.json']+sources
    base.write_new(output,dict(status='OUTWARD_BCH_TAIL_OVERLAP',tail_weights=[38,80],difference_cutoff=80,
        rows=rows,tail_lower=a0,close_pair_probability_upper=base.encode(close),below_one_over_128=True,
        cap_receipts=[str(p.relative_to(base.HERE)) for p in sources],
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))
    print('Pr[wt(U+V)<=80] <1/128; diagnostic bits',math.log2(close.denominator)-math.log2(close.numerator),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
