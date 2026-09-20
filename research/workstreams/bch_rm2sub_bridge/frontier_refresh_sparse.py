"""Size-parameterized refresh-aware Q2..Q7 certificates and independent replay."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import numpy as np
from flint import arb,ctx
import k30_sparse as core
import frontier_sparse as frontier


def run(output,m,first,last,verify=False):
    base=core.base;output=output.resolve();old=base.read(output) if verify else None
    if old:
        assert old['method']=='refresh_four_state_kernel_exclusion'
        m=old['message_exponent'];first,last=old['occupancy_range']
        for name,digest in old['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
    else:assert not output.exists()
    assert 2<=first<=last<=7
    core.restore_or_verify(base.ROOT);ctx.prec=512 if verify else 256
    rows,cutoff=frontier.parameters(m)
    t,s,spectrum,kernel=core.maps.load('t64_s20');counts=core.caps()
    results=[];powers=[]
    for q in range(first,last+1):
        if old:
            saved=old['rows'][q-first];assert saved['occupation']==q
            candidates=[(base.decode(saved['scaled_tilt']),[base.decode(p) for p in saved['p']])]
        else:candidates=[(F(q*n,2),None) for n in (3,5,8,12,18,26,40)]
        best=None
        for a,ps in candidates:
            lam=(arb(a.numerator)/a.denominator)/rows
            region=core.regions(t,s,spectrum,kernel,(-lam).exp(),q,length=rows,reverse=verify)
            if ps is None:
                approx=np.array([[float(v) for v in mat] for mat in region]).reshape(q+1,4,4)
                ps=core.choose_ps(approx,q,counts)
            mat=core.adaptive(region,q,ps,counts,arb)
            for _ in range(8):mat=core.matrix.product(mat,mat)
            upper=math.comb(rows,q)*len(core.BANDS)**q*core.rational((sum(mat[:4],arb(0))*(cutoff*lam).exp()).upper())
            if best is None or upper<best[0]:best=upper,a,ps
        upper,a,ps=best
        assert 0<upper<F(1,1<<44)
        power=upper.numerator.bit_length()-upper.denominator.bit_length()
        if upper>frontier.as_fraction(power):power+=1
        power=max(-80,power)
        assert upper<=frontier.as_fraction(power)
        if old:
            assert upper<=base.decode(saved['upper'])
            assert upper<=frontier.as_fraction(old['upper_powers'][q-first])
        results.append(dict(occupation=q,upper=base.encode(upper),scaled_tilt=base.encode(a),
            p=[base.encode(p) for p in ps],margin_bits=math.log2(upper.denominator)-math.log2(upper.numerator)))
        powers.append(power)
        print('replay' if old else 'producer','Q',q,'margin',results[-1]['margin_bits'],flush=True)
    if old:
        assert old['parameters']==dict(outer_rows=rows,cutoff=cutoff,step_bits=t,state_bits=s)
        assert len(old['rows'])==len(old['upper_powers'])==last-first+1
        assert sum((frontier.as_fraction(p) for p in old['upper_powers']),F(0))==base.decode(old['range_upper'])
        base.write_new(output.with_name(output.stem+'_replay.json'),dict(status='FRONTIER_SPARSE_512_BIT_REPLAY_PASSED',
            method=old['method'],producer_sha256=base.sha(output),message_exponent=m,occupancy_range=[first,last]))
    else:
        hashes=frontier.sources()
        for path in (Path(__file__),Path(core.__file__),Path(core.matrix.__file__),base.HERE/'certify_k30_q1.py'):
            hashes[path.relative_to(base.ROOT).as_posix()]=base.sha(path)
        base.write_new(output,dict(status='FRONTIER_SPARSE_OUTWARD_CERTIFICATE',method='refresh_four_state_kernel_exclusion',
            message_exponent=m,parameters=dict(outer_rows=rows,cutoff=cutoff,step_bits=t,state_bits=s),
            occupancy_range=[first,last],rows=results,upper_powers=powers,
            range_upper=base.encode(sum((frontier.as_fraction(p) for p in powers),F(0))),
            source_sha256=hashes,full_distance_proved=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m',type=int,default=28);p.add_argument('--first',type=int,default=2);p.add_argument('--last',type=int,default=7)
    p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.output,a.m,a.first,a.last,a.verify)
