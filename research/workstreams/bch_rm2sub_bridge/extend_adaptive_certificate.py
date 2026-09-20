"""Reusable outward range extension; never overwrites frozen certificates."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path

from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import general_occupancy as general
import certify_adaptive_range as previous
from general_batch_certificate import density_cost


def build(lower,upper,verify=False,tightened=False):
    assert 1<=lower<=upper<=8192
    import tightened_occupancy as tight
    model=tight if tightened else general
    prefix='tightened' if tightened else 'adaptive'
    suffix='screen' if tightened else 'refined_screen'
    screen_path=base.HERE/'generated'/f'{prefix}_q{lower}_q{upper}_{suffix}.json'
    output=base.HERE/'generated'/f'{prefix}_q{lower}_q{upper}_extension_outward.json'
    if not verify:assert not output.exists()
    old=base.read(output) if verify else None
    if old:
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        for name,digest in old['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    screen=base.read(screen_path)
    for name,digest in screen['source_sha256'].items():assert base.sha(base.HERE/name)==digest
    assert [r['occupation'] for r in screen['rows']]==list(range(lower,upper+1))
    maxima={}
    for row in screen['rows']:
        j=row['witness_tenth'];assert isinstance(j,int) and -120<=j<=0
        maxima[j]=max(maxima.get(j,0),row['occupation'])
    t,s,spectrum=base.load_map(general.NAME);kernel=general.kernel_spectrum();caps=general.q2.deterministic_caps()
    ctx.prec=512 if verify else 256;cache={};total=F(0);rows=[]
    from audit_bch_q1_full_arb import rational
    for row in screen['rows']:
        q=row['occupation'];j=row['witness_tenth']
        if j not in cache:
            lam=(arb(j)/10).exp()
            cache[j]=(model.regions(t,s,spectrum,kernel,(-lam).exp(),arb,maxima[j]),(209716*lam).exp())
        region,correction=cache[j]
        ps=[base.decode(v) for v in row['p']];assert len(ps)==5 and all(0<p<1 for p in ps)
        gamma=[density_cost(band,p,caps) for band,p in zip(general.BANDS,ps)]
        roots=[((arb(c.numerator)/c.denominator).log()/256).exp().upper() for c in gamma]
        matrix=previous.adaptive_kernel(region[:q+1],[arb(p.numerator)/p.denominator for p in ps],roots)
        for _ in range(8):matrix=q1.positive_mul(matrix,matrix)
        bound=math.comb(8192,q)*5**q*rational((sum(matrix[:3],arb(0))*correction).upper())
        assert bound>0
        if old:assert bound<=base.decode(old['rows'][q-lower]['upper'])
        total+=bound
        rows.append(dict(occupation=q,upper=base.encode(bound),margin_bits_diagnostic=math.log2(bound.denominator)-math.log2(bound.numerator),
                         witness_tenth=j,p=row['p']))
        if q%8==0 or q==lower:print('Q',q,'replayed' if verify else 'outward','margin',rows[-1]['margin_bits_diagnostic'],flush=True)
    if old:
        assert total<=base.decode(old['range_upper'])
        print('512-bit range replay passed',lower,upper,flush=True);return
    local,outer=general.q2.source_paths()
    local += [Path(__file__),Path(general.__file__),Path(previous.__file__),base.HERE/'general_batch_certificate.py',
              base.HERE/'certify_q3_compact.py',screen_path]
    local += [base.HERE/name for name in screen['source_sha256']]
    if tightened:local.append(Path(tight.__file__))
    payload=dict(status='OUTWARD_ADAPTIVE_RANGE_EXTENSION',configuration=general.NAME,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,
                        step_bits=t,state_bits=s,distance_cutoff=209716),occupancy_range=[lower,upper],rows=rows,
        range_upper=base.encode(total),range_below_2_to_minus_60=total<F(1,1<<60),all_occupations_certified=False,
        tightened_termination=tightened,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local},
        outer_sha256={str(p.relative_to(base.BCH)):base.sha(p) for p in outer})
    base.write_new(output,payload)
    print('Range',lower,upper,'margin',math.log2(total.denominator)-math.log2(total.numerator),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--lower',required=True,type=int);parser.add_argument('--upper',required=True,type=int)
    parser.add_argument('--verify',action='store_true');parser.add_argument('--tightened',action='store_true')
    args=parser.parse_args();build(args.lower,args.upper,args.verify,args.tightened)
