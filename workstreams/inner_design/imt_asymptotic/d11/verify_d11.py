"""Check composition of the refined 11% dense certificate and old sparse gates."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

from flint import ctx
import outer_majorant as outer
import dense_refined as refined

screen=outer.screen
dense=outer.dense
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import verify_imt_asymptotic as base
import certify_imt_one_two as one_two


def main():
    if not __debug__:raise RuntimeError('Do not use -O.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    paths=[HERE/name for name in ('OUTER_REFINED.json','OUTER_REPLAY.json','DENSE_OUTWARD.json','DENSE_REPLAY.json')]
    paths += [HERE.parent/name for name in ('ASSEMBLY_D1099.json','SPARSE_EXACT.json','SPARSE_EXACT_REPLAY.json','FIXED_LIMIT.json','ONE_TWO_OUTWARD.json')]
    records=[json.loads(p.read_text()) for p in paths]
    for r in records:refined.authenticate(r)
    new,outer_replay,cert,replay,assembly,sparse,sparse_replay,fixed,small=records
    base.reproduce.check_asymptotic()
    expected,new_segments=outer.supports()
    assert new['status']==outer_replay['status']=='proved' and outer_replay['replay_passed'] is True
    assert new['left_segments']==outer_replay['left_segments']==[
        dict(weight=[str(s.omega_lo),str(s.omega_hi)],slope=str(s.slope),intercept=str(s.intercept)) for s in expected[:-1]]
    assert new['central_constant_upper']==str(expected[-1].intercept)
    assert new['new_support_checks']==outer_replay['new_support_checks']
    assert set(new['new_support_checks'])=={str(s.index) for s in new_segments}
    # Fresh high-working-precision leaf replay and exact partition coverage.
    outer.old.iv.dps=100
    for segment in new_segments:
        outer.replay_partition(segment,new['new_support_checks'][str(segment.index)])
    segments=screen.load_segments(paths[0])
    assert len(segments)==39
    for record in (cert,replay):
        assert record['delta']=='11/100'
        assert record['outer_sha256']==screen.sha(paths[0])
        assert record['alpha_range']==['1/10000','1']
        assert record['row_density_range']==['13/125','112/125']
        assert record['boxes']==len(record['leaves'])==1023
        dense.check_geometry(record['leaves'],segments)
        assert record['inner']==assembly['inner']==sparse['inner']==sparse_replay['inner']
    assert cert['precision_bits']==256 and replay['precision_bits']==512
    for a,b in zip(cert['leaves'],replay['leaves']):
        assert all(a[k]==b[k] for k in ('segment','alpha','row_density','rational_witness'))
        assert F(b['exponent_upper'])<=F(a['exponent_upper'])<0
    assert sparse['row_checks']==sparse_replay['row_checks']
    assert sparse['alpha_interval']==['0','1/10000'] and sparse['contraction']=='1-96*alpha'
    assert all(F(r['upper'])<0 for r in sparse['row_checks'])
    ctx.prec=512;num=dense.number
    coefficient=(num(2).log()/2+num(F(5,8)).log()
                 +num(F(3,5)+F(11,100)*F(8,5))/num(1-F(8,5)*F(1,10000))
                 -num(F(96,128))+num(F(1281,100000))+num(2).log()/num(F(39,4)))
    assert dense.upper(coefficient)<=F(sparse['conditional_structured_exponent_coefficient_upper'])<0
    assert assembly['sparse_cutoff']==4096
    assert dense.upper(coefficient+num(512).log()/4096+num(F(4,1000))) < -F(6,1000)
    old_segments=screen.load_segments(screen.FROZEN/'golay_ba3_concave_majorant.json')
    base.check_fixed(fixed,old_segments)
    assert small['delta']=='11/100' and small['block_constant']=='39/4'
    for a,b in zip(small['checks'],one_two.evaluate(512)):
        assert F(b['exponent_per_outer_bit_upper'])<=F(a['exponent_per_outer_bit_upper'])<0
    paths += [Path(__file__),HERE/'PROOF_UPDATE.md',HERE/'test_refined.py']
    payload=dict(status='IMT_D11_PROOF_DRAFT_NUMERICAL_ASSEMBLY_PASSED',delta='11/100',
                 block_constant='39/4',inner=cert['inner'],dense_boxes=1023,
                 dense_maximum_upper=str(max(F(b['exponent_upper']) for b in cert['leaves'])),
                 new_outer_supports=5,outer_segments=39,
                 outer_leaves=sum(map(len,new['new_support_checks'].values())),
                 sparse_coefficient_upper=str(dense.upper(coefficient)),
                 sparse_cutoff=4096,fresh_outer_leaf_replay=True,
                 dense_replay_authenticated=True,analytic_proof_machine_checked=False,
                 independent_analytic_review_pending=True,paper_or_default_changed=False,
                 source_sha256={p.relative_to(screen.ROOT).as_posix():screen.sha(p) for p in paths})
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(payload,stream,indent=2);stream.write('\n')
    print('PASS: 11% numerical assembly, all occupancy regimes. Analytic review remains separate.',flush=True)


if __name__=='__main__':main()
