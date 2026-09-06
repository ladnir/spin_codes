"""Exact full-coverage ledger for the fixed K24 construction."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import bridge as base
from frontier_sparse import parameters
from migrate_legacy_workspace import restore_or_verify


FILES=('frontier_k24_q1_v1.json','frontier_k24_q2_q1024_v1.json',
       'frontier_k24_q1025_q4095_v1.json','frontier_k24_dense_4096_32767_v2.json',
       'frontier_k24_dense_32768_98304_v1.json','frontier_k24_dense_98305_131072_v2.json')


def build():
    restore_or_verify(base.ROOT)
    manifest=base.read(base.HERE/'MIGRATION_MANIFEST.json')
    prefix='workstreams/bch_rm2sub_bridge/generated/joint_shell_'
    expected_caps={r['path'] for r in manifest['files'] if r['path'].startswith(prefix) and r['path'].endswith('/cap.json')}
    actual_caps={p.relative_to(base.ROOT).as_posix() for p in (base.HERE/'generated').glob('joint_shell_*/cap.json')}
    assert actual_caps==expected_caps and len(actual_caps)==46
    rows,cutoff=parameters(24);covered={};hashes={};summaries=[]
    for filename in FILES:
        path=base.HERE/'generated'/filename;data=base.read(path)
        replay_path=path.with_name(path.stem+'_replay.json');replay=base.read(replay_path)
        assert data['message_exponent']==replay['message_exponent']==24
        assert replay['producer_sha256']==base.sha(path)
        assert data['full_distance_proved'] is False
        for name,digest in data['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        if data['status']=='FRONTIER_Q1_OUTWARD_CERTIFICATE':
            assert replay['status']=='FRONTIER_Q1_512_BIT_REPLAY_PASSED' and replay['coefficients_checked']==92
            assert data['parameters']==dict(outer_rows=rows,cutoff=cutoff,step_bits=64,state_bits=20)
            co={int(w):base.decode(v) for w,v in data['coefficient_upper'].items()}
            assert sorted(co)==list(base.WEIGHTS)
            assert base.bch_bound(co)==tuple(base.decode(data[k]) for k in ('Q1_upper','factor','rest'))
            bounds={1:base.decode(data['Q1_upper'])}
        elif data['status']=='FRONTIER_SPARSE_OUTWARD_CERTIFICATE':
            assert replay['status']=='FRONTIER_SPARSE_512_BIT_REPLAY_PASSED'
            assert data['parameters']==dict(outer_rows=rows,cutoff=cutoff,step_bits=64,state_bits=20)
            lo,hi=data['occupancy_range'];assert replay['occupancy_range']==[lo,hi]
            assert len(data['upper_powers'])==hi-lo+1
            assert all(isinstance(v,int) and v<=-70 for v in data['upper_powers'])
            bounds={q:F(1,1<<(-v)) for q,v in zip(range(lo,hi+1),data['upper_powers'])}
        else:
            assert data['status']=='FRONTIER_DENSE_INTERVAL_CERTIFICATE'
            assert replay['status']=='FRONTIER_DENSE_512_BIT_REPLAY_PASSED'
            lo,hi=data['occupancy_range'];assert replay['occupancy_range']==[lo,hi]
            assert data['per_occupancy_upper_power']==-80
            assert replay['leaves']==data['leaves']
            pending=1;leaves=0
            for code in data['tree']:
                assert pending>0;pending-=1
                if code in ('q','p'):pending+=2
                else:assert isinstance(code,int);leaves+=1
            assert pending==0 and leaves==data['leaves']
            bounds={q:F(1,1<<80) for q in range(lo,hi+1)}
        total=sum(bounds.values(),F(0))
        if 'range_upper' in data:assert total==base.decode(data['range_upper'])
        assert not set(covered).intersection(bounds) and all(v>0 for v in bounds.values())
        covered.update(bounds)
        summaries.append(dict(file=filename,occupancy_range=[min(bounds),max(bounds)],upper=base.encode(total)))
        for source in (path,replay_path):hashes[source.relative_to(base.ROOT).as_posix()]=base.sha(source)
    assert sorted(covered)==list(range(1,rows+1))
    total=sum(covered.values(),F(0));higher=total-covered[1]
    assert total<F(1,1<<46)<F(1,1<<40)
    return dict(status='FULL_K24_DISTANCE_CERTIFICATE_LEDGER',full_distance_proved=True,
        configuration='t64_s20',message_bits=1<<24,output_bits=1<<25,outer_rows=rows,
        cutoff=cutoff,minimum_distance_outside_bad_setup=cutoff+1,covered_occupancy_range=[1,rows],
        failure_upper=base.encode(total),higher_occupancy_upper=base.encode(higher),
        margin_bits=math.log2(total.denominator)-math.log2(total.numerator),
        higher_occupancy_margin_bits=math.log2(higher.denominator)-math.log2(higher.numerator),
        scope='Fixed BCH [256,128] and selected t64_s20 map; independent permutations and fresh nonzero field multipliers.',
        certificates=summaries,input_sha256=hashes,audit_source_sha256=base.sha(Path(__file__)))


def run(output,verify=False):
    old=base.read(output) if verify else None
    if not verify:assert not output.exists()
    payload=build()
    if verify:
        assert payload==old
        print('Exact full ledger replay passed',flush=True)
    else:base.write_new(output,payload)
    print('FULL K24 margin',payload['margin_bits'],'higher-occupancy margin',payload['higher_occupancy_margin_bits'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.output,a.verify)
