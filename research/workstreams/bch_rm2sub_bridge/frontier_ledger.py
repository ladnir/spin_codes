"""Exact certificate aggregation with storage proportional to receipts, not K."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import bridge as base
from frontier_sparse import parameters
from migrate_legacy_workspace import restore_or_verify


def validate_tree(tree,lo,hi):
    stack=[(lo,hi)];leaves=0
    for code in tree:
        assert stack,'Trailing tree nodes'
        a,b=stack.pop()
        if code=='q':
            assert a<b,'Cannot split a single integer occupancy'
            mid=(a+b)//2;stack.extend(((mid+1,b),(a,mid)))
        elif code=='p':stack.extend(((a,b),(a,b)))
        else:assert type(code) is int;leaves+=1
    assert not stack,'Incomplete tree'
    return leaves


def validate_coverage(intervals,rows):
    following=1
    for lo,hi in sorted(intervals):
        assert lo==following and lo<=hi<=rows,'Gap, overlap, or invalid occupancy interval'
        following=hi+1
    assert following==rows+1,'Incomplete final coverage'


def build(m,paths):
    restore_or_verify(base.ROOT)
    manifest=base.read(base.HERE/'MIGRATION_MANIFEST.json')
    prefix='workstreams/bch_rm2sub_bridge/generated/joint_shell_'
    expected={r['path'] for r in manifest['files'] if r['path'].startswith(prefix) and r['path'].endswith('/cap.json')}
    actual={p.relative_to(base.ROOT).as_posix() for p in (base.HERE/'generated').glob('joint_shell_*/cap.json')}
    assert actual==expected and len(actual)==46
    rows,cutoff=parameters(m);summaries=[];hashes={};q1=None
    for path in paths:
        path=path.resolve();data=base.read(path);replay_path=path.with_name(path.stem+'_replay.json');replay=base.read(replay_path)
        assert data['message_exponent']==replay['message_exponent']==m
        assert replay['producer_sha256']==base.sha(path) and data['full_distance_proved'] is False
        for name,digest in data['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        status=data['status']
        if status=='FRONTIER_Q1_OUTWARD_CERTIFICATE':
            assert q1 is None and replay['status']=='FRONTIER_Q1_512_BIT_REPLAY_PASSED'
            assert replay['coefficients_checked']==92
            assert data['parameters']==dict(outer_rows=rows,cutoff=cutoff,step_bits=64,state_bits=20)
            co={int(w):base.decode(v) for w,v in data['coefficient_upper'].items()}
            assert sorted(co)==list(base.WEIGHTS)
            assert base.bch_bound(co)==tuple(base.decode(data[k]) for k in ('Q1_upper','factor','rest'))
            lo=hi=1;q1=total=base.decode(data['Q1_upper'])
        elif status=='FRONTIER_SPARSE_OUTWARD_CERTIFICATE':
            assert replay['status']=='FRONTIER_SPARSE_512_BIT_REPLAY_PASSED'
            assert data['parameters']==dict(outer_rows=rows,cutoff=cutoff,step_bits=64,state_bits=20)
            lo,hi=data['occupancy_range'];assert replay['occupancy_range']==[lo,hi]
            powers=data['upper_powers'];assert len(powers)==hi-lo+1
            assert all(type(v) is int and -1000<=v<0 for v in powers)
            total=sum((F(1,1<<(-v)) for v in powers),F(0))
        else:
            assert status=='FRONTIER_DENSE_INTERVAL_CERTIFICATE'
            assert replay['status']=='FRONTIER_DENSE_512_BIT_REPLAY_PASSED'
            lo,hi=data['occupancy_range'];assert replay['occupancy_range']==[lo,hi]
            power=data['per_occupancy_upper_power'];assert type(power) is int and -1000<=power<0
            assert validate_tree(data['tree'],lo,hi)==data['leaves']==replay['leaves']
            total=F(hi-lo+1,1<<(-power))
        assert total>0
        if 'range_upper' in data:assert total==base.decode(data['range_upper'])
        summaries.append(dict(file=path.relative_to(base.ROOT).as_posix(),occupancy_range=[lo,hi],upper=base.encode(total)))
        for source in (path,replay_path):hashes[source.relative_to(base.ROOT).as_posix()]=base.sha(source)
    validate_coverage([r['occupancy_range'] for r in summaries],rows)
    assert q1 is not None
    total=sum((base.decode(r['upper']) for r in summaries),F(0));higher=total-q1
    assert 0<higher and total<F(1,1<<40)
    bits=total.denominator.bit_length()-total.numerator.bit_length()
    if total>=F(1,1<<bits):bits-=1
    assert total<F(1,1<<bits) and total>=F(1,1<<(bits+1))
    margin=lambda v:math.log2(v.denominator)-math.log2(v.numerator)
    return dict(status='FULL_FRONTIER_DISTANCE_CERTIFICATE_LEDGER',full_distance_proved=True,configuration='t64_s20',
        message_exponent=m,message_bits=1<<m,output_bits=1<<(m+1),outer_rows=rows,cutoff=cutoff,
        minimum_distance_outside_bad_setup=cutoff+1,covered_occupancy_range=[1,rows],
        failure_upper=base.encode(total),higher_occupancy_upper=base.encode(higher),margin_bits=margin(total),
        higher_occupancy_margin_bits=margin(higher),strict_certified_bits=bits,
        scope='Fixed BCH [256,128], selected t64_s20, independent permutations and fresh nonzero field multipliers.',
        certificates=sorted(summaries,key=lambda r:r['occupancy_range']),input_sha256=hashes,
        audit_source_sha256=base.sha(Path(__file__)))


def run(output,m,paths,verify=False):
    old=base.read(output) if verify else None
    if old:
        m=old['message_exponent'];paths=[base.ROOT/r['file'] for r in old['certificates']]
    else:assert not output.exists() and paths
    payload=build(m,paths)
    if old:assert payload==old;print('Exact interval-ledger replay passed',flush=True)
    else:base.write_new(output,payload)
    print('FULL K',m,'margin',payload['margin_bits'],'higher margin',payload['higher_occupancy_margin_bits'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m',type=int,default=26);p.add_argument('--certificates',type=Path,nargs='+')
    p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.output,a.m,a.certificates,a.verify)
