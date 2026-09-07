"""Balanced subdivision of the frozen two-tilt bound and partition checker."""
import argparse
from fractions import Fraction as F
import heapq
import math
from pathlib import Path
import time

from flint import ctx

import search_ladder_dense_tilt as old

base,core,dense = old.base,old.core,old.dense
POLICY = 'two-tilt-balanced-input-density-width-v1'


def split_axis(node,checker):
    lo,hi,a,b = old.geometry(node)
    if lo == hi:
        return 'v'
    nlo,nhi = checker.density(a,b)
    # Compare the two contributions to theta=Q*nu/L variation exactly.
    return 'q' if (hi-lo)*(nlo+nhi) >= (hi+lo)*(nhi-nlo) else 'v'


def run(directory,m,minimum,seconds,nodes,verify=False):
    directory = directory.resolve()
    spec = dense.labels.identity.instance(m)
    core.require(1 <= minimum <= spec['rows'] and seconds > 0 and nodes >= 0,'Invalid run geometry/budget')
    directory.mkdir(parents=True,exist_ok=True)
    paths = sorted(directory.glob('cover_*.json'))
    saved = base.read(paths[-1]) if paths else None
    core.require(not verify or saved is not None,'No cover to replay')
    if saved:
        core.authenticate(saved,base.ROOT)
        core.require(saved['instance'] == spec and saved['minimum'] == minimum and saved['policy'] == POLICY,
                     'Changed cover instance or subdivision policy')
    ctx.prec = 512 if verify else 256
    ps = [base.decode(p) for p in saved['probabilities']] if saved else dense.labels.row_probabilities(.5)
    checker = dense.Checker(spec,ps)
    if saved:
        leaves,splits = saved['leaves'],saved['splits']
        old.check_partition(leaves,splits,minimum,spec['rows'])
    else:
        root = old.box(minimum,spec['rows'],F(0),F(1),checker.witness(minimum,spec['rows'],F(0),F(1)))
        root['power'] = old.retained_power(checker,root)
        leaves,splits = {'':root},{}
    if verify:
        for i,(key,node) in enumerate(leaves.items()):
            core.require(old.retained_power(checker,node) <= node['power'],'Outward leaf replay failed: '+key)
            if i % 250 == 0:
                print('replay',i+1,'/',len(leaves),flush=True)
        unresolved = sum(node['power'] > -80 for node in leaves.values())
        core.require(unresolved == saved['unresolved'],'Wrong unresolved count')
        if not unresolved:
            total = sum((F(2)**node['power'] for node in leaves.values()),F(0))
            core.require(total == base.decode(saved['union_upper']),'Wrong exact union')
        base.write_new(directory/('replay_'+paths[-1].stem+'.json'),
                       dict(status='TWO_TILT_512_BIT_REPLAY_PASSED',producer_sha256=base.sha(paths[-1]),
                            source_sha256=old.provenance.sources(),complete_dense_certificate=not unresolved))
        print('Replayed',len(leaves),'leaves; unresolved',unresolved,flush=True)
        return
    heap = [(-node['power'],key) for key,node in leaves.items() if node['power'] > -80]
    heapq.heapify(heap)
    started = time.monotonic()
    used = 0
    while heap and used < nodes and time.monotonic()-started < seconds:
        _,key = heapq.heappop(heap)
        node = leaves[key]
        candidate = dict(node,witness=checker.witness(*old.geometry(node)))
        candidate['power'] = old.retained_power(checker,candidate)
        if candidate['power'] < node['power']:
            node = leaves[key] = candidate
        used += 1
        if node['power'] > -80:
            axis = split_axis(node,checker)
            branch = old.children(node,axis)
            del leaves[key]
            splits[key] = axis
            for i,child in enumerate(branch):
                child['power'] = old.retained_power(checker,child)
                child_key = key+str(i)
                leaves[child_key] = child
                if child['power'] > -80:
                    heapq.heappush(heap,(-child['power'],child_key))
        if used % 100 == 0:
            print('refined',used,'leaves',len(leaves),'unresolved',len(heap),
                  'worst exponent',-heap[0][0] if heap else None,flush=True)
    coverage = old.check_partition(leaves,splits,minimum,spec['rows'])
    unresolved = sum(node['power'] > -80 for node in leaves.values())
    result = dict(status='TWO_TILT_DENSE_RANGE_CERTIFIED' if not unresolved else 'TWO_TILT_COMPLETE_COVER_UNRESOLVED',
                  policy=POLICY,instance=spec,minimum=minimum,probabilities=[base.encode(p) for p in ps],
                  leaves=leaves,splits=splits,coverage=coverage,unresolved=unresolved,
                  refinements_this_run=used,seconds=time.monotonic()-started,
                  full_distance_proved=False,source_sha256=old.provenance.sources())
    if not unresolved:
        total = sum((F(2)**node['power'] for node in leaves.values()),F(0))
        result.update(union_upper=base.encode(total),margin_bits=math.log2(total.denominator)-math.log2(total.numerator))
    output = directory/f'cover_{len(paths):04d}.json'
    base.write_new(output,result)
    print('saved',output,'unresolved',unresolved,'margin',result.get('margin_bits'),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--m',type=int,default=20)
    p.add_argument('--minimum',type=int,default=512)
    p.add_argument('--seconds',type=float,default=180)
    p.add_argument('--nodes',type=int,default=2000)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.directory,a.m,a.minimum,a.seconds,a.nodes,a.verify)
