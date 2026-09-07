"""Memory-bounded complete rectangle covers for the two-tilt dense bound.

Every checkpoint is a full partition, including explicitly unresolved boxes.
Only a cover with no unresolved boxes receives an exact dense-range union.
"""
import argparse
from fractions import Fraction as F
import heapq
import math
from pathlib import Path
import time

from flint import arb,ctx

import ladder_dense_tilt as dense
import bch256_dense_types as provenance

base,core = dense.base,dense.core


def geometry(node):
    return node['lo'],node['hi'],base.decode(node['vlo']),base.decode(node['vhi'])


def box(lo,hi,vlo,vhi,witness):
    return dict(lo=lo,hi=hi,vlo=base.encode(vlo),vhi=base.encode(vhi),witness=witness)


def children(node,axis):
    lo,hi,a,b = geometry(node)
    witness = node['witness']
    if axis == 'q':
        core.require(lo < hi,'Cannot split singleton occupancy')
        middle = (lo+hi)//2
        return [box(lo,middle,a,b,witness),box(middle+1,hi,a,b,witness)]
    core.require(axis == 'v' and a < b,'Invalid density split')
    middle = (a+b)/2
    return [box(lo,hi,a,middle,witness),box(lo,hi,middle,b,witness)]


def check_partition(leaves,splits,minimum,rows):
    pending = [('',box(minimum,rows,F(0),F(1),None))]
    seen_leaves,seen_splits = set(),set()
    while pending:
        key,expected = pending.pop()
        if key in leaves:
            core.require(key not in splits and geometry(leaves[key]) == geometry(expected),'Wrong leaf geometry')
            seen_leaves.add(key)
        else:
            core.require(key in splits,'Missing partition branch')
            seen_splits.add(key)
            pending.extend((key+str(i),child) for i,child in enumerate(children(expected,splits[key])))
    core.require(seen_leaves == set(leaves) and seen_splits == set(splits),'Extraneous partition nodes')
    return dict(complete=True,occupancy_interval=[minimum,rows],leaves=len(leaves),splits=len(splits))


def retained_power(checker,node):
    value = checker.bound(*geometry(node),node['witness'])
    bits = dense.labels.intervals.rational((value/arb(2).log()).upper())
    return max(-200,-(-bits.numerator//bits.denominator))


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
        core.require(saved['instance'] == spec and saved['minimum'] == minimum,'Changed cover instance')
    ctx.prec = 512 if verify else 256
    ps = [base.decode(p) for p in saved['probabilities']] if saved else dense.labels.row_probabilities(.5)
    checker = dense.Checker(spec,ps)
    if saved:
        leaves,splits = saved['leaves'],saved['splits']
        check_partition(leaves,splits,minimum,spec['rows'])
    else:
        root = box(minimum,spec['rows'],F(0),F(1),checker.witness(minimum,spec['rows'],F(0),F(1)))
        root['power'] = retained_power(checker,root)
        leaves,splits = {'':root},{}
    if verify:
        for i,(key,node) in enumerate(leaves.items()):
            core.require(retained_power(checker,node) <= node['power'],'Outward leaf replay failed: '+key)
            if i % 250 == 0:
                print('replay',i+1,'/',len(leaves),flush=True)
        unresolved = sum(node['power'] > -80 for node in leaves.values())
        core.require(unresolved == saved['unresolved'],'Wrong unresolved count')
        if not unresolved:
            total = sum((F(2)**node['power'] for node in leaves.values()),F(0))
            core.require(total == base.decode(saved['union_upper']),'Wrong exact union')
        output = paths[-1].with_name(paths[-1].stem+'_replay.json')
        # Replay files must not match the numeric checkpoint pattern on resume.
        output = directory/('replay_'+paths[-1].stem+'.json')
        base.write_new(output,dict(status='TWO_TILT_512_BIT_REPLAY_PASSED',producer_sha256=base.sha(paths[-1]),
                       source_sha256=provenance.sources(),complete_dense_certificate=not unresolved))
        print('Replayed',len(leaves),'leaves; unresolved',unresolved,flush=True)
        return
    heap = [(-node['power'],key) for key,node in leaves.items() if node['power'] > -80]
    heapq.heapify(heap)
    started = time.monotonic()
    used = 0
    while heap and used < nodes and time.monotonic()-started < seconds:
        _,key = heapq.heappop(heap)
        node = leaves[key]
        new_witness = checker.witness(*geometry(node))
        candidate = dict(node,witness=new_witness)
        candidate['power'] = retained_power(checker,candidate)
        if candidate['power'] < node['power']:
            node = leaves[key] = candidate
        used += 1
        if node['power'] > -80:
            choices = []
            for axis in ('q','v') if node['lo'] < node['hi'] else ('v',):
                branch = children(node,axis)
                for child in branch:
                    child['power'] = retained_power(checker,child)
                choices.append((max(c['power'] for c in branch),sum(c['power'] for c in branch),axis,branch))
            _,_,axis,branch = min(choices,key=lambda c:c[:2])
            del leaves[key]
            splits[key] = axis
            for i,child in enumerate(branch):
                child_key = key+str(i)
                leaves[child_key] = child
                if child['power'] > -80:
                    heapq.heappush(heap,(-child['power'],child_key))
        if used % 100 == 0:
            print('refined',used,'leaves',len(leaves),'unresolved',len(heap),
                  'worst exponent',-heap[0][0] if heap else None,flush=True)
    coverage = check_partition(leaves,splits,minimum,spec['rows'])
    unresolved = sum(node['power'] > -80 for node in leaves.values())
    result = dict(status='TWO_TILT_DENSE_RANGE_CERTIFIED' if not unresolved else 'TWO_TILT_COMPLETE_COVER_UNRESOLVED',
                  instance=spec,minimum=minimum,probabilities=[base.encode(p) for p in ps],leaves=leaves,splits=splits,
                  coverage=coverage,unresolved=unresolved,refinements_this_run=used,seconds=time.monotonic()-started,
                  full_distance_proved=False,source_sha256=provenance.sources())
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
