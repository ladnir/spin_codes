"""Bounded candidate dense cover; old trees supply geometry, never bounds."""
import argparse
from fractions import Fraction as F
import heapq
from pathlib import Path
import time
from flint import ctx
import candidate
import certify_dense
import search_ladder_dense_fixed as splitting
model = candidate.model
geometry = certify_dense.geometry


def run(output,name,exponent,minimum,nodes,seconds,seed,verify):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        name,exponent,minimum = (saved[k] for k in ('candidate','exponent','minimum'))
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    checker = candidate.Checker(exponent,name)
    def bound(node):
        return certify_dense.power_bound(checker.bound(*geometry.geometry(node),node['witness']))
    if saved:
        assert saved['instance'] == checker.engine.identity()
        leaves,splits = saved['leaves'],saved['splits']
    elif seed:
        prior = model.base.read(seed)
        model.authenticate(prior)
        assert prior.get('exponent',prior.get('message_exponent')) == exponent
        assert prior['minimum'] == minimum
        leaves,splits = prior['leaves'],prior['splits']
        geometry.check_partition(leaves,splits,minimum,checker.rows)
        for node in leaves.values():
            node['power'] = bound(node)
    else:
        root = geometry.box(minimum,checker.rows,F(0),F(1),checker.witness(minimum,checker.rows,F(0),F(1)))
        root['power'] = bound(root)
        leaves,splits = {'':root},{}
    geometry.check_partition(leaves,splits,minimum,checker.rows)
    used = 0
    if saved:
        for key,node in leaves.items():
            assert bound(node) <= node['power'],key
    else:
        heap = [(-node['power'],key) for key,node in leaves.items() if node['power'] > -80]
        heapq.heapify(heap)
        print('initial cover',name,exponent,'leaves',len(leaves),'weak',len(heap),flush=True)
        started = time.monotonic()
        while heap and used < nodes and time.monotonic()-started < seconds:
            _,key = heapq.heappop(heap)
            node = leaves[key]
            proposal = dict(node,witness=checker.witness(*geometry.geometry(node)))
            proposal['power'] = bound(proposal)
            if proposal['power'] < node['power']:
                node = leaves[key] = proposal
            used += 1
            if node['power'] > -80:
                lo,hi,a,b = geometry.geometry(node)
                if lo == hi and b-a <= F(1,1 << 18):
                    continue
                axis = splitting.split_axis(node,checker)
                splits[key] = axis
                del leaves[key]
                for i,child in enumerate(geometry.children(node,axis)):
                    child['power'] = bound(child)
                    child_key = key+str(i)
                    leaves[child_key] = child
                    if child['power'] > -80:
                        heapq.heappush(heap,(-child['power'],child_key))
            if used%10 == 0:
                print('cover',exponent,used,'leaves',len(leaves),'weak',len(heap),flush=True)
    coverage = geometry.check_partition(leaves,splits,minimum,checker.rows)
    unresolved = sum(node['power'] > -80 for node in leaves.values())
    total = sum((F(2)**node['power'] for node in leaves.values()),F(0)) if not unresolved else None
    if saved:
        assert unresolved == saved['unresolved']
        if total is not None:
            assert total == model.base.decode(saved['upper'])
        model.base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='CANDIDATE_DENSE_COVER_REPLAY_PASSED',producer_sha256=model.base.sha(output),
            complete_dense_certificate=not unresolved,full_distance_proved=False))
    else:
        sources = candidate.sources()
        if seed:
            sources[seed.relative_to(model.ROOT).as_posix()] = model.base.sha(seed)
        result = dict(status='OUTWARD_CANDIDATE_DENSE_COVER',candidate=name,exponent=exponent,
            minimum=minimum,instance=checker.engine.identity(),leaves=leaves,splits=splits,
            coverage=coverage,unresolved=unresolved,refinements=used,source_sha256=sources,full_distance_proved=False)
        if total is not None:
            result['upper'] = model.base.encode(total)
        model.base.write_new(output,result)
    print('cover done',exponent,'unresolved',unresolved,flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--candidate',choices=('weight5_seed0','weight5_seed1','balanced'),default='weight5_seed0')
    p.add_argument('--m',type=int,choices=(16,18,20),default=20)
    p.add_argument('--minimum',type=int,default=512)
    p.add_argument('--nodes',type=int,default=200)
    p.add_argument('--seconds',type=float,default=180)
    p.add_argument('--seed',type=Path)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    assert a.nodes >= 0 and a.seconds > 0
    run(a.output.resolve(),a.candidate,a.m,a.minimum,a.nodes,a.seconds,a.seed.resolve() if a.seed else None,a.verify)
