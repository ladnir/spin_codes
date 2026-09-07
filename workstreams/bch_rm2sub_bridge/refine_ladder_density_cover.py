"""Exact density-factor refinement of a fully replayed rectangle cover.

Unresolved source leaves are permitted: they are finite valid upper bounds.
The new union covers the same full dense domain with D_L instead of L+1.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import search_ladder_dense_fixed as search
import poisson_density_factor as factor

base,core = search.base,search.core


def run(source,output,verify=False):
    receipt = base.read(output) if verify else None
    if receipt:
        core.authenticate(receipt,base.ROOT)
        source = base.ROOT/receipt['source_cover']
        core.require(base.sha(source) == receipt['source_cover_sha256'],'Source cover changed')
    else:
        core.require(source is not None and not output.exists(),'Provide a source and fresh output')
    saved = base.read(source)
    replay_path = source.with_name('replay_'+source.stem+'.json')
    replay = base.read(replay_path)
    core.authenticate(saved,base.ROOT)
    core.authenticate(replay,base.ROOT)
    spec = saved['instance']
    core.require(spec == search.dense.labels.identity.instance(spec['message_exponent']) and
                 saved['policy'] == search.POLICY,'Wrong source instance or bounding policy')
    core.require(replay['status'] == 'TWO_TILT_512_BIT_REPLAY_PASSED' and
                 replay['producer_sha256'] == base.sha(source),'Missing numerical replay')
    coverage = search.old.check_partition(saved['leaves'],saved['splits'],saved['minimum'],spec['rows'])
    powers = [node['power'] for node in saved['leaves'].values()]
    core.require(all(type(p) is int and -200 <= p <= 100000 for p in powers),'Invalid retained bound')
    source_sum = sum((F(2)**p for p in powers),F(0))
    refined = factor.density_factor(spec['rows'])
    ratio = F(refined,spec['rows']+1)**256
    total = source_sum*ratio
    core.require(total <= F(2)**-60,'Refined dense range has not reached the 60-bit allocation')
    result = dict(status='REFINED_DENSITY_DENSE_RANGE_CERTIFIED',instance=spec,
                  interval=[saved['minimum'],spec['rows']],coverage=coverage,
                  source_cover=source.resolve().relative_to(base.ROOT).as_posix(),source_cover_sha256=base.sha(source),
                  density_factor=refined,reference_density_factor=spec['rows']+1,
                  source_union_upper=base.encode(source_sum),exact_ratio=base.encode(ratio),union_upper=base.encode(total),
                  margin_bits=math.log2(total.denominator)-math.log2(total.numerator),full_distance_proved=False)
    if receipt:
        core.require(all(receipt[k] == v for k,v in result.items()),'Exact refinement replay mismatch')
        base.write_new(output.with_name(output.stem+'_replay.json'),
                       dict(status='REFINED_DENSITY_EXACT_REPLAY_PASSED',producer_sha256=base.sha(output),
                            source_cover_replay_sha256=base.sha(replay_path),source_sha256=search.old.provenance.sources()))
    else:
        pins = search.old.provenance.sources()
        for path in (source,replay_path,base.HERE/'POISSON_DENSITY_REFINEMENT.md',base.HERE/'DENSE_TWO_TILT.md'):
            pins[path.resolve().relative_to(base.ROOT).as_posix()] = base.sha(path)
        base.write_new(output,dict(result,source_sha256=pins))
    print(result['status'],'range',result['interval'],'margin',result['margin_bits'],flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.source,a.output,a.verify)
