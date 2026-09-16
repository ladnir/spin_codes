"""Combine the verified sparse prefix and a verified refined dense suffix."""
import argparse
import math
from pathlib import Path

import numpy as np
import parameter_range_bridge as bridge
import parameter_scoped_cover as scoped


def combine(prefix_paths,dense_paths,output):
    if output.exists():
        raise FileExistsError(output)
    bank = bridge.bank
    model,base = bank.model,bank.base
    prefix,pv,pins = bridge.union.pair(*prefix_paths,'VERIFIED_BINARY64_IMT_RANGE_BRIDGED_FULL_PARAMETER_BOUND')
    dense,dv,extra = bridge.union.pair(*dense_paths,'VERIFIED_'+scoped.STATUS)
    pins.update(extra)
    assert prefix['status'] == 'BINARY64_IMT_RANGE_BRIDGED_FULL_PARAMETER_BOUND'
    assert dense['status'] == scoped.STATUS
    assert prefix['geometry'] == pv['geometry'] == dense['geometry'] == dv['geometry'] == [128,64,20,20]
    assert prefix['inner'] == dense['inner']
    through = prefix['through']
    assert [dense['dense']['occupation_min'],dense['dense']['occupation_max']] == dv['covered_occupancies'] == [through+1,16384]
    assert abs(dv['margin_bits']-dense['dense']['margin_bits']) < 1e-7
    assert prefix['covered_occupancies'] == pv['covered_occupancies'] == [1,16384]

    def read_pair(key,status):
        record,replay,source = bridge.union.pair(*(base.grid.ROOT/p for p in key),status)
        pins.update(source)
        return record,replay
    q1,qv = read_pair(prefix['inputs']['q1'],'VERIFIED_BINARY64_IMT_Q1_GRID')
    assert q1['maps']['t64_s20'] == prefix['inner'] and qv['geometries_checked'] == 130
    first = [row for row in q1['rows'] if [row['b'],row['t'],row['s'],row['exponent']] == prefix['geometry']]
    assert len(first) == 1
    terms = [-first[0]['q1_margin_bits']*math.log(2)]
    small,sv = read_pair(prefix['inputs']['compositions'],'VERIFIED_BINARY64_IMT_SPARSE_COMPOSITIONS')
    assert small['inner'] == prefix['inner'] and small['geometry'] == sv['geometry'] == prefix['geometry']
    assert [r['q'] for r in small['results']] == [2,3,4] and sv['covered_occupancies'] == [2,4]
    terms.extend(r['log_bound'] for r in small['results'])
    next_q = 5
    for pair in prefix['inputs']['ranges']:
        record,rv = read_pair(pair,'VERIFIED_BINARY64_IMT_SPARSE_RANGE')
        assert record['inner'] == prefix['inner'] and record['geometry'] == rv['geometry'] == prefix['geometry']
        assert [record['first'],record['last']] == rv['covered_occupancies']
        assert record['first'] == next_q
        end = min(record['last'],through)
        assert end >= next_q and len(record['log_bounds']) == record['last']-record['first']+1
        terms.extend(record['log_bounds'][:end-next_q+1])
        next_q = end+1
    assert next_q == through+1
    prefix_margin = -float(np.logaddexp.reduce(terms))/math.log(2)
    terms.append(dense['dense']['log_union_upper'])
    assert all(math.isfinite(v) for v in terms)
    full = -float(np.logaddexp.reduce(terms))/math.log(2)
    pins.update(base.grid.ladder.candidate.sources())
    model.base.write_new(output,dict(status='REPLAYED_BINARY64_IMT_SCOPED_FULL_PARAMETER_BOUND',
        geometry=prefix['geometry'],inner=prefix['inner'],covered_occupancies=[1,16384],
        sparse_through=through,sparse_prefix_margin_bits=prefix_margin,
        dense_margin_bits=dense['dense']['margin_bits'],full_margin_bits=full,
        q1_margin_bits=first[0]['q1_margin_bits'],full_bound_evaluated=True,full_bound_useful=full>0,
        full_distance_proved=False,source_sha256=pins))
    print('scoped full margin',full,'sparse prefix',prefix_margin,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prefix',type=Path,nargs=2,required=True)
    p.add_argument('--dense',type=Path,nargs=2,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    combine([x.resolve() for x in a.prefix],[x.resolve() for x in a.dense],a.output.resolve())
