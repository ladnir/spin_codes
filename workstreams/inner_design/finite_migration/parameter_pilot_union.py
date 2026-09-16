"""Combine authenticated all-occupancy numerical bounds for the IMT pilot."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import parameter_no_constant as maps

HERE = Path(__file__).resolve().parent
model = maps.grid.ladder.model


def pair(producer_path,replay_path,status):
    producer,replay = (model.base.read(p) for p in (producer_path,replay_path))
    model.authenticate(producer); model.authenticate(replay)
    assert replay['status'] == status
    assert replay['producer_sha256'] == model.base.sha(producer_path)
    return producer,replay,{p.relative_to(model.ROOT).as_posix():model.base.sha(p)
                            for p in (producer_path,replay_path)}


def run(output,dense_path,dense_replay):
    assert not output.exists()
    geometry = [128,64,20,20]
    sources = {}
    q1,q1v,paths = pair(HERE/'PARAMETER_NO_CONSTANT_Q1_v1.json',HERE/'PARAMETER_NO_CONSTANT_Q1_VERIFIED_v1.json',
                        'VERIFIED_BINARY64_IMT_Q1_GRID')
    sources.update(paths)
    sparse,sv,paths = pair(HERE/'PARAMETER_NO_CONSTANT_SPARSE_v1.json',HERE/'PARAMETER_NO_CONSTANT_SPARSE_VERIFIED_v1.json',
                          'VERIFIED_BINARY64_IMT_SPARSE_COMPOSITIONS')
    sources.update(paths)
    ranges,rv,paths = pair(HERE/'PARAMETER_NO_CONSTANT_RANGE_v1.json',HERE/'PARAMETER_NO_CONSTANT_RANGE_VERIFIED_v1.json',
                          'VERIFIED_BINARY64_IMT_SPARSE_RANGE')
    sources.update(paths)
    dense,dv,paths = pair(dense_path,dense_replay,'VERIFIED_BINARY64_IMT_DENSE_COVER')
    sources.update(paths)
    record = json.loads(json.dumps(maps.inner(64,20)))
    assert q1['maps']['t64_s20'] == record
    assert q1v['geometries_checked'] == 130 and q1['geometry_complete']
    for result in (sparse,ranges,dense):
        assert result['geometry'] == geometry and result['inner'] == record
    for receipt in (sv,rv,dv):
        assert receipt['geometry'] == geometry
    first = [r for r in q1['rows'] if [r['b'],r['t'],r['s'],r['exponent']] == geometry]
    assert len(first) == 1
    assert [r['q'] for r in sparse['results']] == [2,3,4]
    assert sv['covered_occupancies'] == [2,4]
    assert [ranges['first'],ranges['last']] == rv['covered_occupancies'] == [5,64]
    assert len(ranges['log_bounds']) == 60
    assert [dense['dense']['occupation_min'],dense['dense']['occupation_max']] == dv['covered_occupancies'] == [65,16384]
    terms = [-first[0]['q1_margin_bits']*math.log(2),*[r['log_bound'] for r in sparse['results']],
             *ranges['log_bounds'],dense['dense']['log_union_upper']]
    assert all(math.isfinite(v) for v in terms)
    margin = -float(np.logaddexp.reduce(terms))/math.log(2)
    sources.update(maps.grid.ladder.candidate.sources())
    model.base.write_new(output,dict(status='REPLAYED_BINARY64_IMT_FULL_PARAMETER_BOUND',
        geometry=geometry,inner=record,covered_occupancies=[1,16384],
        full_margin_bits=margin,q1_margin_bits=first[0]['q1_margin_bits'],
        full_bound_evaluated=True,full_bound_useful=margin>0,full_distance_proved=False,
        arithmetic='nearest binary64; not an outward certificate',source_sha256=sources))
    print('full pilot numerical margin',margin,flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--dense',type=Path,required=True)
    p.add_argument('--dense-replay',type=Path,required=True)
    a = p.parse_args()
    run(a.output.resolve(),a.dense.resolve(),a.dense_replay.resolve())
