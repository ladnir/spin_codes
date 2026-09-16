"""Join replayed sparse occupancies to a clipped, replayed dense partition.

All values are nearest binary64 diagnostics. Source authentication, exact
integer coverage, and arithmetic replay are separate checks; none upgrades
the result to an outward certificate.
"""
import argparse
import copy
import json
import math
from pathlib import Path

import numpy as np
import parameter_activation_bank as bank
import parameter_pilot_union as union
import type_box_coverage

HERE = Path(__file__).resolve().parent


def clip(box,total,minimum):
    lo,hi = list(box['lower']),list(box['upper'])
    hi[0] = min(hi[0],total-minimum)
    if any(a>b for a,b in zip(lo,hi)):
        return None
    corners = bank.base.typed.vertices(lo,hi,total)
    if not len(corners):
        return None
    return dict(lower=corners.min(axis=0).tolist(),upper=corners.max(axis=0).tolist(),
                witness=copy.deepcopy(box['witness']))


class Moments:
    def __init__(self):
        self.original = bank.previous.Dense(64,20,128,20,[])
        self.activation = bank.activation.Dense(64,20,128,20,[])
        self.cache = {}

    def bound(self,box):
        checker = self.original
        w = box['witness']
        proposal = np.array(w['proposal'])
        assert proposal.shape == (3,) and np.all(np.isfinite(proposal)) and np.all(proposal>0)
        assert abs(float(sum(proposal))-1) < 1e-14
        assert w['probability_bank'] == 0 and math.isfinite(w['tilt'])
        theta = float(proposal@checker.ps)
        key = (w.get('family'),theta,w['tilt'])
        if key not in self.cache:
            if key[0] == 'syndrome_activation':
                if len(self.activation.epoch_cache)>4:
                    self.activation.epoch_cache.clear()
                value = self.activation.moment(theta,w['tilt'])
            elif key[0] is None:
                value = checker.moment(theta,w['tilt'])
            else:
                value = checker.replay_moment(theta,w)
            self.cache[key] = value
        corners = bank.base.typed.vertices(box['lower'],box['upper'],checker.length)
        value = float(max(bank.base.typed.point_logs(corners,checker.length,128,checker.log_gammas,
            proposal,self.cache[key],checker.cutoff,math.exp(w['tilt']))))
        value += bank.base.typed.lattice_log_count(box['lower'],box['upper'])
        assert math.isfinite(value)
        return value


def compute(inputs,through):
    model,base = bank.model,bank.base
    geometry = [128,64,20,20]
    sources = {}
    def pair(paths,status):
        result,replay,pins = union.pair(*(base.grid.ROOT/p for p in paths),status)
        sources.update(pins)
        return result,replay
    dense,dv = pair(inputs['dense'],'VERIFIED_BINARY64_IMT_DENSE_COVER')
    assert dense['status'] == 'BINARY64_IMT_ACTIVATION_BANK_DENSE_COVER'
    assert dense['geometry'] == dv['geometry'] == geometry
    q1,qv = pair(inputs['q1'],'VERIFIED_BINARY64_IMT_Q1_GRID')
    sparse,sv = pair(inputs['compositions'],'VERIFIED_BINARY64_IMT_SPARSE_COMPOSITIONS')
    all_ranges = [pair(paths,'VERIFIED_BINARY64_IMT_SPARSE_RANGE') for paths in inputs['ranges']]
    record = json.loads(json.dumps(bank.previous.maps.inner(64,20)))
    assert dense['inner'] == sparse['inner'] == q1['maps']['t64_s20'] == record
    assert sparse['geometry'] == sv['geometry'] == geometry
    assert qv['geometries_checked'] == 130 and q1['geometry_complete']
    assert [r['q'] for r in sparse['results']] == [2,3,4]
    assert sv['covered_occupancies'] == [2,4]
    first = [r for r in q1['rows'] if [r['b'],r['t'],r['s'],r['exponent']] == geometry]
    assert len(first) == 1
    terms = [-first[0]['q1_margin_bits']*math.log(2),*[r['log_bound'] for r in sparse['results']]]
    next_q = 5
    sparse_components = []
    for saved,replay in all_ranges:
        assert saved['status'] == 'BINARY64_IMT_SPARSE_RANGE'
        assert saved['geometry'] == replay['geometry'] == geometry and saved['inner'] == record
        assert [saved['first'],saved['last']] == replay['covered_occupancies']
        assert len(saved['log_bounds']) == saved['last']-saved['first']+1
        assert saved['first'] == next_q
        last = min(saved['last'],through)
        assert last >= saved['first']
        selected = saved['log_bounds'][:last-saved['first']+1]
        assert all(math.isfinite(x) for x in selected)
        terms.extend(selected)
        sparse_components.append(dict(first=saved['first'],last=last,
            margin_bits=-float(np.logaddexp.reduce(selected))/math.log(2)))
        next_q = last+1
    assert next_q == through+1
    assert 65 <= through < 16384
    with bank.previous.maps.use():
        moments = Moments()
        checker = moments.original
        assert dense['dense']['bands'] == checker.bands
        for expected,actual in zip(dense['dense']['probability_banks'],checker.probability_banks,strict=True):
            scale,p,costs = actual
            assert expected['scale'] == scale
            np.testing.assert_allclose(expected['probabilities'],p,atol=1e-14,rtol=0)
            np.testing.assert_allclose(expected['log_density_costs'],costs,atol=1e-12,rtol=0)
        assert dv['covered_occupancies'] == [65,checker.length]
        parents = dense['dense']['selected_boxes']
        type_box_coverage.check(parents,checker.length,65,3)
        boxes = [child for parent in parents if (child:=clip(parent,checker.length,through+1)) is not None]
        count = type_box_coverage.check(boxes,checker.length,through+1,3)
        for box in boxes:
            box['own_log_bound'] = moments.bound(box)
        dense_log = float(np.logaddexp.reduce([b['own_log_bound'] for b in boxes]))
        terms.append(dense_log)
    margin = -float(np.logaddexp.reduce(terms))/math.log(2)
    sources.update(base.grid.ladder.candidate.sources())
    return dict(status='BINARY64_IMT_RANGE_BRIDGED_FULL_PARAMETER_BOUND',inputs=inputs,
        through=through,geometry=geometry,inner=record,covered_occupancies=[1,16384],
        sparse_components=sparse_components,q1_margin_bits=first[0]['q1_margin_bits'],
        dense=dict(occupation_min=through+1,occupation_max=16384,selected_boxes=boxes,
                   integer_types_checked=str(count),log_union_upper=dense_log,margin_bits=-dense_log/math.log(2)),
        full_margin_bits=margin,full_bound_evaluated=True,full_bound_useful=margin>0,
        full_distance_proved=False,source_sha256=sources)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--dense',type=Path,nargs=2)
    parser.add_argument('--range',dest='ranges',type=Path,nargs=2,action='append')
    parser.add_argument('--through',type=int)
    parser.add_argument('--verify',type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    def relative(path):
        return path.resolve().relative_to(bank.base.grid.ROOT).as_posix()
    if args.verify:
        if args.dense or args.ranges or args.through is not None:
            parser.error('--verify uses the recorded inputs, not replacement arguments')
        old = bank.model.base.read(args.verify.resolve())
        bank.model.authenticate(old)
        assert old['status'] == 'BINARY64_IMT_RANGE_BRIDGED_FULL_PARAMETER_BOUND'
        result = compute(old['inputs'],old['through'])
        assert result['inner'] == old['inner'] and result['covered_occupancies'] == old['covered_occupancies']
        assert result['dense']['integer_types_checked'] == old['dense']['integer_types_checked']
        error = abs(result['full_margin_bits']-old['full_margin_bits'])
        for a,b in zip(result['dense']['selected_boxes'],old['dense']['selected_boxes'],strict=True):
            assert (a['lower'],a['upper'],a['witness']) == (b['lower'],b['upper'],b['witness'])
            error = max(error,abs(a['own_log_bound']-b['own_log_bound']))
        assert error < 1e-7,error
        result = dict(status='VERIFIED_BINARY64_IMT_RANGE_BRIDGED_FULL_PARAMETER_BOUND',
            producer_sha256=bank.model.base.sha(args.verify.resolve()),geometry=result['geometry'],
            covered_occupancies=result['covered_occupancies'],full_margin_bits=result['full_margin_bits'],
            maximum_log_error=error,full_distance_proved=False,source_sha256=result['source_sha256'])
    else:
        if not args.dense or not args.ranges or args.through is None:
            parser.error('production requires --dense, --range, and --through')
        inputs = dict(dense=list(map(relative,args.dense)),
            q1=[relative(HERE/name) for name in ('PARAMETER_NO_CONSTANT_Q1_v1.json','PARAMETER_NO_CONSTANT_Q1_VERIFIED_v1.json')],
            compositions=[relative(HERE/name) for name in ('PARAMETER_NO_CONSTANT_SPARSE_v1.json','PARAMETER_NO_CONSTANT_SPARSE_VERIFIED_v1.json')],
            ranges=[[relative(HERE/name) for name in ('PARAMETER_NO_CONSTANT_RANGE_v1.json','PARAMETER_NO_CONSTANT_RANGE_VERIFIED_v1.json')],
                    *[list(map(relative,pair)) for pair in args.ranges]])
        result = compute(inputs,args.through)
    bank.model.base.write_new(args.output.resolve(),result)
    print(result['status'],'full margin',result['full_margin_bits'],flush=True)
