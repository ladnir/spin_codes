"""Reconstruct all selected moments and the integer activation-bank cover."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import parameter_activation_bank as producer
import type_box_coverage


def verify(path, output):
    if output.exists():
        raise FileExistsError(output)
    base,model = producer.base,producer.model
    saved = model.base.read(path)
    model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_ACTIVATION_BANK_DENSE_COVER'
    assert saved['geometry'] == [128,64,20,20]
    with producer.previous.maps.use():
        original = producer.previous.Dense(64,20,128,20,[])
        activation = producer.activation.Dense(64,20,128,20,[])
        assert saved['inner'] == json.loads(json.dumps(original.record))
        assert original.record == activation.record
        dense = saved['dense']
        assert dense['occupation_min'] == 65 and dense['occupation_max'] == original.length
        assert dense['bands'] == original.bands
        for bank,expected in zip(original.probability_banks,dense['probability_banks'],strict=True):
            scale,p,costs = bank
            assert scale == expected['scale']
            np.testing.assert_allclose(p,expected['probabilities'],rtol=0,atol=1e-14)
            np.testing.assert_allclose(costs,expected['log_density_costs'],rtol=0,atol=1e-12)
        boxes = dense['selected_boxes']
        count = type_box_coverage.check(boxes,original.length,65,3)
        assert str(count) == saved['integer_types_checked']
        cache,values,errors = {},[],[]
        for i,box in enumerate(boxes):
            w = box['witness']
            proposal = np.array(w['proposal'])
            assert proposal.shape == (3,) and np.all(np.isfinite(proposal))
            assert np.all(proposal > 0) and abs(float(sum(proposal))-1) < 1e-14
            assert w['probability_bank'] == 0 and math.isfinite(w['tilt'])
            theta = float(proposal@original.ps)
            if w.get('family') == 'syndrome_activation':
                key = (theta,w['tilt'])
                if key not in cache:
                    if len(activation.epoch_cache) > 4:
                        activation.epoch_cache.clear()
                    cache[key] = activation.moment(*key)
                moment = cache[key]
            elif 'family' in w:
                moment = original.replay_moment(theta,w)
            else:
                moment = original.moment(theta,w['tilt'])
            vertices = base.typed.vertices(box['lower'],box['upper'],original.length)
            value = float(max(base.typed.point_logs(vertices,original.length,128,original.log_gammas,
                proposal,moment,original.cutoff,math.exp(w['tilt']))))
            value += base.typed.lattice_log_count(box['lower'],box['upper'])
            assert math.isfinite(value)
            errors.append(abs(value-box['own_log_bound']))
            values.append(value)
            if (i+1) % 128 == 0:
                print('replayed activation boxes',i+1,flush=True)
        error = max(errors)
        assert error < 1e-7, error
        union = float(np.logaddexp.reduce(values))
        assert abs(union-dense['log_union_upper']) < 1e-7
        assert abs(-union/math.log(2)-saved['margin_bits']) < 1e-7
        sources = base.grid.ladder.candidate.sources()
        sources[path.relative_to(base.grid.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output,dict(status='VERIFIED_BINARY64_IMT_DENSE_COVER',
            producer_sha256=model.base.sha(path),geometry=saved['geometry'],
            covered_occupancies=[65,original.length],boxes_checked=len(boxes),
            integer_types_checked=str(count),margin_bits=-union/math.log(2),
            maximum_log_error=error,full_distance_proved=False,source_sha256=sources))
        print('verified activation bank margin',-union/math.log(2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    verify(args.input.resolve(),args.output.resolve())
