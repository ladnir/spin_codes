"""Replay the joint mixed-family IMT pilot and its exact integer partition."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import parameter_joint_mixed_dense as producer
import type_box_coverage


def verify(path,output):
    base = producer.base
    model = base.grid.ladder.model
    saved = model.base.read(path)
    model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_JOINT_MIXED_DENSE_COVER'
    cls = producer.Dense
    b,t,s,exponent = saved['geometry']
    with producer.maps.use():
        checker = cls(t,s,b,exponent,[])
        assert saved['inner'] == json.loads(json.dumps(checker.record))
        result = saved['dense']
        assert result['occupation_max'] == checker.length and result['bands'] == checker.bands
        for bank,expected in zip(checker.probability_banks,result['probability_banks'],strict=True):
            scale,p,costs = bank
            assert scale == expected['scale']
            np.testing.assert_allclose(p,expected['probabilities'],rtol=0,atol=1e-14)
            np.testing.assert_allclose(costs,expected['log_density_costs'],rtol=0,atol=1e-12)
        boxes = result['selected_boxes']
        count = type_box_coverage.check(boxes,checker.length,result['occupation_min'],len(checker.ps))
        values,errors = [],[]
        for box in boxes:
            lo,hi,w = box['lower'],box['upper'],box['witness']
            corners = base.typed.vertices(lo,hi,checker.length)
            proposal = np.array(w['proposal'])
            assert math.isfinite(w['tilt']) and w['probability_bank'] == 0
            theta = float(proposal@checker.ps)
            moment = (checker.replay_moment(theta,w) if 'family' in w
                      else checker.moment(theta,w['tilt']))
            vertices = base.typed.point_logs(corners,checker.length,b,checker.log_gammas,proposal,moment,
                                             checker.cutoff,math.exp(w['tilt']))
            value = float(max(vertices))+base.typed.lattice_log_count(lo,hi)
            assert math.isfinite(value)
            errors.append(abs(value-box['own_log_bound']))
            values.append(value)
        error = max(errors)
        assert error < 1e-7
        union = float(np.logaddexp.reduce(values))
        assert abs(union-result['log_union_upper']) < 1e-7
        assert abs(-union/math.log(2)-saved['margin_bits']) < 1e-7
        model.base.write_new(output,dict(status='VERIFIED_BINARY64_IMT_DENSE_COVER',
            producer_sha256=model.base.sha(path),geometry=saved['geometry'],
            covered_occupancies=[result['occupation_min'],checker.length],
            boxes_checked=len(boxes),integer_types_checked=str(count),
            margin_bits=-union/math.log(2),maximum_log_error=error,
            full_distance_proved=False,source_sha256=base.grid.ladder.candidate.sources()))
        print('verified joint dense cover',len(boxes),'boxes, margin',-union/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    verify(a.input.resolve(),a.output.resolve())
