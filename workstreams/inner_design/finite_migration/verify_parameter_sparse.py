"""Replay every retained small-Q composition witness for the IMT pilot."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

import numpy as np
import parameter_sparse_balanced as balanced


def verify(path,output):
    full,prior = balanced.full,balanced.prior
    model = full.grid.ladder.model
    saved = model.base.read(path)
    model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_BALANCED_COMPOSITIONS'
    b,t,s,exponent = saved['geometry']
    record,bs,kernel,caps,low = full.prepare(t,s)
    assert saved['inner'] == json.loads(json.dumps(record))
    counts,_ = full.grid.outer(b)
    rows = (1 << exponent)//(b//2)
    maximum = max(r['q'] for r in saved['results'])
    assert [r['q'] for r in saved['results']] == list(range(2,maximum+1))
    region_cache = {}
    checked = 0
    for result in saved['results']:
        problem = prior.composition.SparseComposition(counts,b,result['q'],tail_bands=6,singleton_prefix=3)
        assert problem.bands == result['bands'] and problem.indices.tolist() == result['indices']
        assert len(result['component_log_bounds']) == len(result['witnesses']) == len(problem.indices)
        groups = defaultdict(list)
        for i,w in enumerate(result['witnesses']):
            assert w is not None and w['scale'] in (0.,0.5,1.)
            groups[w['tilt'],w['scale']].append(i)
        for (tilt,scale),indices in groups.items():
            if tilt not in region_cache:
                epoch = full.g.epochs(record['spectrum'],kernel,record['feedback_columns'],caps,low,
                                      math.exp(tilt),maximum=min(t,maximum))
                region_cache[tilt] = full.g.regions(epoch,t,rows,maximum)
            regions = region_cache[tilt]
            probabilities,cost = balanced.bank(problem,scale)
            probabilities,cost = probabilities[indices],cost[indices]
            matrices = np.full((len(indices),regions.shape[1],regions.shape[2]),-np.inf)
            for j in range(problem.occupation+1):
                np.logaddexp(matrices,probabilities[:,j,None,None]+regions[j],out=matrices)
            values = np.minimum(problem.trivial[indices],prior.terminal(matrices,b)+cost+(b*rows//10)*math.exp(tilt))
            expected = np.array(result['component_log_bounds'])[indices]
            assert np.all(np.isfinite(values)) and np.max(np.abs(values-expected)) < 1e-8
            checked += len(indices)
        total = problem.aggregate(np.array(result['component_log_bounds']),rows)
        assert abs(total-result['log_bound']) < 1e-8
        assert abs(-total/math.log(2)-result['margin_bits']) < 1e-8
    model.base.write_new(output,dict(status='VERIFIED_BINARY64_IMT_SPARSE_COMPOSITIONS',
        producer_sha256=model.base.sha(path),geometry=saved['geometry'],covered_occupancies=[2,maximum],
        compositions_checked=checked,full_distance_proved=False,source_sha256=full.grid.ladder.candidate.sources()))
    print('verified',checked,'compositions at Q2..'+str(maximum),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    verify(a.input.resolve(),a.output.resolve())
