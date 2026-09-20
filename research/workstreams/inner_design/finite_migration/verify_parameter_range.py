"""Replay the retained IMT sparse-range witnesses in the log domain."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

import numpy as np
import parameter_sparse_range as producer


def verify(path,output):
    full = producer.full
    model = full.grid.ladder.model
    saved = model.base.read(path)
    model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_SPARSE_RANGE'
    b,t,s,exponent = saved['geometry']
    record,_,kernel,caps,low = full.prepare(t,s)
    assert saved['inner'] == json.loads(json.dumps(record))
    first,last = saved['first'],saved['last']
    rows = (1 << exponent)//(b//2)
    assert 2 <= first <= last < rows and rows % t == 0
    assert len(saved['log_bounds']) == len(saved['witnesses']) == last-first+1
    counts,_ = full.grid.outer(b)
    groups = defaultdict(list)
    for q,witness in zip(range(first,last+1),saved['witnesses']):
        assert witness is not None and witness['scale'] in (0.,0.5,1.)
        assert math.isfinite(witness['tilt'])
        groups[witness['tilt']].append((q,witness['scale']))
    replay = np.full(last-first+1,np.nan)
    for tilt,entries in groups.items():
        maximum = max(q for q,_ in entries)
        epoch = full.g.epochs(record['spectrum'],kernel,record['feedback_columns'],caps,low,
                              math.exp(tilt),maximum=min(t,maximum))
        # Independent arithmetic path: do not use the producer's positive
        # scaled convolution. Retain all coefficients in logarithmic form.
        regions = full.g.regions(epoch,t,rows,maximum)
        for scale in sorted({scale for _,scale in entries}):
            targets = {q for q,a in entries if a == scale}
            bands,roots,active,inactive = full.balanced.density_roots(counts,b,scale)
            envelope = full.balanced.Envelope(roots,active,inactive)
            current = regions.copy()
            for q in range(1,max(targets)+1):
                current = envelope.apply(current[:-1],current[1:])
                if q in targets:
                    replay[q-first] = (math.log(math.comb(rows,q))+q*math.log(len(bands))
                        +(b*rows//10)*math.exp(tilt)+full.g.terminal(current[0],b))
        print('replayed sparse range tilt',tilt,flush=True)
    assert np.all(np.isfinite(replay))
    error = float(max(abs(replay-np.array(saved['log_bounds']))))
    assert error < 1e-8
    margin = -float(np.logaddexp.reduce(replay))/math.log(2)
    assert abs(margin-saved['margin_bits']) < 1e-8
    model.base.write_new(output,dict(status='VERIFIED_BINARY64_IMT_SPARSE_RANGE',
        producer_sha256=model.base.sha(path),geometry=saved['geometry'],
        covered_occupancies=[first,last],occupancies_checked=last-first+1,
        margin_bits=margin,maximum_log_error=error,arithmetic='nearest binary64 log-domain replay',
        full_distance_proved=False,source_sha256=full.grid.ladder.candidate.sources()))
    print('verified sparse range',first,last,'margin',margin,'error',error,flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    verify(a.input.resolve(),a.output.resolve())
