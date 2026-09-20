"""Authenticate spectrum-only setup events and check event containment."""
import json

import activation_occupation as markov
import random_spectrum_variance as variance
import run_complete_q1_grid as grid

MARKOV_METHOD='simultaneous Markov shell caps; one reused full-rank subspace'
VARIANCE_METHOD='simultaneous minimum of Markov and exact-variance Chebyshev shell caps; one reused full-rank subspace'


def implies(stronger,weaker):
    if (stronger['block_bits'],stronger['dimension'])!=(weaker['block_bits'],weaker['dimension']):return False
    return all(stronger['counts'].get(w,0)<=weaker['counts'].get(w,0) for w in range(1,stronger['block_bits']+1))


def load_registered_events():
    catalog=json.loads((grid.HERE/'catalog.json').read_text())
    result={}
    for spec in catalog['activation_sources']:
        receipt=json.loads((grid.HERE/spec['manifest']).read_text())
        for name,digest in receipt['source_sha256'].items():
            if not name.endswith('_caps.json') or digest in result:continue
            path=grid.pilot.ROOT/name
            if grid.pilot.sha(path)!=digest:raise ValueError('spectrum event hash changed')
            payload=json.loads(path.read_text())
            b,d,h=payload['block_bits'],payload['dimension'],payload['setup_failure_bits']
            if payload['method']==MARKOV_METHOD:expected=markov.random_spectrum_caps(b,d,h)
            elif payload['method']==VARIANCE_METHOD:expected=variance.caps(b,d,h)
            else:raise ValueError('unknown setup event semantics')
            counts={int(w):int(n) for w,n in payload['counts'].items()}
            if counts!=expected:raise ValueError('spectrum event does not match its declared probability bound')
            result[digest]=dict(block_bits=b,dimension=d,setup_failure_bits=h,counts=counts,source=name)
    return result


def containment(events):
    return {(a,b):implies(first,second) for a,first in events.items() for b,second in events.items()}
