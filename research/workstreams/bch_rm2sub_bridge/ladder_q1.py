"""Size-specific Q1 coefficients and exact BCH aggregation through K=2^24.

This certifies Q1 only. Higher occupancies require separate coverage.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

from flint import arb, ctx

import ladder_instance as identity
import certificate_search_backend as frozen
import bch256_dense_types as dependencies

core, base, sparse, refresh = frozen.core, frozen.base, frozen.sparse, frozen.refresh


def compute(spec, witnesses, precision):
    core.require(spec == identity.instance(spec['message_exponent']), 'Wrong ladder instance')
    core.require(witnesses and all(type(a) is int and 1 <= a <= 10000 for a in witnesses), 'Invalid Q1 tilt')
    ctx.prec = precision
    t, s, spectrum, _ = core.inputs.load(spec['configuration'])
    rows, cutoff = spec['rows'], spec['cutoff']
    best = {w: F(rows) for w in base.WEIGHTS}
    for a in witnesses:
        lam = arb(a)/(100*rows)
        reg = refresh.region(*refresh.epoch(t, s, spectrum, (-lam).exp(), arb), rows//t, arb, linear=precision >= 512)
        coefficients = refresh.coefficients(*reg, 256, arb)
        factor = rows*(cutoff*lam).exp()
        for w in best:
            best[w] = min(best[w], sparse.rational((coefficients[w]*factor).upper()))
        print('m', spec['message_exponent'], 'Q1 precision', precision, 'scaled tilt', a, flush=True)
    return best


def run(output, m=20, verify=False):
    saved = base.read(output) if verify else None
    if saved is not None:
        core.authenticate(saved, base.ROOT)
        core.require(saved['status'] == 'LADDER_Q1_OUTWARD_PRODUCER', 'Wrong producer')
        m, witnesses = saved['instance']['message_exponent'], saved['scaled_tilts']
    else:
        core.require(not output.exists(), 'Use fresh output')
        witnesses = [264, 295, 332, 376]
    spec = identity.instance(m)
    if saved is not None:
        core.require(spec == saved['instance'], 'Wrong replay instance')
    best = compute(spec, witnesses, 512 if verify else 256)
    if saved is not None:
        recorded = {int(w): base.decode(v) for w, v in saved['coefficient_upper'].items()}
        core.require(set(best) == set(recorded) and all(best[w] <= recorded[w] for w in best), 'Coefficient replay failed')
        best = recorded
    upper, factor, rest = base.bch_bound(best)
    if saved is not None:
        core.require((upper, factor, rest) == tuple(base.decode(saved[k]) for k in ('upper', 'factor', 'rest')), 'BCH aggregation mismatch')
        result = dict(status='LADDER_Q1_512_BIT_LINEAR_REPLAY_PASSED', producer_sha256=base.sha(output),
            instance=spec, coefficients_checked=len(best), source_sha256=dependencies.sources())
        base.write_new(output.with_name(output.stem+'_replay.json'), result)
    else:
        result = dict(status='LADDER_Q1_OUTWARD_PRODUCER', instance=spec, scaled_tilts=witnesses,
            coefficient_upper={str(w): base.encode(v) for w, v in best.items()}, upper=base.encode(upper),
            factor=base.encode(factor), rest=base.encode(rest), source_sha256=dependencies.sources(), full_distance_proved=False)
        base.write_new(output, result)
    print(result['status'], 'Q1 ONLY margin', math.log2(upper.denominator)-math.log2(upper.numerator), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, default=20)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output, a.m, a.verify)
