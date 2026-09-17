"""512-bit replay with native positive matrix products for the Q1 recurrence.

The producer's scalar recurrence is left intact. This evaluator batches each
degree-by-state layer into an Arb matrix; region products remain linear.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

from flint import arb, arb_mat, ctx
import mixing_rounds as study

model = study.model


def moments(zero, one, n, length):
    z = arb_mat([list(zero[i*n:(i+1)*n]) for i in range(n)])
    a = arb_mat([list(one[i*n:(i+1)*n]) for i in range(n)])
    current = arb_mat([[arb(j==0) for j in range(n)]])
    for count in range(1,length+1):
        inactive, active = current*z, current*a
        current = arb_mat([[(inactive[w,j] if w<count else arb(0))+
                            (active[w-1,j] if w else arb(0))
                            for j in range(n)] for w in range(count+1)])
    return [sum((current[w,j] for j in range(n)),arb(0))/math.comb(length,w)
            for w in range(length+1)]


def run(source, output):
    if output.exists():
        raise FileExistsError('Use a fresh replay output')
    ctx.prec = 512
    saved = model.base.read(source)
    model.authenticate(saved)
    assert saved['status'] == 'OUTWARD_IMT_MIXING_ROUNDS_Q1'
    single = model.independent.single
    original = single.moments
    single.moments = moments
    try:
        for row in saved['cells']:
            engine = study.Engine(row['exponent'],row['rounds'])
            assert engine.identity() == row['instance']
            best = None
            for tilt in row['tilts']:
                current = study.audit.coefficients(engine,F(tilt),linear=True)
                best = current if best is None else study.audit.minimum(best,current)
            accepted = {int(w):model.base.decode(v) for w,v in row['coefficients'].items()}
            assert set(best) == set(accepted) and all(best[w] <= accepted[w] for w in best)
            upper,_,rest = model.base.bch_bound(accepted)
            assert upper == model.base.decode(row['q1_upper'])
            assert study.audit.bits(upper) == row['q1_margin_bits']
            assert float(rest/upper) == row['remaining_shell_fraction']
            print('verified',row['exponent'],row['rounds'],row['q1_margin_bits'],flush=True)
    finally:
        single.moments = original
    model.base.write_new(output,dict(status='IMT_MIXING_ROUNDS_Q1_512_BIT_NATIVE_REPLAY_PASSED',
        producer_sha256=model.base.sha(source),cells_checked=len(saved['cells']),
        region_evaluation='linear products',coefficient_evaluation='positive native Arb matrices',
        full_distance_proved=False,source_sha256=study.ladder.candidate.sources()))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    run(args.source.resolve(),args.output.resolve())
