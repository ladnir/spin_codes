"""Replayable comparison of all-one fraction cutoffs at the difficult point."""
import argparse
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx
import constant_point as point

model = point.model
DENOMINATORS = (32, 64, 128, 256, 512, 1024, 2048)


def run(seed, output, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_CONSTANT_THRESHOLD_POINTS'
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    prior = model.base.read(seed)
    model.authenticate(prior)
    checker = point.short_dense.Checker(16)
    assert prior['instance'] == checker.engine.identity()
    q, coordinate, witness = prior['q'], F(prior['coordinate']), prior['witness']
    nu = checker.density(coordinate, coordinate)[0]
    old = checker.fixed_input_bound(q, q, coordinate, coordinate, witness)
    original = checker.comparison(q, q, coordinate, coordinate, ctx.prec)
    results = []
    for i, denominator in enumerate(DENOMINATORS):
        maximum, xi = F(1, denominator), F(-3335, 16)
        ratio = point.constant_density.factor(checker.rows, q, q, nu, nu, checker.ps, maximum)
        ordinary = model.up(old+256*(ratio.log()-original.log()))
        exceptional = point.evaluate(checker, q, coordinate, witness, maximum, xi)
        upper = model.up((ordinary.exp()+exceptional.exp()).log())
        record = dict(maximum=str(maximum), xi=str(xi),
            ordinary_margin=float(-ordinary/arb(2).log()),
            exceptional_margin=float(-exceptional/arb(2).log()),
            margin_bits=float(-upper/arb(2).log()), log_upper=model.base.encode(model.exact(upper)))
        if saved:
            previous = saved['results'][i]
            assert previous['maximum'] == str(maximum) and previous['xi'] == str(xi)
            assert upper <= model.number(model.base.decode(previous['log_upper']))
        results.append(record)
        print('threshold', ctx.prec, denominator, record['margin_bits'], flush=True)
    if saved:
        assert len(saved['results']) == len(results)
        model.base.write_new(output.with_name(output.stem+'_replay.json'), dict(
            status='IMT_CONSTANT_THRESHOLD_POINTS_512_REPLAY_PASSED',
            producer_sha256=model.base.sha(output), full_distance_proved=False))
    else:
        sources = point.short_dense.mixed_dense.ladder.candidate.sources()
        sources[seed.relative_to(model.ROOT).as_posix()] = model.base.sha(seed)
        model.base.write_new(output, dict(status='IMT_CONSTANT_THRESHOLD_POINTS',
            instance=checker.engine.identity(), q=q, coordinate=str(coordinate),
            witness=witness, results=results, full_distance_proved=False, source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.seed.resolve(), a.output.resolve(), a.verify)
