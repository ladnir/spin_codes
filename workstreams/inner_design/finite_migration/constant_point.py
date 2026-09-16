"""Two all-one-count slices at a fixed dense point; no full cover claim."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

from flint import arb, ctx
from scipy.optimize import minimize_scalar
import constant_density
import fixed_input
import short_dense

model = fixed_input.model


def evaluate(checker, q, coordinate, witness, maximum, xi):
    nu = checker.density(coordinate, coordinate)[0]
    theta = F(q, checker.rows)*nu
    r, eta = (model.base.decode(witness[k]) for k in ('fixed_r', 'eta'))
    x = theta*(1-r)/(r*(1-theta))
    costs = fixed_input.scalar.tilted_costs(checker.bands, checker.ps, checker.caps, model.number(x).log())
    assert xi <= 0
    log_s = sum(((g-model.number(eta*p+(xi if i == len(costs)-1 else 0))).exp()
                 for i, (g, p) in enumerate(zip(costs, checker.ps+(F(1),)))), arb(0)).log().upper()
    scalar = fixed_input.scalar.scalar_bound(checker.rows, q, q, nu, nu, log_s, eta, model.number(x))
    # This branch bounds all-one fractions >= maximum; xi<=0 makes the
    # smallest allowed count the largest scalar penalty.
    old_log_s = sum(((g-model.number(eta*p)).exp()
                     for g, p in zip(costs, checker.ps+(F(1),))), arb(0)).log().upper()
    old_scalar = fixed_input.scalar.scalar_bound(checker.rows, q, q, nu, nu, old_log_s, eta, model.number(x))
    old = checker.fixed_input_bound(q, q, coordinate, coordinate, witness)
    return model.up(old + scalar-old_scalar+model.number(q*xi*maximum))


def run(seed, output, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_CONSTANT_COUNT_SPLIT_POINT'
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    prior = model.base.read(seed)
    model.authenticate(prior)
    q, coordinate, witness = prior['q'], F(prior['coordinate']), prior['witness']
    checker = short_dense.Checker(16)
    assert prior['instance'] == checker.engine.identity()
    nu = checker.density(coordinate, coordinate)[0]
    maximum = F(1, 32)
    old = checker.fixed_input_bound(q, q, coordinate, coordinate, witness)
    original = checker.comparison(q, q, coordinate, coordinate, ctx.prec)
    constrained = constant_density.factor(checker.rows, q, q, nu, nu, checker.ps, maximum)
    ordinary = model.up(old+256*(constrained.log()-original.log()))
    if saved:
        xi = F(saved['xi'])
    else:
        result = minimize_scalar(lambda x: float(evaluate(checker, q, coordinate, witness, maximum,
                                  F(round(16*x), 16))), bounds=(-512, 0), method='bounded')
        xi = F(round(16*result.x), 16)
    exceptional = evaluate(checker, q, coordinate, witness, maximum, xi)
    upper = model.up((ordinary.exp()+exceptional.exp()).log())
    result = dict(instance=checker.engine.identity(), q=q, coordinate=str(coordinate),
                  witness=witness, maximum=str(maximum), xi=str(xi),
                  ordinary_margin=float(-ordinary/arb(2).log()),
                  exceptional_margin=float(-exceptional/arb(2).log()),
                  margin_bits=float(-upper/arb(2).log()),
                  density_factor=float(constrained), old_density_factor=float(original),
                  log_upper=model.base.encode(model.exact(upper)))
    print({k:v for k,v in result.items() if k not in ('instance','witness','log_upper')}, flush=True)
    if saved:
        assert saved['maximum'] == str(maximum)
        assert upper <= model.number(model.base.decode(saved['log_upper']))
        model.base.write_new(output.with_name(output.stem+'_replay.json'), dict(
            status='IMT_CONSTANT_COUNT_SPLIT_POINT_512_REPLAY_PASSED',
            producer_sha256=model.base.sha(output), full_distance_proved=False))
    else:
        sources = short_dense.mixed_dense.ladder.candidate.sources()
        sources[seed.relative_to(model.ROOT).as_posix()] = model.base.sha(seed)
        model.base.write_new(output, dict(status='IMT_CONSTANT_COUNT_SPLIT_POINT',
            **result, full_distance_proved=False, source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.seed.resolve(), a.output.resolve(), a.verify)
