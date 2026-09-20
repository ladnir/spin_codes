"""Split low-reference assignments from those with at least one high band."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

from flint import arb,ctx
import constant_density
import zero_constant_dense

model = zero_constant_dense.model


def bound(checker,q,coordinate,witness,cutoff,xi=None):
    nu = checker.density(coordinate,coordinate)[0]
    theta = F(q,checker.rows)*nu
    r,eta = (model.base.decode(witness[k]) for k in ('fixed_r','eta'))
    x = theta*(1-r)/(r*(1-theta))
    scalar_module = zero_constant_dense.parent.fixed_input.scalar
    costs = scalar_module.tilted_costs(checker.bands,checker.ps,checker.caps,model.number(x).log())
    low,high = arb(0),arb(0)
    for cost,p in zip(costs,checker.ps+(F(1),)):
        term = (cost-model.number(eta*p)).exp()
        if p <= cutoff:
            low += term
        else:
            high += term
    assert low > 0 and high > 0
    if xi is None:
        # At a point, minimizing Q log(low+high exp(-xi))+xi is
        # equivalent to giving the high-band class tilted mass 1/Q.
        optimum = float((high*(q-1)/low).log())
        xi = F(min(0,round(64*optimum)),64)
    assert xi <= 0
    def scalar(total):
        return scalar_module.scalar_bound(checker.rows,q,q,nu,nu,total.log().upper(),eta,model.number(x))
    original_scalar = scalar(low+high)
    old = checker.fixed_input_bound(q,q,coordinate,coordinate,witness)
    exceptional = model.up(old+scalar(low+high*(-model.number(xi)).exp())-original_scalar+model.number(xi))
    ps = tuple(p for p in checker.ps if p <= cutoff)
    if nu > max(ps):
        return exceptional,xi,None,exceptional
    ratio = constant_density.factor(checker.rows,q,q,nu,nu,ps,F(0))
    original = checker.comparison(q,q,coordinate,coordinate,ctx.prec)
    ordinary = model.up(old+scalar(low)-original_scalar+256*(ratio.log()-original.log()))
    total = model.up((ordinary.exp()+exceptional.exp()).log())
    return total,xi,ordinary,exceptional


def run(seed,output,verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_HIGH_BAND_SPLIT_POINTS'
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    prior = model.base.read(seed)
    model.authenticate(prior)
    checker = zero_constant_dense.Checker(16)
    assert prior['instance'] == checker.engine.identity()
    records = []
    for i,point in enumerate(prior['points']):
        q,coordinate,witness = point['q'],F(point['coordinate']),point['witness']
        previous = saved['points'][i] if saved else None
        cutoff = F(1,2)
        upper,xi,ordinary,exceptional = bound(checker,q,coordinate,witness,cutoff,
                                             F(previous['xi']) if previous else None)
        if previous:
            assert previous['q'] == q and previous['coordinate'] == str(coordinate)
            assert previous['cutoff'] == str(cutoff) and previous['witness'] == witness
            assert upper <= model.number(model.base.decode(previous['log_upper']))
        record = dict(q=q,coordinate=str(coordinate),witness=witness,cutoff=str(cutoff),xi=str(xi),
            margin_bits=float(-upper/arb(2).log()),old_margin=point['margin_bits'],
            ordinary_margin=float(-ordinary/arb(2).log()) if ordinary is not None else None,
            exceptional_margin=float(-exceptional/arb(2).log()),log_upper=model.base.encode(model.exact(upper)))
        records.append(record)
        print('high-band split',ctx.prec,q,str(coordinate),record['margin_bits'],
              'ordinary',record['ordinary_margin'],'exceptional',record['exceptional_margin'],flush=True)
    if saved:
        assert len(saved['points']) == len(records)
        model.base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='IMT_HIGH_BAND_SPLIT_POINTS_512_REPLAY_PASSED',
            producer_sha256=model.base.sha(output),full_distance_proved=False))
    else:
        sources = zero_constant_dense.parent.short_dense.mixed_dense.ladder.candidate.sources()
        sources[seed.relative_to(model.ROOT).as_posix()] = model.base.sha(seed)
        model.base.write_new(output,dict(status='IMT_HIGH_BAND_SPLIT_POINTS',points=records,
            instance=checker.engine.identity(),full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.seed.resolve(),a.output.resolve(),a.verify)
