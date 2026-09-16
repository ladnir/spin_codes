"""Bounded candidate-specific Q1, sparse, and dense checks, with 512-bit replay."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from flint import arb, ctx
import candidate
import sparse_bch
model = candidate.model


def q1(engine, old=None):
    tilts = old['tilts'] if old else [str(F.from_float(math.log(a/(100*engine.length)))) for a in (235,264,295,332,376,425)]
    best = {w: arb(engine.length) for w in model.base.WEIGHTS}
    for z in tilts:
        lam = model.number(F(z)).exp()
        co = engine.q1(F(z), linear=bool(old))
        factor = engine.length*(engine.cutoff*lam).exp()
        for w in best:
            best[w] = min(best[w], model.up(co[w]*factor))
    if old:
        bounds = old['coefficients']
        assert set(bounds) == set(map(str, best))
        assert all(v <= model.unpack(bounds[str(w)]) for w, v in best.items())
    else:
        bounds = {str(w): model.pack(v) for w, v in best.items()}
    upper, _, _ = model.base.bch_bound({int(w): model.exact(model.unpack(v)) for w, v in bounds.items()})
    result = dict(tilts=tilts, coefficients=bounds, upper=model.base.encode(upper),
                  margin_bits=math.log2(upper.denominator)-math.log2(upper.numerator))
    if old:
        assert result == old
    return result


def sparse(engine, q, old=None):
    if old:
        trials = [(F(old['tilt']), list(map(model.base.decode, old['probabilities'])))]
    else:
        prediction = (-76+6*math.log2(q)-10*(engine.exponent-20)*math.log(2))/10
        trials = [(F.from_float(prediction+d), None) for d in (-.5, 0, .5)]
    best = None
    for tilt, ps in trials:
        region = engine.region(tilt, q)
        ps = sparse_bch.choose(engine, region, q) if ps is None else ps
        value = model.up(engine.adaptive(region, ps, q)*(engine.cutoff*model.number(tilt).exp()).exp())
        if best is None or value < best[0]:
            best = value, tilt, ps
    value, tilt, ps = best
    if old:
        assert value <= model.unpack(old['upper'])
        bound = old['upper']
    else:
        bound = model.pack(value)
    return dict(q=q, tilt=str(tilt), probabilities=list(map(model.base.encode, ps)), upper=bound,
                margin_bits=float(-model.unpack(bound).log()/arb(2).log()))


def dense_points(checker, old=None):
    if old:
        jobs = old
    else:
        prior = model.base.read(model.HERE/'DENSE_RETUNED_POINTS.json')['points']
        jobs = [dict(q=r['occupation'], v=r['coordinate'], witness=r['witness'])
                for r in prior if r['message_exponent'] == checker.engine.exponent]
        # Add the previously missed intermediate-density bottleneck at every length.
        gap = model.base.read(model.HERE/'DENSE_GAP_POINT.json')
        q = max(2, round(gap['occupation']*checker.rows/8192))
        jobs += [dict(q=q, v=str(v), witness=gap['witness']) for v in (F(1,8), F(21,128), F(1,4))]
    results = []
    for row in jobs:
        q, v = row['q'], F(row['v'])
        witness = row['witness']
        value = checker.bound(q, q, v, v, witness)
        if not old and value > -80*arb(2).log():
            proposal = checker.witness(q, q, v, v)
            other = checker.bound(q, q, v, v, proposal)
            if other < value:
                value, witness = other, proposal
        if old:
            assert model.up(value.exp()) <= model.unpack(row['upper'])
            bound = row['upper']
        else:
            bound = model.pack(model.up(value.exp()))
        results.append(dict(q=q, v=str(v), witness=witness, upper=bound,
                            margin_bits=float(-model.unpack(bound).log()/arb(2).log())))
        print('dense', checker.engine.exponent, q, v, results[-1]['margin_bits'], flush=True)
    return results


def run(output, name, exponents, mode, occupations, verify):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        name, exponents, mode, occupations = (saved[k] for k in ('candidate','exponents','mode','occupations'))
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    results = []
    for i, exponent in enumerate(exponents):
        engine = candidate.Engine(exponent, name)
        old = saved['results'][i] if saved else None
        result = dict(instance=engine.identity(), audit=engine.audit())
        if old:
            assert result['instance'] == old['instance'] and result['audit'] == old['audit']
        if mode in ('q1','all'):
            result['q1'] = q1(engine, old['q1'] if old else None)
            print('Q1', name, exponent, result['q1']['margin_bits'], flush=True)
        if mode in ('sparse','all'):
            result['sparse'] = []
            for j, q in enumerate(occupations):
                row = sparse(engine, q, old['sparse'][j] if old else None)
                result['sparse'].append(row)
                print('sparse', name, exponent, q, row['margin_bits'], flush=True)
            engine.region.cache_clear()
            engine.epoch.cache_clear()
        if mode in ('dense','all'):
            checker = candidate.Checker(exponent, name)
            result['dense'] = dense_points(checker, old['dense'] if old else None)
        results.append(result)
    if saved:
        model.base.write_new(output.with_name(output.stem+'_replay.json'), dict(
            status='CANDIDATE_SCREEN_512_BIT_REPLAY_PASSED', producer_sha256=model.base.sha(output), full_distance_proved=False))
    else:
        dependencies = candidate.sources()
        for filename in ('DENSE_RETUNED_POINTS.json','DENSE_GAP_POINT.json'):
            path = model.HERE/filename
            dependencies[path.relative_to(model.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output, dict(status='OUTWARD_CANDIDATE_SCREEN_NOT_FULL_CERTIFICATE',
            candidate=name, exponents=exponents, mode=mode, occupations=occupations,
            results=results, source_sha256=dependencies, full_distance_proved=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--candidate', choices=('weight5_seed0','weight5_seed1','balanced'), default='weight5_seed0')
    p.add_argument('--m', type=int, nargs='+', choices=(16,18,20), default=[16,18,20])
    p.add_argument('--mode', choices=('q1','sparse','dense','all'), default='all')
    p.add_argument('--occupations', type=int, nargs='+', default=[2,4,16,64])
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    assert a.occupations and all(2 <= q <= 512 for q in a.occupations)
    run(a.output.resolve(), a.candidate, a.m, a.mode, a.occupations, a.verify)
