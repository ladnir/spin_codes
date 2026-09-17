"""Refine Q1 tilt witnesses without changing maps or frozen certificates.

The search uses a bounded grid and refines near its best complete outer bound.
Each weight takes the minimum of independently valid outward bounds. Replay
uses linear region products at 512 bits. Existing sparse/dense contributions
are retained exactly; no new all-occupancy transfer is asserted.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

from flint import ctx
import ladder

model = ladder.model
HERE = Path(__file__).resolve().parent


def bits(value):
    return math.log2(value.denominator) - math.log2(value.numerator)


def coefficients(engine, tilt, linear=False):
    moment = engine.q1(tilt, linear=linear)
    factor = engine.length * (engine.cutoff * model.number(tilt).exp()).exp()
    return {w: model.exact(model.up(moment[w] * factor)) for w in model.base.WEIGHTS}


def minimum(left, right):
    if set(left) != set(right):
        raise ValueError('Weight coverage differs')
    return {w: min(left[w], right[w]) for w in left}


def run(output, verify=False):
    ctx.prec = 512 if verify else 256
    if not verify and output.exists():
        raise FileExistsError('Use a fresh output path')
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
    original = model.base.read(HERE / 'Q1_LOWER.json')
    full = model.base.read(HERE / 'FULL_M16_VERIFIED.json')
    row = next(r for r in original['results'] if r['instance']['message_bits'] == 65536)
    engine = ladder.Engine(16)
    assert engine.identity() == row['instance'] == full['instance']
    if saved:
        assert saved['instance'] == engine.identity()
    old = {int(w): model.exact(model.unpack(v)) for w, v in row['q1']['coefficients'].items()}
    best = dict(old)
    old_upper = model.base.decode(row['q1']['upper'])
    assert model.base.bch_bound(old)[0] == old_upper
    witnesses = []
    scores = []

    def evaluate(tilt):
        nonlocal best
        tilt = F(tilt)
        current = coefficients(engine, tilt, linear=verify)
        best = minimum(best, current)
        upper, _, _ = model.base.bch_bound(current)
        witnesses.append(str(tilt))
        score = bits(upper)
        scores.append((score, tilt))
        print('tilt', float(tilt), 'single_bits', score,
              'combined_bits', bits(model.base.bch_bound(best)[0]), flush=True)

    if saved:
        for tilt in saved['tilts']:
            evaluate(tilt)
        accepted = {int(w): model.base.decode(v) for w, v in saved['coefficients'].items()}
        assert set(accepted) == set(best)
        assert all(best[w] <= accepted[w] <= old[w] for w in best)
        best = accepted
    else:
        # Include and substantially extend the old scaled-lambda grid, 2.35..4.25.
        for j in range(-12, 13):
            evaluate(F.from_float(math.log(3 / engine.length) + j / 8))
        center = max(scores)[1]
        for j in range(-7, 8):
            if j:
                evaluate(center + F(j, 64))

    upper, factor, rest = model.base.bch_bound(best)
    old_parts = [model.base.decode(v) for v in full['component_upper']]
    assert old_parts[0] == old_upper
    assert sum(old_parts, F()) == model.base.decode(full['union_upper'])
    union = upper + sum(old_parts[1:], F())
    diagnostic = dict(old_q1_bits=bits(old_upper), q1_bits=bits(upper),
                      old_full_bits=bits(sum(old_parts, F())), full_bits=bits(union),
                      unchanged_other_occupancies_limit_bits=bits(sum(old_parts[1:], F())),
                      remaining_shell_fraction=float(rest / upper))
    result = dict(status='IMT_K16_Q1_TILT_REFINEMENT', instance=engine.identity(),
                  tilts=witnesses, coefficients={str(w): model.base.encode(v) for w, v in best.items()},
                  q1_upper=model.base.encode(upper), union_upper=model.base.encode(union),
                  retained_other_upper=model.base.encode(sum(old_parts[1:], F())),
                  diagnostic=diagnostic)
    if saved:
        for key in ('q1_upper', 'union_upper', 'retained_other_upper', 'diagnostic'):
            assert result[key] == saved[key], key
        replay = output.with_name(output.stem + '_replay.json')
        model.base.write_new(replay, dict(status='IMT_K16_Q1_TILT_512_BIT_LINEAR_REPLAY_PASSED',
                             producer_sha256=model.base.sha(output), diagnostic=diagnostic))
    else:
        result['source_sha256'] = {**ladder.candidate.sources(),
            **{p.relative_to(model.ROOT).as_posix(): model.base.sha(p)
               for p in (HERE / 'Q1_LOWER.json', HERE / 'FULL_M16_VERIFIED.json')}}
        model.base.write_new(output, result)
    print(diagnostic, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    run(args.output.resolve(), args.verify)
