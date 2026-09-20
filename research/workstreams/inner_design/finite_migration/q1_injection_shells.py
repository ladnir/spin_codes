"""Retain arbitrary shell-supported mass after a singleton activates IMT.

Only the Q1 analysis changes. Fixed maps, mixer distribution, and encoding do
not. New coordinates H_v hold arbitrary mass on wt(Aq)=v, not uniform mass.
They retain their shell on a lazy zero-input update, but fall back to the
original arbitrary-state coordinate after a lazy one-input update.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

from flint import ctx
import q1_slack as audit

ladder = audit.ladder
model = ladder.model
HERE = Path(__file__).resolve().parent


def extend(zero, one, spectrum, injection_weights, cancellation_weights, z, t):
    """Positive extension; also accepts exact rational entries for toy tests."""
    levels = sorted(spectrum)
    old_n = len(levels) + 2
    n = old_n + len(levels)
    assert len(zero) == len(one) == old_n * old_n
    assert len(injection_weights) == len(cancellation_weights) == t
    assert set(injection_weights) <= set(levels)
    m = sum(spectrum.values())
    matrices = [[z * 0 for _ in range(n*n)] for _ in range(2)]
    for source, dest in zip((zero, one), matrices):
        for i in range(old_n):
            for j in range(old_n):
                dest[i*n+j] = source[i*old_n+j]
    rz, ra = matrices
    # Nonzero, distinct feedback columns make singleton activation exact.
    assert one[0] == 0
    ra[1] = z * 0
    counts = Counter(injection_weights)
    for i, v in enumerate(levels):
        h = old_n + i
        ra[h] = z * counts[v] / t
        d0 = z**v
        d1 = (v * z**(v-1) + (t-v) * z**(v+1)) / t
        rz[h*n+h] = d0 / 2
        # At most one singleton input cancels a fixed state. Its emission
        # weight is audited directly for each feedback column in this shell.
        cancel = max((z**w for iw, w in zip(injection_weights, cancellation_weights)
                      if iw == v), default=z*0) / t
        ra[h*n] = cancel / 2 + d1 / (2*m)
        ra[h*n+1] = d1 / 2
        for j, w in enumerate(levels):
            rz[h*n+j+2] = d0 * spectrum[w] / (2*m)
            ra[h*n+j+2] = d1 * spectrum[w] / (2*m)
    return tuple(rz), tuple(ra), n


class Engine(ladder.Engine):
    def __init__(self, exponent=16):
        super().__init__(exponent)
        self.original = ladder.Engine(exponent)
        self.n += len(self.levels)
        assert len(set(self.columns)) == len(self.columns) and all(self.columns)
        image = lambda q: sum(((q & a).bit_count() & 1) << i
                              for i, a in enumerate(self.a_columns))
        images = [image(q) for q in self.columns]
        self.injection_weights = [v.bit_count() for v in images]
        self.cancellation_weights = [(v ^ (1 << i)).bit_count() for i, v in enumerate(images)]

    @lru_cache(maxsize=64)
    def epoch(self, log_lam):
        zero, one = self.original.epoch(log_lam)[:2]
        z = (-model.number(log_lam).exp()).exp()
        rz, ra, n = extend(zero, one, self.spectrum, self.injection_weights,
                           self.cancellation_weights, z, len(self.columns))
        assert n == self.n
        return tuple(map(model.up, rz)), tuple(map(model.up, ra))


def run(output, verify=False):
    ctx.prec = 512 if verify else 256
    if not verify and output.exists():
        raise FileExistsError('Use a fresh output path')
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
    original = model.base.read(HERE / 'Q1_LOWER.json')
    row = next(r for r in original['results'] if r['instance']['message_bits'] == 65536)
    full = model.base.read(HERE / 'FULL_M16_VERIFIED.json')
    engine = Engine()
    assert engine.identity() == row['instance'] == full['instance']
    best = {int(w): model.exact(model.unpack(v)) for w, v in row['q1']['coefficients'].items()}
    tilts = saved['tilts'] if saved else row['q1']['tilts']
    for tilt in tilts:
        current = audit.coefficients(engine, F(tilt), linear=verify)
        best = audit.minimum(best, current)
        print('tilt', float(F(tilt)), 'q1_bits', audit.bits(model.base.bch_bound(best)[0]), flush=True)
    if saved:
        accepted = {int(w): model.base.decode(v) for w, v in saved['coefficients'].items()}
        assert set(best) == set(accepted) and all(best[w] <= accepted[w] for w in best)
        best = accepted
    upper, _, rest = model.base.bch_bound(best)
    parts = [model.base.decode(v) for v in full['component_upper']]
    assert parts[0] == model.base.decode(row['q1']['upper'])
    assert sum(parts, F()) == model.base.decode(full['union_upper'])
    union = upper + sum(parts[1:], F())
    diagnostic = dict(old_q1_bits=audit.bits(parts[0]), q1_bits=audit.bits(upper),
                      old_full_bits=audit.bits(sum(parts, F())), full_bits=audit.bits(union),
                      retained_other_occupancies_limit_bits=audit.bits(sum(parts[1:], F())),
                      remaining_shell_fraction=float(rest / upper))
    result = dict(status='IMT_K16_Q1_INJECTION_SHELL_REFINEMENT', instance=engine.identity(),
                  tilts=tilts, injection_weights=dict(sorted(Counter(engine.injection_weights).items())),
                  coefficients={str(w): model.base.encode(v) for w, v in best.items()},
                  q1_upper=model.base.encode(upper), union_upper=model.base.encode(union),
                  retained_other_upper=model.base.encode(sum(parts[1:], F())), diagnostic=diagnostic)
    if saved:
        for key in ('instance', 'q1_upper', 'union_upper', 'retained_other_upper', 'diagnostic'):
            assert result[key] == saved[key], key
        model.base.write_new(output.with_name(output.stem + '_replay.json'), dict(
            status='IMT_K16_Q1_INJECTION_SHELL_512_BIT_LINEAR_REPLAY_PASSED',
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
