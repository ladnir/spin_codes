"""BCH-256 Q1 certificates only; 512-bit replay uses linear epoch iteration."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from flint import arb, ctx
import bch_model as model


def run(output, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'OUTWARD_BCH256_ASYMMETRIC_Q1_ONLY'
    else:
        assert not output.exists(), 'Use a fresh output path'
    ctx.prec = 512 if verify else 256
    results = []
    for index, exponent in enumerate((16, 18, 20)):
        engine = model.Engine(exponent)
        old = saved['results'][index] if saved else None
        if old:
            assert old['instance'] == engine.identity()
            tilts = old['log_tilts']
        else:
            tilts = [str(F.from_float(math.log(a/(100*engine.length)))) for a in (235,264,295,332,376,425)]
        best = {w: arb(engine.length) for w in model.base.WEIGHTS}
        for z in tilts:
            lam = model.number(F(z)).exp()
            co = engine.q1(F(z), linear=verify)
            factor = engine.length*(engine.cutoff*lam).exp()
            for w in best:
                best[w] = min(best[w], model.up(co[w]*factor))
        if old:
            bounds = old['coefficient_upper']
            assert set(map(int, bounds)) == set(best)
            assert all(v <= model.unpack(bounds[str(w)]) for w, v in best.items())
        else:
            bounds = {str(w): model.pack(v) for w, v in best.items()}
        exact = {int(w): model.exact(model.unpack(v)) for w, v in bounds.items()}
        upper, factor, rest = model.base.bch_bound(exact)
        assert upper < F(1, 1 << 40)
        result = dict(instance=engine.identity(), log_tilts=tilts, coefficient_upper=bounds,
                      upper=model.base.encode(upper), factor=model.base.encode(factor), rest=model.base.encode(rest),
                      margin_bits=math.log2(upper.denominator)-math.log2(upper.numerator),
                      covered_occupancies=[1, 1], full_distance_proved=False)
        if old:
            assert result == old
        results.append(result)
        print(ctx.prec, 'K=2^'+str(exponent), 'Q1 margin', result['margin_bits'], flush=True)
    if saved:
        model.base.write_new(output.with_name(output.stem+'_replay.json'),
                             dict(status='BCH256_ASYMMETRIC_Q1_512_BIT_LINEAR_REPLAY_PASSED',
                                  producer_sha256=model.base.sha(output), full_distance_proved=False))
    else:
        model.base.write_new(output, dict(status='OUTWARD_BCH256_ASYMMETRIC_Q1_ONLY', precision_bits=256,
                                         results=results, source_sha256=model.sources()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=model.HERE/'Q1_CERTIFICATE.json')
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    run(args.output.resolve(), args.verify)
