"""Selected-point search with fresh outer authentication, not a certificate."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
from time import monotonic
from flint import arb, ctx
import dense_model


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path,
        default=Path('tmp/packed-hill/dense-r4-ekr-d10-complete-p256.json'))
    parser.add_argument('--feedback', choices=('weight5', 'bch16'), default='bch16')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--feedback-basis-seed', type=int)
    parser.add_argument('--refresh', choices=('uniform', 'transvections'), default='uniform')
    parser.add_argument('--updates', type=int, default=8)
    parser.add_argument('--means', nargs='+', default=['.032', '.104'])
    parser.add_argument('--scales', nargs='+', default=['1', '1/2', '1/4'])
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--variance-bins', type=int, default=16)
    parser.add_argument('--distance', type=Q, default=Q(1, 10))
    parser.add_argument('--outward', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if (args.output.exists() or any(not 0 < Q(p) < 1 for p in args.means)
            or any(Q(s) <= 0 for s in args.scales)):
        parser.error('fresh output, interior points and positive scales required')
    start = monotonic()
    model, scope, source = dense_model.fresh_model(args.source,
        feedback=args.feedback, seed=args.seed, refresh=args.refresh,
        updates=args.updates, precision=args.precision,
        variance_bins=args.variance_bins, distance=args.distance,
        feedback_basis_seed=args.feedback_basis_seed)
    record = dict(schema='packed-gl32-s16-dense-points-1', scope=scope,
        source=source, precision=args.precision, points=[],
        whole_code_certificate=False,
        proof_status='Selected comparison points only; no interval or sparse coverage.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def save():
        record['elapsed_seconds'] = monotonic()-start
        args.output.write_text(json.dumps(record, indent=2)+'\n')
    save()
    for mean in map(Q, args.means):
        if not model.root[0] <= mean <= model.root[1]:
            raise ValueError('comparison point outside the authenticated root domain')
        point = dict(mean=str(mean), trials=[])
        record['points'].append(point)
        def checkpoint(row, best):
            point['trials'].append(row)
            point['best_proposal'] = float(best[0])
            save()
            print('S16 TRIAL', str(mean), row, flush=True)
        ctx.prec = args.precision
        best = model.candidates((mean, mean), tuple(map(Q, args.scales)), checkpoint)
        point['witness'] = best[1]
        if args.outward:
            ctx.prec = args.precision
            upper = model.outward((mean, mean), best[1])
            if ctx.prec != args.precision or not upper.is_finite() or not upper > 0:
                raise ArithmeticError('positive finite outward bound at requested precision required')
            point['upper'] = list(map(int, upper.upper().man_exp()))
            point['log2_upper'] = str(upper.log()/arb(2).log())
            print('S16 CHECKED', str(mean), point['log2_upper'], flush=True)
        save()


if __name__ == '__main__':
    main()
