"""Bounded S16 point diagnostics using the retained outer comparison mixture.

No saved inner operator or numerical endpoint is reused. The outer mixture
is loaded from the named R4 search scope and is not freshly authenticated;
results are candidate-screening evidence, not new distance certificates.
The model remains t=128 so the existing 32-packet regional geometry applies.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

from flint import arb, ctx
import s16_maps

PARENT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PARENT))
sys.path.insert(0, str(PARENT/'gf16_packets'))
import birth_classes
import refresh_kernel
import scalar_cover as sc
import shared_mixture


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('tmp/packed-hill/dense-r4-ekr-d10-complete-p256.json'))
    parser.add_argument('--feedback', choices=('weight5', 'bch16'), default='weight5')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--updates', type=int, nargs='+', default=[4, 8])
    parser.add_argument('--means', nargs='+', default=['.032', '.104'])
    parser.add_argument('--kernel', choices=('refresh', 'birth'), default='birth')
    parser.add_argument('--regional', action='store_true')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--outward', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or any(not 1 <= r <= 32 for r in args.updates):
        parser.error('fresh output and 1..32 updates required')
    raw = args.source.read_bytes()
    source = json.loads(raw)
    scope = source['scope']
    if (scope['K'] != 1 << 20 or scope['N'] != 1 << 21 or scope['block_width'] != 8
            or scope['minimum_groups'] != 33 or scope['threshold'] != 209715):
        raise ValueError('expected fixed K20 width-eight 10% dense comparison scope')
    mixture = [(Q(row['mass']), Q(row['activity'])) for row in scope['mixture']]
    images, columns, maps = s16_maps.candidate(args.feedback, args.seed)
    inner = birth_classes if args.kernel == 'birth' else refresh_kernel
    print('S16 preparing exact profiles', args.feedback, args.kernel, flush=True)
    data = inner.prepare(images, columns, 16, args.updates[0])
    print('S16 profiles ready', len(data['records']), 'characters;',len(data['histograms']), 'images', flush=True)
    record = dict(schema='s16-dense-candidate-points-1', maps=maps,
        source=dict(path=str(args.source.resolve()), sha256=hashlib.sha256(raw).hexdigest()),
        outer_scope=scope, precision=args.precision, kernel=args.kernel, points=[],
        proof_status='Selected comparison points only, with retained outer mixture. No outer premise authentication, complete domain, or sparse coverage.',
        whole_code_certificate=False, regional=args.regional)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def save():
        args.output.write_text(json.dumps(record, indent=2)+'\n')
    save()
    for updates in args.updates:
        selected = dict(data, updates=updates)
        ctx.prec = args.precision
        model = sc.Model(shared_mixture.as_components(mixture), selected, 209715, 33, Q(3, 16),
            inner=inner, variance_shuffle=True, variance_bins=16,
            regional_count=args.regional and args.kernel == 'birth')
        for mean in map(Q, args.means):
            cell = (mean, mean)
            best = model.propose_with(cell, model.tilt)
            trials = [dict(method='iid', log2_proposal=float(best[0]))]
            print('S16 POINT', args.feedback, 'r',updates,'mean',str(mean),'iid',best[0],flush=True)
            import variance_partition
            variance = variance_partition.propose(model, cell, best)
            trials.append(dict(method='variance', log2_proposal=float(variance[0])))
            if variance[0] < best[0]:
                best = variance
            print('S16 POINT', args.feedback, 'r',updates,'mean',str(mean),'variance',variance[0],flush=True)
            if args.regional:
                if args.kernel != 'birth':
                    raise ValueError('regional refinement requires the birth-class kernel')
                import regional_count
                witness = copy.deepcopy(variance[1])
                witness.update(regional_tilted_atom=True, regional_fine_tilts=True,
                    regional_tilted_variance=True, regional_direct_counts=True, regional_exact_zero=True,
                    regional_feedback_classes_from=3, regional_feedback_classes_through=32,
                    regional_feedback_uniform_classes=True, regional_feedback_uniform_replace=True)
                # No *_actual map-dependent adapter: the legacy lazy-density
                # and joint-return entry points authenticate the 19-state maps.
                regional = regional_count.propose(model, cell, witness)
                trials.append(dict(method='regional_without_legacy_map_refinements', log2_proposal=float(regional[0])))
                if regional[0] < best[0]:
                    best = regional
                print('S16 POINT', args.feedback, 'r',updates,'mean',str(mean),'regional',regional[0],flush=True)
            point = dict(updates=updates, mean=str(mean), trials=trials,
                         log2_proposal=float(best[0]), witness=best[1])
            if args.outward:
                ctx.prec = args.precision
                upper = model.outward(cell, best[1])
                if not upper.is_finite() or not upper > 0:
                    raise ArithmeticError('positive finite outward point bound required')
                point.update(upper=list(map(int, upper.upper().man_exp())),
                             log2_upper=str(upper.log()/arb(2).log()))
            record['points'].append(point)
            save()


if __name__ == '__main__':
    main()
