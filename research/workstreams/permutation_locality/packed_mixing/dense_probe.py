"""Outward point diagnostics for canonical full-block mixing, not a certificate.

The forward outer is four BCH words followed by independent uniform GL32
maps on canonical four-row/eight-column blocks. An independent shared column
shuffle follows. No separate GF16 packet randomizer is required. All groups
have independent mixer setup; averaging one local enumerator is not valid
for a mixer table shared across all groups.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import sys

from canonical_counts import canonical_rank_caps, total_caps, transport_cdf
from local_models import full_block


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance', default='.1')
    parser.add_argument('--means', nargs='*', default=['.032'])
    parser.add_argument('--updates', type=int, choices=(2, 3, 4), default=2)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--last-lp', type=int, default=104)
    parser.add_argument('--refined', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('new output path required; earlier diagnostics are preserved')
    if not 0 < Q(args.distance) < Q(1, 2) or args.precision < 128:
        parser.error('distance in (0,1/2) and precision >=128 required')
    if not 0 <= args.last_lp <= 256 or any(not 0 <= Q(x) <= 1 for x in args.means):
        parser.error('last-lp in0..256 and means in[0,1] required')

    parent = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(parent))
    sys.path.insert(0, str(parent / 'gf16_packets'))
    from bch_joint_support import authenticated_caps
    from shortened_bound import dimension_caps
    from flint import arb, ctx
    import shared_mixture
    from positive_prune import prune

    spectrum = authenticated_caps()
    if next(w for w in range(1, 257) if spectrum[w]) != 38:
        raise ArithmeticError('unexpected BCH minimum-distance premise')
    dimensions = dimension_caps(last_lp=args.last_lp)
    if args.refined:
        from shortening_polynomial import improve_dimensions
        from dual_shortening import improve_dimensions as dual_dimensions
        dimensions = dual_dimensions(improve_dimensions(dimensions))
    block_caps = total_caps(canonical_rank_caps(dimensions))
    caps = transport_cdf(block_caps, full_block(8))
    centers = sorted({Q(u, 256) for u in range(1, 257, 4)} | {Q(1)})
    # The expected CDF itself bounds each expected shell. Its differences do
    # NOT: using them as shell bounds would invalidate this comparison.
    mixture = shared_mixture.envelope(caps, centers, zero_bits=64, cost_tilt=Q(1, 4))
    mixture, pruning = prune(caps, mixture, Q(1, 4))
    shared_mixture.verify(caps, mixture)
    components = shared_mixture.as_components(mixture)
    print('CANONICAL GL32: exact expected-shell majorant checked;',
          len(mixture), 'components', flush=True)
    record = dict(schema='canonical-packed-gl32-dense-points-1',
        proof_status='Expected outer caps and selected inner point diagnostics only; not a distance certificate.',
        distance=args.distance, updates=args.updates, precision=args.precision,
        last_lp=args.last_lp, refined=args.refined, dimensions=dimensions,
        minimum_groups=33, maximum_groups=2048, base_tilt='3/16',
        variance_shuffle=True, variance_bins=16, regional_count=True,
        uncovered_scope='Occupancies 1..32 and all untested mean cells remain uncovered.',
        block_cdf_caps=list(map(str, block_caps)),
        expected_output_cdf_caps=list(map(str, caps)),
        mixture=[dict(mass=str(c), activity=str(p)) for c, p in mixture],
        pruning=pruning, probes=[])

    def save():
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(record, indent=2) + '\n')

    save()
    if args.means:
        import birth_classes
        import scalar_cover as sc
        data = birth_classes.actual(args.updates)
        ctx.prec = args.precision
        model = sc.Model(components, data, int(Q(args.distance) * sc.N), 33, Q(3, 16),
            inner=birth_classes, variance_shuffle=True, variance_bins=16, regional_count=True)
        model.proposal_stop_bits = 42
        for mean in map(Q, args.means):
            cell = (mean, mean)
            if model.empty(cell):
                record['probes'].append(dict(mean=str(mean), empty=True))
                print('CANONICAL GL32 POINT', mean, 'empty comparison cell', flush=True)
                save()
                continue
            score, witness = model.proposal(cell)
            ctx.prec = args.precision
            upper = model.outward(cell, witness)
            row = dict(mean=str(mean), proposal=score, witness=witness,
                upper=[int(x) for x in upper.upper().man_exp()],
                log2_upper=str(upper.log() / arb(2).log()))
            record['probes'].append(row)
            print('CANONICAL GL32 POINT', args.distance, mean, row['log2_upper'], flush=True)
            save()


if __name__ == '__main__':
    main()
