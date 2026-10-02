"""Screen lower distances for the measured shared-row, two-update encoder.

Every run regenerates the shared-row counting bounds and checks the positive
mixture exactly. Point probes are diagnostics, never a complete certificate.
A cover is a separate operation; its cell contributions still need to be
summed and combined with a disjoint sparse prefix before claiming a margin.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

from flint import arb, ctx
import birth_classes
import scalar_cover
import shared_mixture


def cutoff(distance):
    distance = Q(distance)
    if not 0 < distance < Q(1, 2):
        raise ValueError('relative distance in (0,1/2) required')
    return int(distance * scalar_cover.N)


def save(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + '\n')


def starting_mixture(record):
    """Extract rational coefficients only; no saved bound is a proof input."""
    if record.get('schema') == 'shared-gf16-relaxed-dense-1':
        rows = record.get('mixture')
    elif (record.get('schema') == 'shared-relaxed-alternative-mixtures-1'
            and isinstance(record.get('rows'), list) and len(record['rows']) == 1):
        rows = record['rows'][0].get('mixture')
    else:
        raise ValueError('one unambiguous shared-route comparison record required')
    if (not isinstance(rows, list) or not rows
            or any(not isinstance(row, dict) or set(row) != {'mass', 'activity'}
                   or any(type(row[k]) not in (str, int) for k in row) for row in rows)):
        raise ValueError('nonempty exact rational component list required')
    result = [dict(row) for row in rows]
    for row in result:
        if not Q(row['mass']) > 0 or not 0 < Q(row['activity']) <= 1:
            raise ValueError('positive mass and valid activity required')
    return result


def reuse_exact_census(previous, current):
    """Reuse only integer census data for the identical live inner object.

    Weighted operators, probability bounds, variance duals, and all proposal
    scratch stay local to each model. No census is loaded from a receipt.
    """
    if previous.data is not current.data:
        raise ValueError('census sharing requires identical actual inner data')
    for name in ('joint_return_cache', 'fiber_density_data'):
        if hasattr(previous, name):
            setattr(current, name, getattr(previous, name))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distances', nargs='+', default=['.05', '.07', '.08', '.09'])
    parser.add_argument('--means', nargs='+', default=['.001', '.008', '.032', '.064', '.128', '.25', '.5', '.75'])
    parser.add_argument('--minimum-groups', type=int, default=33)
    parser.add_argument('--cost-tilt', default='1/4')
    parser.add_argument('--base-tilt', default='3/16')
    parser.add_argument('--zero-bits', type=int, default=64)
    parser.add_argument('--prune', action='store_true', help='Exactly verify reductions to the positive comparison coefficients')
    parser.add_argument('--count-witnesses', nargs='+', type=Path, default=[],
        help='Regenerate and strengthen counts with exact dual witnesses; requires --prune')
    parser.add_argument('--count-refinement-iterations', type=int, default=0)
    parser.add_argument('--starting-mixture', type=Path,
        help='Componentwise prune a prior comparison; requires --count-witnesses and --prune')
    parser.add_argument('--variance-bins', type=int, default=16)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--target-bits', type=int, default=36)
    parser.add_argument('--max-cells', type=int, default=0, help='0 means point probes only; positive requests a separate cover')
    parser.add_argument('--max-depth', type=int, default=22)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    distances = list(map(Q, args.distances))
    means = list(map(Q, args.means))
    if (args.precision < 128 or args.target_bits < 1 or args.max_cells < 0
            or not 1 <= args.minimum_groups <= scalar_cover.G
            or not 1 <= args.variance_bins <= 64 or args.max_depth < 1
            or len(set(distances)) != len(distances)
            or any(not 0 <= mean <= 1 for mean in means)
            or not 0 <= args.count_refinement_iterations <= 30
            or (args.count_witnesses and not args.prune)
            or (args.count_refinement_iterations and not args.count_witnesses)
            or (args.starting_mixture and not args.count_witnesses)):
        parser.error('invalid precision, search limits, distances, or mean fractions')
    thresholds = [cutoff(d) for d in distances]
    if args.output.exists():
        parser.error('refusing to overwrite an existing run; choose a fresh output path')
    ctx.prec = args.precision
    count_metadata = {}
    starting_raw = None
    initial = None
    if args.starting_mixture:
        starting_raw = args.starting_mixture.read_bytes()
        initial = starting_mixture(json.loads(starting_raw))
    if args.count_witnesses:
        from shared_relaxed_alternative_counts import build_components
        components, caps, mixture, count_metadata = build_components(args.count_witnesses,
            iterations=args.count_refinement_iterations, zero_bits=args.zero_bits,
            cost_tilt=Q(args.cost_tilt), starting_mixture=initial)
        pruning = count_metadata['pruning']
    else:
        components, caps, mixture = shared_mixture.actual_components(
            coupled=True, zero_bits=args.zero_bits, cost_tilt=Q(args.cost_tilt))
        pruning = None
        if args.prune:
            from positive_prune import prune
            mixture, pruning = prune(caps, mixture, Q(args.cost_tilt))
            components = shared_mixture.as_components(mixture)
    data = birth_classes.actual(2)
    # Inner-census builders use their own precision. Restore the requested
    # context only after those builders finish.
    ctx.prec = args.precision
    fingerprint = hashlib.sha256(json.dumps([str(c) for c in caps]).encode()).hexdigest()
    record = dict(schema='shared-gf16-relaxed-dense-1', updates=2,
        K=1 << 20, N=scalar_cover.N, minimum_groups=args.minimum_groups,
        precision=args.precision, base_tilt=str(Q(args.base_tilt)),
        cost_tilt=str(Q(args.cost_tilt)), zero_bits=args.zero_bits,
        pruned=args.prune, pruning=pruning,
        variance_bins=args.variance_bins, target_bits=args.target_bits,
        max_cells=args.max_cells, max_depth=args.max_depth,
        cap_sha256=fingerprint,
        mixture=[dict(mass=str(c), activity=str(p)) for c, p in mixture],
        scope='shared four-row coordinate shuffle; independent GF16* packet labels; two IMT updates',
        proof_status='diagnostics or dense-only cover; no whole-code claim', results=[])
    if count_metadata:
        record.update(count_witnesses=count_metadata['count_witnesses'],
            count_refinement_iterations=count_metadata['count_refinement_iterations'])
    if starting_raw is not None:
        record['starting_mixture_source'] = dict(path=str(args.starting_mixture),
            sha256=hashlib.sha256(starting_raw).hexdigest(),
            scope='Search initialization only; saved output coefficients are independently checked against fresh caps.')
    print('EXACT SHARED MAJORANT', len(mixture), 'components', fingerprint, flush=True)
    previous_model = None
    for distance, threshold in zip(distances, thresholds):
        model = scalar_cover.Model(components, data, threshold, args.minimum_groups,
            Q(args.base_tilt), inner=birth_classes, variance_shuffle=True,
            variance_bins=args.variance_bins, regional_count=True)
        if previous_model is not None:
            reuse_exact_census(previous_model, model)
        previous_model = model
        model.proposal_stop_bits = args.target_bits + 2
        row = dict(distance=str(distance), threshold=threshold,
            root=list(map(str, model.root)), probes=[])
        record['results'].append(row)
        for mean in means:
            cell = (mean, mean)
            if model.empty(cell):
                row['probes'].append(dict(mean=str(mean), empty=True))
                continue
            score, witness = model.proposal(cell)
            ctx.prec = args.precision
            upper = model.outward(cell, witness)
            if not upper > 0:
                raise ArithmeticError('positive verified upper bound required')
            result = dict(mean=str(mean), proposal=score, witness=witness,
                upper=[int(v) for v in upper.upper().man_exp()],
                log2_upper=str(upper.log()/arb(2).log()))
            row['probes'].append(result)
            print('SHARED R2', 'distance', distance, 'mean', mean,
                'log2 upper', result['log2_upper'], flush=True)
            save(args.output, record)
        if args.max_cells:
            def checkpoint(state):
                row['cover'] = state
                save(args.output, record)
            row['cover'] = scalar_cover.run(model, args.max_cells,
                args.max_depth, args.target_bits, checkpoint=checkpoint)
            save(args.output, record)
    save(args.output, record)


if __name__ == '__main__':
    main()
