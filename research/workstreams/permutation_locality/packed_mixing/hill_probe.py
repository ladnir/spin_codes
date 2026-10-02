"""Fresh pointwise dense hill climb; never a full distance certificate.

Use a saved comparison mixture only after authenticating all its premises.
Rebuild the output threshold, variance partition, and local operators. Search
a broader output-tilt range than the closure driver's short fallback list.
Only the best candidate at each requested point is evaluated outward.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import sys

from flint import arb, ctx
import cell_search

SCHEMA = 'packed-canonical-gl32-dense-hill-points-1'


def retarget(source, distance):
    """Keep count premises, but never import old points or accepted bounds."""
    cell_search.dense.validate_record(source)
    if source['comparison'] != 'direct-expected-shell-majorant':
        raise ValueError('direct expected-shell comparison required')
    result = {key: copy.deepcopy(source[key]) for key in cell_search.SCOPE_FIELDS}
    result.update(schema=cell_search.dense.SCHEMA, distance=str(Q(distance)),
        threshold=cell_search.dense.cutoff(distance),
        cover=dict(leaves={}, unresolved={'': dict(cell=copy.deepcopy(source['root']))}))
    cell_search.dense.validate_record(result)
    return result


def refined_model(source, distance, precision, variance_bins, *, updates=2, intersection=False):
    """Rebuild the entire comparison from fresh H5/incidence count premises.

    The source supplies only validated construction/scope fields. Its old
    mixture, endpoints, partitions, and point witnesses are not reused.
    This context deliberately has a new schema, rejected by the legacy
    dense replay that authenticates the older count interface.
    """
    scope = retarget(source, distance)
    if type(updates) is not int or updates not in (2, 3, 4) or type(intersection) is not bool:
        raise ValueError('two through four actual inner updates and boolean intersection option required')
    if intersection:
        from outer_hill_intersection import authenticated_bch_cdf
    else:
        from outer_hill_incidence import authenticated_bch_cdf
    from monotone import transport_shells
    from local_models import full_block
    caps, premises = authenticated_bch_cdf()
    comparison = transport_shells(premises['canonical_cdf'], full_block(8))
    parent = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(parent))
    sys.path.insert(0, str(parent/'gf16_packets'))
    mixture, info = cell_search.dense.exact_mixture(comparison)
    import birth_classes
    import shared_mixture
    import scalar_cover as sc
    data = birth_classes.actual(updates)
    ctx.prec = precision
    model = sc.Model(shared_mixture.as_components(mixture), data, scope['threshold'],
        scope['minimum_groups'], Q(3, 16), inner=birth_classes,
        variance_shuffle=True, variance_bins=variance_bins, regional_count=True)
    scope.update(schema='packed-canonical-gl32-hill-dense-context-1',
        updates=updates, ensemble=f'canonical-gl32-width8-shared4-r{updates}',
        root=list(map(str, model.root)), expected_cdf_sha256=cell_search.dense.fingerprint(caps),
        comparison_caps_sha256=cell_search.dense.fingerprint(comparison),
        outer_premises=premises, variance_bins=variance_bins,
        mixture=[dict(mass=str(c), activity=str(p)) for c, p in mixture],
        mixture_verification=info,
        cover=dict(leaves={}, unresolved={'': dict(cell=list(map(str, model.root)))}))
    return model, scope


def candidates(model, cell, scales, variants, callback=None):
    """Force regional refinement even when an earlier candidate passes."""
    import variance_partition as variance
    import regional_count as regional
    if not scales or any(Q(x) <= 0 for x in scales):
        raise ValueError('positive output-tilt scales required')
    if not variants or any(x not in ('classified', 'uniform') for x in variants):
        raise ValueError('classified or uniform regional candidates required')
    base = model.propose_with(cell, model.tilt)
    variance_candidate = variance.propose(model, cell, base)
    # Keep the cheap alternatives as well: stronger local envelopes need
    # not compensate for every comparison loss at every output tilt.
    best = min((base, variance_candidate), key=lambda candidate: candidate[0])
    best = (best[0], copy.deepcopy(best[1]))
    witness = copy.deepcopy(variance_candidate[1])
    witness.update(regional_tilted_atom=True, regional_fine_tilts=True,
        regional_tilted_variance=True, regional_direct_counts=True,
        regional_exact_zero=True, regional_lazy_density_through=6,
        regional_joint_return_through=3, regional_feedback_classes_from=3,
        regional_feedback_classes_through=32)
    _, witness = regional.prepare_witness(model, cell, witness)
    lam = Q(witness['parameters'][0])
    trials = []
    for scale in map(Q, scales):
        for variant in variants:
            trial = copy.deepcopy(witness)
            trial['parameters'][0] = str(lam*scale)
            if variant == 'uniform':
                trial.update(regional_feedback_uniform_classes=True,
                    regional_feedback_uniform_replace=True)
            score, checked = regional.propose(model, cell, trial)
            if not math.isfinite(score):
                raise ArithmeticError('finite screening score required')
            row = dict(scale=str(scale), variant=variant, output_tilt=str(lam*scale),
                proposal=float(score))
            trials.append(row)
            if best is None or score < best[0]:
                best = (score, checked)
            if callback:
                callback(row, best)
    return trials, best


def check_point(model, cell, best, precision):
    """A fresh outward bound at a point is not coverage of its neighbours."""
    ctx.prec = precision
    upper = model.outward(cell, best[1])
    if ctx.prec != precision or not upper.is_finite() or not upper > 0:
        raise ArithmeticError('finite positive fresh outward bound required')
    return dict(witness=best[1], proposal=float(best[0]),
        upper=[int(v) for v in upper.upper().man_exp()],
        log2_upper=str(upper.log()/arb(2).log()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--distance', default='.1')
    parser.add_argument('--means', nargs='+', default=['.032', '.104'])
    parser.add_argument('--tilt-scales', nargs='+', default=['1/8', '1/4', '1/2', '3/4', '1'])
    parser.add_argument('--variants', nargs='+', choices=('classified', 'uniform'),
        default=['classified', 'uniform'])
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--variance-bins', type=int, default=16)
    parser.add_argument('--outer-refinement', action='store_true',
        help='Fresh H5 exclusion and incidence counts, followed by a new positive mixture')
    parser.add_argument('--intersection', action='store_true',
        help='Also use the separately checked EKR count refinement')
    parser.add_argument('--updates', type=int, choices=(2, 3, 4), default=2,
        help='Actual inner updates; changing this changes the construction')
    parser.add_argument('--quick', action='store_true',
        help='Use the ordinary adaptive proposer, stopping after a sufficiently strong bound')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    temporary = args.output.with_name(args.output.name+'.tmp')
    if (args.output.exists() or temporary.exists() or args.precision < 128
            or not 1 <= args.variance_bins <= 64
            or any(not 0 < Q(x) < 1 for x in args.means)
            or any(Q(x) <= 0 for x in args.tilt_scales)):
        parser.error('fresh output, interior means, positive tilts and valid precision/bins required')
    if (args.intersection or args.updates != 2) and not args.outer_refinement:
        parser.error('intersection or changed inner requires the freshly refined construction factory')
    raw = args.source.read_bytes()
    if args.outer_refinement:
        model, scope = refined_model(json.loads(raw), args.distance, args.precision, args.variance_bins,
            updates=args.updates, intersection=args.intersection)
    else:
        scope = retarget(json.loads(raw), args.distance)
        model = cell_search.fresh_model(scope, args.precision, 52)
    # Rebuilding the partition below is essential: changing this field does
    # not refine a variance partition already stored inside a witness.
    model.variance_bins = args.variance_bins
    model.proposal_stop_bits = 54
    record = dict(schema=SCHEMA, scope=scope,
        source=dict(path=str(args.source.resolve()), sha256=hashlib.sha256(raw).hexdigest()),
        precision=args.precision, variance_bins=args.variance_bins,
        count_refinement=('H5 exact + incidence + EKR' if args.intersection else
            'H5 exact + incidence' if args.outer_refinement else 'original closure counts'),
        changed_construction=args.updates != 2, search='adaptive' if args.quick else 'forced regional grid',
        partial_only=True, whole_code_certificate=False,
        proof_status='Selected points only. Floating trials are not proof bounds.',
        points=[])

    def save():
        if args.source.read_bytes() != raw:
            raise ValueError('source changed during search')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_text(json.dumps(record, indent=2)+'\n')
        temporary.replace(args.output)

    save()
    for mean in map(Q, args.means):
        if not model.root[0] <= mean <= model.root[1]:
            raise ValueError('requested point outside global comparison domain')
        ctx.prec = args.precision
        point = dict(mean=str(mean), cell=[str(mean), str(mean)], trials=[])
        record['points'].append(point)

        def checkpoint(row, best):
            point['trials'].append(row)
            point['best_proposal'] = float(best[0])
            save()
            print('HILL SCREEN', str(mean), row, flush=True)

        if args.quick:
            best = model.proposal((mean, mean))
        else:
            _, best = candidates(model, (mean, mean), args.tilt_scales, args.variants, checkpoint)
        point['checked'] = check_point(model, (mean, mean), best, args.precision)
        save()
        print('HILL CHECKED', str(mean), point['checked']['log2_upper'], flush=True)


if __name__ == '__main__':
    main()
