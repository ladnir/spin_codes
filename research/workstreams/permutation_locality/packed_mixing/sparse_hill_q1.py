"""Fresh occupancy-one tilt diagnostics, with canonical-block contributions.

Only q=1 is covered. The canonical H caps count fixed outer inputs; local
GL32 randomness gives the exact support law conditional on H. Independent
column/region permutations then permit the existing fixed-support transfer.
No numerical endpoint from an older receipt is used.
"""
import argparse
from fractions import Fraction as Q
import importlib.util
import json
import os
from pathlib import Path
import sys

from flint import arb, ctx
from local_models import checked_pmf, convolve, full_block

DEFAULT_TILTS = ['.00032', '.001', '.0032', '.00016', '.00020', '.00024',
    '.00028', '.00036', '.00040', '.00048', '.00064', '.00080', '.00128']


def aq(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def up(value):
    mantissa, exponent = value.upper().man_exp()
    return arb(int(mantissa))*arb(2)**int(exponent)


def endpoint(value):
    return [int(v) for v in up(value).man_exp()]


def support_laws(local, maximum_h):
    local = checked_pmf(local)
    if local[0] != 0 or type(maximum_h) is not int or maximum_h < 1:
        raise ValueError('positive local support and positive H range required')
    laws = [(Q(1),)]
    for _ in range(maximum_h):
        laws.append(convolve(laws[-1], local))
    return laws


def canonical_terms(caps, laws, weights, groups=2048):
    """Directed Abel fold over H, retaining the terms for diagnosis.

    For each H first average the valid support-specific event bounds over
    the exact GL32 law. The suffix maximum makes these averages decreasing;
    only then may differences of cumulative caps be used in the Abel sum.
    """
    caps = tuple(map(Q, caps))
    if (len(caps) != len(laws) or len(caps) < 2 or caps[0] != 0
            or any(a < 0 or a > b for a, b in zip(caps, caps[1:]))
            or type(groups) is not int or groups < 1
            or len(weights) != len(laws[-1])
            or any(not w.is_finite() or w < 0 for w in weights)):
        raise ValueError('matching cumulative caps, complete laws, and nonnegative finite weights required')
    averages = []
    for law in laws:
        checked_pmf(law)
        averages.append(up(sum((aq(p)*weights[u] for u, p in enumerate(law) if p), arb(0))))
    tails = [arb(0)]*len(averages)
    for h in range(len(averages)-1, -1, -1):
        tails[h] = max(averages[h], tails[h+1] if h+1 < len(tails) else arb(0))
    terms = [up(groups*aq(caps[h]-caps[h-1])*tails[h]) for h in range(1, len(caps))]
    return up(sum(terms, arb(0))), terms


def _sparse_interface():
    # Do not occupy legacy global names such as sparse or dense_cover.
    path = Path(__file__).resolve().with_name('sparse.py')
    spec = importlib.util.spec_from_file_location('packed_sparse_hill_interface', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def authenticated_counts(sparse, h5_exact=False):
    """Intersect only canonical rank CDFs with the freshly checked H5 count."""
    if type(h5_exact) is not bool:
        raise ValueError('h5_exact must be bool')
    counts, premises = sparse.authenticated_bch_cdf(104, True, 8)
    if h5_exact:
        from outer_hill_octets import check_production
        from canonical_counts import transport_cdf, total_caps
        proof = check_production()
        if (proof['dimension'], proof['length'], proof['block_width'], proof['minimum_distance']) != (128, 256, 8, 38):
            raise ArithmeticError('production H5 proof has different code premises')
        ranks = [row[:] for row in premises['canonical_rank_cdfs']]
        for row, bound in zip(ranks, proof['rank_cumulative_caps_at5']):
            row[5] = min(row[5], bound)
            for h in range(4, -1, -1):
                row[h] = min(row[h], row[h+1])
        premises = dict(premises, canonical_rank_cdfs=ranks,
            canonical_cdf=list(total_caps(ranks)), canonical_h5_exact=proof)
        counts = transport_cdf(premises['canonical_cdf'], full_block(8))
    return sparse.checked_counts(counts)[0], premises


def run(distances, tilts, precision, output, h5_exact=False):
    if (output.exists() or type(precision) is not int or precision < 128
            or not distances or any(not 0 < Q(d) < Q(1, 2) for d in distances)
            or not tilts or len(set(map(Q, tilts))) != len(tilts) or any(Q(t) <= 0 for t in tilts)):
        raise ValueError('new output, exact positive tilts/distances, and precision >=128 required')
    sparse = _sparse_interface()
    counts, premises = authenticated_counts(sparse, h5_exact)
    laws = support_laws(full_block(8), 32)
    record = dict(schema='packed-gl32-q1-hill-1', precision=precision, updates=2,
        count_sha256=sparse.count_hash(counts), count_premises=premises,
        h5_exact=h5_exact, minimum_possible_support=next(u for u, v in enumerate(counts) if v),
        joint_return_through=1, lazy_density_through=1, tilts=list(tilts),
        completed_tilts=[], process_id=os.getpid(), results=[],
        proof_status='Fresh occupancy-one bounds only. No sparse-prefix or whole-code claim.')
    best = {str(Q(d)): [arb(1)]*257 for d in distances}
    choices = {d: [None]*257 for d in best}

    def save():
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(record, indent=2)+'\n')

    save()
    import occupancy_birth_classes as birth
    import occupancy_model
    from occupancy_memory import rounded
    import return_moment
    import lazy_density
    data = birth.actual(2)
    ctx.prec = precision
    census = return_moment.actual_census(data, 1)
    for tilt in tilts:
        ctx.prec = precision
        local = birth.outward(data, tilt)
        z = (-aq(Q(tilt))).exp()
        local = return_moment.refine_class_returns(data, local, census, z)
        local = lazy_density.refine_actual(data, local, z, 1)
        exact = occupancy_model.placement(local, rounding=rounded, maximum_groups=1)
        moments = sparse.single_group.support_moments(exact[0], exact[1])
        record['results'] = []
        for distance, weights in best.items():
            threshold = int(Q(distance)*(1 << 21))
            factor = (aq(Q(tilt))*threshold).exp()
            for u, moment in enumerate(moments):
                bound = up(factor*moment)
                if bound < weights[u]:
                    weights[u], choices[distance][u] = bound, tilt
            original = up(2048*sparse.fold_cdf(counts, weights))
            direct, terms = canonical_terms(premises['canonical_cdf'], laws, weights)
            if ctx.prec != precision or not 0 < direct or not direct.is_finite():
                raise ArithmeticError('finite positive fresh q1 bound at requested precision required')
            record['results'].append(dict(distance=distance, threshold=threshold,
                cdf_upper=endpoint(original), direct_upper=endpoint(direct),
                margin=str(-direct.log()/arb(2).log()), canonical_H_terms=[endpoint(v) for v in terms],
                support_choices=choices[distance][:], support_probability_uppers=[endpoint(v) for v in weights]))
            print('Q1 HILL tilt', tilt, 'distance', distance, 'margin', record['results'][-1]['margin'], flush=True)
        record['completed_tilts'].append(tilt)
        save()
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distances', nargs='+', default=['.095', '.096', '.10'])
    parser.add_argument('--tilts', nargs='+', default=DEFAULT_TILTS)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--h5-exact', action='store_true', help='Fresh exact canonical H<=5 enumeration before transport')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print('Q1 HILL PID', os.getpid(), flush=True)
    run(args.distances, args.tilts, args.precision, args.output, args.h5_exact)


if __name__ == '__main__':
    main()
