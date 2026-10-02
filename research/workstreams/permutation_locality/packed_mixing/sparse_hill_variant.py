"""Fresh sparse GL32 bounds for an explicitly selected inner/count variant.

No R2 endpoint or saved operator is imported. In particular, increasing the
number of random inner updates is not assumed to improve the low-weight
tail: additional linear mixing can also cause cancellation. Each requested
variant builds its own operators and authenticates its own count premises.

Successful rows cover their complete support domains. A complete q1..32
receipt still excludes dense occupancies and is not a whole-code certificate.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path

from flint import arb, ctx
import sparse_hill_q1 as q1
from outer_hill_incidence import authenticated_bch_cdf as incidence_counts
from outer_hill_intersection import authenticated_bch_cdf as intersection_counts


REFINEMENTS = {
    'incidence': 'packed-canonical-full32-h5-incidence-expected-cdf-1',
    'ekr': 'packed-canonical-full32-h5-incidence-ekr-expected-cdf-1',
}
DEFAULT_TILTS = ['.00016', '.00032', '.00064', '.001', '.0016', '.0032',
                 '.0064', '.008', '.012', '.016', '.024', '.032']


def validate(occupancies, distance, tilts, precision, target_bits, max_splits,
             output, updates, refinement):
    if (not occupancies or len(set(occupancies)) != len(occupancies)
            or any(type(q) is not int or not 1 <= q <= 32 for q in occupancies)
            or not 0 < Q(distance) < Q(1, 2)
            or not tilts or len(set(map(Q, tilts))) != len(tilts) or any(Q(t) <= 0 for t in tilts)
            or type(precision) is not int or precision < 128
            or type(target_bits) is not int or target_bits < 1
            or type(max_splits) is not int or max_splits < 0
            or type(updates) is not int or updates not in (2, 3, 4)
            or refinement not in REFINEMENTS or output.exists()):
        raise ValueError('distinct q=1..32, explicit supported variant, positive tilts, valid limits, and fresh output required')


def authenticated_counts(sparse, refinement):
    if refinement not in REFINEMENTS:
        raise ValueError('explicit supported count refinement required')
    counts, premises = (intersection_counts if refinement == 'ekr' else incidence_counts)()
    if (premises.get('schema') != REFINEMENTS[refinement]
            or premises.get('block_width') != 8 or premises.get('block_rows') != 4
            or premises.get('groups_per_outer_word') != 32):
        raise ValueError('fresh count premises do not match the requested GL32 variant')
    counts, floor = sparse.checked_counts(counts)
    return counts, floor, premises


def one_group(sparse, counts, operators, threshold, precision):
    """Fresh all-support placement moments, with no hardcoded support floor."""
    counts, floor = sparse.checked_counts(counts)
    if (not operators or type(precision) is not int or precision < 128
            or type(threshold) is not int or not 0 <= threshold < 1 << 21):
        raise ValueError('fresh operators, valid precision, and integer threshold required')
    ctx.prec = precision
    best, choices = [arb(1)]*257, [None]*257
    for (tilt, penalty), (exact, _) in operators.items():
        if penalty != '1' or Q(tilt) <= 0 or len(exact) < 2:
            raise ValueError('unweighted positive-tilt operators through occupancy1 required')
        ctx.prec = precision
        moments = sparse.single_group.support_moments(exact[0], exact[1])
        if len(moments) != 257:
            raise ValueError('every support0..256 requires a fresh moment')
        factor = (sparse.aq(Q(tilt))*threshold).exp()
        for u, moment in enumerate(moments):
            if not moment.is_finite() or not moment >= 0:
                raise ArithmeticError('finite nonnegative placement moments required')
            bound = sparse.up(factor*moment)
            if bound < best[u]:
                best[u], choices[u] = bound, tilt
    upper = sparse.up(2048*sparse.fold_cdf(counts, best))
    if ctx.prec != precision or not upper.is_finite() or not upper > 0:
        raise ArithmeticError('finite positive fresh q1 bound at requested precision required')
    return upper, dict(method='all-support exact placement', support_min=floor,
        support_choices=choices, support_probability_uppers=[sparse.endpoint(v) for v in best])


def _support_cover(*args, **kwargs):
    # Imported after the isolated sparse interface has established legacy
    # module paths. This does not change any legacy function or global.
    from support_cover import cover
    return cover(*args, **kwargs)


def _support_replay(*args, **kwargs):
    from support_cover import replay
    return replay(*args, **kwargs)


def validate_record(record):
    """Select witness-bearing search rows only after strict variant checks."""
    if (record.get('schema') != 'packed-gl32-sparse-hill-variant-1'
            or record.get('K') != 1 << 20 or record.get('N') != 1 << 21
            or record.get('group_count') != 2048 or record.get('outer') != 'BCH256128'
            or record.get('block_rows') != 4 or record.get('block_columns') != 8
            or record.get('mixing') != 'independent uniform GL32 per group/block; no additional GF16 stage'
            or record.get('routing') != 'independent shared column shuffle per group and independent regional shuffles'
            or record.get('count_refinement') not in REFINEMENTS
            or record.get('joint_return_through') != 2 or record.get('lazy_density_through') != 4):
        raise ValueError('matching sparse variant search scope required')
    inner = record.get('inner', {})
    updates = inner.get('updates')
    if (type(updates) is not int or updates not in (2, 3, 4)
            or inner != dict(t=128, s=19, updates=updates)
            or not 0 < Q(record['distance']) < Q(1, 2)
            or type(record.get('threshold')) is not int
            or record['threshold'] != (Q(record['distance'])*(1 << 21)).__floor__()
            or type(record.get('support_min')) is not int or not 1 <= record['support_min'] <= 256):
        raise ValueError('explicit supported inner variant and exact threshold required')
    requested, rows, tilts = record.get('requested'), record.get('results'), record.get('tilts')
    if (not isinstance(requested, list) or not requested or len(set(requested)) != len(requested)
            or any(type(q) is not int or not 1 <= q <= 32 for q in requested)
            or not isinstance(rows, list) or len({row['occupancy'] for row in rows}) != len(rows)
            or any(type(row['occupancy']) is not int or row['occupancy'] not in requested
                   or type(row.get('passed')) is not bool for row in rows)
            or not isinstance(tilts, list) or not tilts or len(set(map(Q, tilts))) != len(tilts)
            or any(not isinstance(t, str) or Q(t) <= 0 for t in tilts)):
        raise ValueError('distinct requested occupancies and positive exact witness tilts required')
    selected = [row for row in rows if row['passed']]
    for row in selected:
        details = row['details']
        if details.get('support_min') != record['support_min']:
            raise ValueError('support witness must retain the complete variant domain')
        if row['occupancy'] == 1:
            choices = details.get('support_choices')
            if (details.get('method') != 'all-support exact placement'
                    or not isinstance(choices, list) or len(choices) != 257
                    or any(t is not None and t not in tilts for t in choices)
                    or not any(t is not None for t in choices)):
                raise ValueError('complete q1 support witnesses required')
        elif (details.get('method') != 'full support-box CDF cover'
              or not details.get('leaves')
              or any(leaf['tilt'] not in tilts for leaf in details['leaves'])):
            raise ValueError('declared support-box witnesses required')
    return selected


def replay_records(sources, output, precision=384, target_bits=46,
                   distance=None, updates=None, refinement=None):
    """Rebuild variant counts/operators and replay witnesses; ignore saved bounds."""
    if (not sources or output.exists() or type(precision) is not int or precision < 128
            or type(target_bits) is not int or target_bits < 1):
        raise ValueError('source search receipts, fresh output, and valid replay limits required')
    records, entries, seen, tilts, metadata = [], [], set(), set(), []
    for source in sources:
        source_bytes = source.read_bytes()
        record = json.loads(source_bytes)
        selected = validate_record(record)
        metadata.append(dict(path=str(source.resolve()), sha256=hashlib.sha256(source_bytes).hexdigest()))
        for entry in selected:
            occupancy = entry['occupancy']
            if occupancy in seen:
                raise ValueError('duplicate successful occupancy across variant receipts')
            seen.add(occupancy)
            details = entry['details']
            if occupancy == 1:
                tilts.update(t for t in details['support_choices'] if t is not None)
            else:
                tilts.update(leaf['tilt'] for leaf in details['leaves'])
            entries.append(entry)
        records.append(record)
    if not entries:
        raise ValueError('at least one successful sparse witness required')
    updates = records[0]['inner']['updates'] if updates is None else updates
    refinement = records[0]['count_refinement'] if refinement is None else refinement
    if (type(updates) is not int or updates not in (2, 3, 4) or refinement not in REFINEMENTS
            or any(r['inner']['updates'] != updates or r['count_refinement'] != refinement for r in records)):
        raise ValueError('all sources must match the explicitly requested update/count variant')
    distance = min(Q(r['distance']) for r in records) if distance is None else Q(distance)
    if not 0 < distance < Q(1, 2) or any(distance > Q(r['distance']) for r in records):
        raise ValueError('only the same or a lower distance can be replayed')
    sparse = q1._sparse_interface()
    counts, floor, premises = authenticated_counts(sparse, refinement)
    digest = sparse.count_hash(counts)
    if any(r['count_sha256'] != digest or r['count_premises'] != premises
           or r['support_min'] != floor for r in records):
        raise ValueError('saved counts or support floor differ from fresh variant premises')
    threshold = (distance*(1 << 21)).__floor__()
    tilts = sorted(tilts, key=Q)
    args = sparse.old_args.build_args(max(seen), tilts, precision, 0, target_bits, None, updates)
    args.exact_feedback = True
    args.joint_return_through = 2
    args.lazy_density_through = 4
    operators = sparse.occupancy_birth_classes.build_operators(args)
    if not operators:
        raise ValueError('fresh nonempty replay operators required')
    terminal = sparse.np.ones(next(iter(operators.values()))[0][0].nrows())
    result = dict(schema='packed-gl32-sparse-hill-variant-replay-1', precision=precision,
        distance=str(distance), threshold=threshold, target_bits=target_bits,
        count_sha256=digest, count_premises=premises, count_refinement=refinement,
        scope={key: records[0][key] for key in ('K', 'N', 'outer', 'group_count',
            'block_rows', 'block_columns', 'mixing', 'routing', 'inner', 'support_min')},
        sources=metadata, tilts=tilts, joint_return_through=2, lazy_density_through=4,
        requested=sorted(seen), uncovered_sparse=sorted(set(range(1, 33))-seen),
        results=[], complete_requested=False, complete_sparse_prefix=False,
        aggregate_upper=None, aggregate_target_passed=False,
        note='Fresh variant-specific sparse replay. Dense occupancies and whole-code certification remain outside this receipt.')
    sparse.save(output, result)
    ctx.prec = precision
    total = arb(0)
    for entry in sorted(entries, key=lambda row: row['occupancy']):
        occupancy = entry['occupancy']
        ctx.prec = precision
        if occupancy == 1:
            upper, _ = one_group(sparse, counts, operators, threshold, precision)
        else:
            upper = _support_replay(entry['details'], operators, counts, terminal,
                groups=occupancy, cutoff=threshold, precision=precision, target_bits=target_bits)
        if not upper.is_finite() or not 0 < upper < arb(2)**(-target_bits) or ctx.prec != precision:
            raise ArithmeticError('fresh variant sparse replay failed')
        total = sparse.up(total+upper)
        result['results'].append(dict(occupancy=occupancy, passed=True, upper=sparse.endpoint(upper)))
        if len(result['results']) == len(entries):
            for source in metadata:
                if hashlib.sha256(Path(source['path']).read_bytes()).hexdigest() != source['sha256']:
                    raise ValueError('source search receipt changed during fresh replay')
        update_scope(result, total, sparse)
        sparse.save(output, result)
        print('FRESH SPARSE VARIANT R', updates, refinement, 'q', occupancy,
              'margin', -upper.log()/arb(2).log(), flush=True)
    return result


def update_scope(record, total, sparse):
    """Separate complete coverage from the exact aggregate target test."""
    record['complete_requested'] = (
        len(record['results']) == len(record['requested'])
        and all(row['passed'] for row in record['results']))
    record['uncovered_sparse'] = sorted(set(range(1, 33))-
        {row['occupancy'] for row in record['results'] if row['passed']})
    record['complete_sparse_prefix'] = record['complete_requested'] and not record['uncovered_sparse']
    record['successful_aggregate_upper'] = sparse.endpoint(total) if total > 0 else None
    record['aggregate_upper'] = sparse.endpoint(total) if record['complete_requested'] else None
    record['aggregate_target_passed'] = bool(record['complete_requested']
        and total > 0 and total < arb(2)**(-record['target_bits']))


def run(occupancies, distance, tilts, precision, target_bits, max_splits, output,
        updates=4, refinement='ekr'):
    validate(occupancies, distance, tilts, precision, target_bits, max_splits,
             output, updates, refinement)
    sparse = q1._sparse_interface()
    counts, floor, premises = authenticated_counts(sparse, refinement)
    threshold = (Q(distance)*(1 << 21)).__floor__()
    record = dict(schema='packed-gl32-sparse-hill-variant-1', K=1 << 20, N=1 << 21,
        group_count=2048, outer='BCH256128', block_rows=4, block_columns=8,
        mixing='independent uniform GL32 per group/block; no additional GF16 stage',
        routing='independent shared column shuffle per group and independent regional shuffles',
        inner=dict(t=128, s=19, updates=updates), distance=str(Q(distance)), threshold=threshold,
        precision=precision, target_bits=target_bits, tilts=list(tilts), max_splits=max_splits,
        joint_return_through=2, lazy_density_through=4,
        count_sha256=sparse.count_hash(counts), count_premises=premises, support_min=floor,
        requested=list(occupancies), results=[], process_id=os.getpid(),
        count_refinement=refinement, complete_requested=False, complete_sparse_prefix=False,
        uncovered_sparse=list(range(1, 33)), aggregate_upper=None,
        aggregate_target_passed=False,
        note='Fresh bounds only for successful listed occupancies. Dense occupancies and whole-code certification remain outside this receipt.')
    sparse.save(output, record)
    args = sparse.old_args.build_args(max(occupancies), tilts, precision, max_splits,
                                     target_bits, None, updates)
    args.exact_feedback = True
    args.analytic_gradient = True
    args.joint_return_through = 2
    args.lazy_density_through = 4
    # build_operators calls birth.actual(args.updates), including fresh local
    # return/lazy refinements, rather than adapting or importing R2 matrices.
    operators = sparse.occupancy_birth_classes.build_operators(args)
    if not operators:
        raise ValueError('fresh nonempty operator family required')
    terminal = sparse.np.ones(next(iter(operators.values()))[0][0].nrows())
    ctx.prec = precision
    total = arb(0)
    for occupancy in occupancies:
        args.groups = occupancy
        ctx.prec = precision
        if occupancy == 1:
            upper, details = one_group(sparse, counts, operators, threshold, precision)
        else:
            upper, details = _support_cover(args, operators, counts, terminal,
                                           cutoff=threshold, support_min=floor)
        passed = bool(upper is not None and upper.is_finite()
                      and upper > 0 and upper < arb(2)**(-target_bits))
        record['results'].append(dict(occupancy=occupancy, passed=passed,
            upper=sparse.endpoint(upper) if upper is not None else None, details=details))
        if passed:
            total = sparse.up(total+upper)
            print('SPARSE VARIANT R', updates, refinement, 'q', occupancy,
                  'margin', -upper.log()/arb(2).log(), flush=True)
        else:
            print('SPARSE VARIANT R', updates, refinement, 'q', occupancy,
                  'not closed within split budget', flush=True)
        update_scope(record, total, sparse)
        sparse.save(output, record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occupancies', type=int, nargs='+')
    parser.add_argument('--distance')
    parser.add_argument('--updates', type=int, choices=(2, 3, 4))
    parser.add_argument('--refinement', choices=tuple(REFINEMENTS))
    parser.add_argument('--precision', type=int)
    parser.add_argument('--target-bits', type=int, default=46)
    parser.add_argument('--max-splits', type=int)
    parser.add_argument('--tilts', nargs='+')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', type=Path, nargs='+', help='Freshly replay variant search witnesses; never saved bounds')
    args = parser.parse_args()
    print('SPARSE VARIANT PID', os.getpid(), flush=True)
    if args.replay:
        if args.occupancies is not None or args.max_splits is not None or args.tilts is not None:
            parser.error('replay rejects search occupancy, split, and tilt options')
        replay_records(args.replay, args.output, 384 if args.precision is None else args.precision, args.target_bits,
                       args.distance, args.updates, args.refinement)
    else:
        run(list(range(1, 33)) if args.occupancies is None else args.occupancies,
            '.10' if args.distance is None else args.distance,
            DEFAULT_TILTS if args.tilts is None else args.tilts,
            256 if args.precision is None else args.precision, args.target_bits,
            64 if args.max_splits is None else args.max_splits, args.output,
            4 if args.updates is None else args.updates,
            'ekr' if args.refinement is None else args.refinement)


if __name__ == '__main__':
    main()
