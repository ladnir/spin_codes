"""Fresh whole-code replay for the shared-four-row, GF16, two-update route.

Saved numerical bounds and acceptance labels are not proof inputs. Rebuild
the exact shared-support comparison, replay every dense cell, and regenerate
the complementary sparse occupancies before summing exact dyadic endpoints.
The minimum margin is a caller choice; it does not change these checks.
"""
import argparse
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
from types import SimpleNamespace

from flint import arb, ctx
import birth_classes
import scalar_cover as sc
import shared_mixture
import shared_sparse
from audit_complete import dyadic


SCHEMA = 'shared-gf16-relaxed-complete-replay-1'
DENSE_SCHEMA = 'shared-gf16-relaxed-dense-1'
ENSEMBLE = 'shared4-gf16-r2'
TILTS = ('.00032', '.001', '.0032', '.0064', '.008', '.012', '.016',
         '.024', '.032', '.048', '.064', '.096')


def require_uncached_verification():
    import regional_count
    if getattr(regional_count.outward, '_shared_search_region_cache', False):
        raise ArithmeticError('full verification must run outside the optional search-cache context')


def _exact_endpoint(value):
    value = Q(value)
    denominator = value.denominator
    if value < 0 or denominator & (denominator-1):
        raise ValueError('nonnegative exact dyadic endpoint required')
    return [value.numerator, 1-denominator.bit_length()]


def compact_endpoint(value, bits=None):
    """Round a nonnegative exact dyadic sum upward to a compact endpoint.

    The sum is formed exactly before this final rounding. This avoids giant
    JSON integers when accepted terms have widely separated exponents.
    """
    bits=ctx.prec+32 if bits is None else bits
    if type(bits) is not int or bits<32:
        raise ValueError('integer endpoint precision at least 32 required')
    n,e=_exact_endpoint(value)
    if n.bit_length()>bits:
        shift=n.bit_length()-bits
        n=(n+(1 << shift)-1) >> shift
        e+=shift
    if n:
        shift=(n&-n).bit_length()-1
        n>>=shift;e+=shift
    return [n,e]


def endpoint(value):
    """Compact outward endpoint; exact addition precedes final rounding."""
    return compact_endpoint(value)


def validate_record(record, index=0, complete=False):
    if (not isinstance(record, dict) or record.get('schema') != DENSE_SCHEMA
            or type(record.get('updates')) is not int or record['updates'] != 2
            or record.get('K') != 1 << 20 or record.get('N') != 1 << 21
            or type(record.get('minimum_groups')) is not int
            or not 1 <= record['minimum_groups'] <= 2048
            or type(record.get('zero_bits')) is not int or not 0 <= record['zero_bits'] <= 512
            or type(record.get('variance_bins')) is not int or not 1 <= record['variance_bins'] <= 64
            or type(record.get('pruned', False)) is not bool
            or not isinstance(record.get('results'), list)
            or type(index) is not int or not 0 <= index < len(record['results'])):
        raise ValueError('matching shared-row two-update construction metadata required')
    if Q(record['base_tilt']) <= 0 or not 0 < Q(record['cost_tilt']) <= 1:
        raise ValueError('valid positive comparison tilts required')
    row = record['results'][index]
    if (not isinstance(row, dict) or type(row.get('threshold')) is not int
            or not 0 < Q(row['distance']) < Q(1, 2)
            or row['threshold'] != int(Q(row['distance'])*(1 << 21))
            or not isinstance(row.get('root'), list) or len(row['root']) != 2
            or not 0 <= Q(row['root'][0]) <= Q(row['root'][1]) <= 1):
        raise ValueError('exact consistent cutoff, distance, and root required')
    if complete:
        cover = row.get('cover', {})
        leaves = cover.get('leaves')
        if (record.get('screen_only', False) or row.get('screen_only', False)
                or cover.get('unresolved') != {} or not isinstance(leaves, dict) or not leaves
                or any(not isinstance(v, dict) or not isinstance(v.get('witness'), dict)
                       for v in leaves.values())):
            raise ValueError('complete dense partition with actual witnesses required')
        sc.partition(SimpleNamespace(root=(Q(0), Q(1))), leaves, {})
    return row


def checked_mixture(record, caps, original):
    """Authenticate a stored pruned candidate with exact shell inequalities.

    The LP and its saved success labels are not inputs to verification.
    Restricting the candidate to reductions of existing components also
    preserves the original empty-input constraints and component geometry.
    """
    if type(record.get('pruned', False)) is not bool:
        raise ValueError('explicit boolean pruning selector required')
    encoded = [dict(mass=str(c), activity=str(p)) for c, p in original]
    if not record.get('pruned', False):
        if record.get('mixture') != encoded:
            raise ValueError('saved comparison differs from the original authenticated envelope')
        return original
    saved = record.get('mixture')
    if (not isinstance(saved, list) or not saved
            or any(not isinstance(row, dict) or set(row) != {'mass', 'activity'}
                   or any(type(row[k]) not in (str, int) for k in ('mass', 'activity'))
                   for row in saved)):
        raise ValueError('nonempty list of exact rational component witnesses required')
    try:
        candidate = [(Q(row['mass']), Q(row['activity'])) for row in saved]
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError('valid exact rational component witnesses required') from error
    prior = {Q(p):Q(c) for c,p in original}
    if (len(prior) != len(original) or len({p for _,p in candidate}) != len(candidate)
            or any(p not in prior or not 0 < c <= prior[p] for c,p in candidate)):
        raise ValueError('distinct original activities and positive nonenlarged masses required')
    shared_mixture.verify(caps, candidate)
    return candidate


def checked_comparison(record):
    """Rebuild default counts or authenticate and replay additional duals."""
    sources = record.get('count_witnesses', [])
    iterations = record.get('count_refinement_iterations', 0)
    if (not isinstance(sources, list) or type(iterations) is not int
            or not 0 <= iterations <= 30 or (not sources and iterations)):
        raise ValueError('count witnesses and a valid exact refinement count required')
    if sources:
        if (record.get('pruned') is not True
                or any(not isinstance(v, dict) or set(v) != {'path', 'sha256'}
                       or not isinstance(v['path'], str) or not v['path']
                       or not isinstance(v['sha256'], str) or len(v['sha256']) != 64
                       or any(c not in '0123456789abcdef' for c in v['sha256']) for v in sources)):
            raise ValueError('pruned comparison with explicit authenticated counting sources required')
        paths = [Path(v['path']) for v in sources]
        if (len({p.resolve() for p in paths}) != len(paths)
                or any(hashlib.sha256(p.read_bytes()).hexdigest() != v['sha256']
                       for p, v in zip(paths, sources))):
            raise ValueError('counting witness duplicate or changed file')
        from shared_relaxed_alternative_counts import build_components
        components, caps, mixture, metadata = build_components(paths,
            iterations=iterations, step=8, zero_bits=record['zero_bits'],
            cost_tilt=Q(record['cost_tilt']), saved_mixture=record['mixture'])
        if ([v['sha256'] for v in metadata['count_witnesses']]
                != [v['sha256'] for v in sources]):
            raise ValueError('counting witness changed while regenerating the comparison')
    else:
        components, caps, mixture = shared_mixture.actual_components(
            coupled=True, zero_bits=record['zero_bits'], cost_tilt=Q(record['cost_tilt']))
        mixture = checked_mixture(record, caps, mixture)
        if record.get('pruned', False):
            components = shared_mixture.as_components(mixture)
    fingerprint = hashlib.sha256(json.dumps([str(c) for c in caps]).encode()).hexdigest()
    if record.get('cap_sha256') != fingerprint:
        raise ValueError('saved cap hash differs from freshly authenticated shared bounds')
    return components


def build_model(record, precision=384, index=0):
    """Regenerate model geometry, comparison coefficients, and inner maps."""
    row = validate_record(record, index)
    if type(precision) is not int or precision < 128:
        raise ValueError('integer precision at least 128 required')
    ctx.prec = precision
    components = checked_comparison(record)
    data = birth_classes.actual(2)
    ctx.prec = precision
    model = sc.Model(components, data, row['threshold'], record['minimum_groups'],
        Q(record['base_tilt']), inner=birth_classes, variance_shuffle=True,
        variance_bins=record['variance_bins'], regional_count=True)
    expected = dict(G=2048, REGIONS=256, PACKETS=1 << 19, EPOCHS=1 << 14, N=1 << 21)
    if (any(getattr(sc, key) != value for key, value in expected.items())
            or any(data.get(key) != value for key, value in dict(bits=19, windows=32, updates=2).items())
            or tuple(map(Q, row['root'])) != model.root):
        raise ValueError('fresh model must have the declared shared-route geometry and root')
    return model


def replay_dense(record, precision=384, index=0):
    require_uncached_verification()
    row = validate_record(record, index, complete=True)
    model = build_model(record, precision, index)
    leaves = row['cover']['leaves']
    cells = sc.partition(model, leaves, {})
    total = Q(0)
    checked = []
    for path in leaves:
        ctx.prec = precision
        upper = model.outward(cells[path], leaves[path]['witness'])
        if ctx.prec != precision or not upper > 0:
            raise ArithmeticError('positive fresh outward bound at requested precision required')
        bound = [int(v) for v in upper.upper().man_exp()]
        total += dyadic(bound)
        checked.append(dict(path=path, upper=bound,
                            output_tilt=str(Q(leaves[path]['witness']['parameters'][0]))))
        print('SHARED DENSE REPLAY', len(checked), 'of', len(leaves), 'path', path, flush=True)
    return dict(complete=True, leaves=len(leaves), upper=compact_endpoint(total), checked=checked,
                threshold=row['threshold'], ensemble=ENSEMBLE,
                aggregation='Exact dyadic sum, then upward rounding to precision+32 mantissa bits.')


def fixed_witness_upper(fresh_dense, threshold, precision=384):
    """Retarget a freshly replayed dense bound without rebuilding matrices.

    This low-level helper requires the result of replay_dense in the current
    proof workflow, not an unauthenticated receipt. For every fixed witness,
    the cutoff D occurs only in exp(lambda*D). All matrix and comparison
    terms are unchanged. The result covers the dense occupancies only.
    """
    if (fresh_dense.get('complete') is not True or fresh_dense.get('ensemble') != ENSEMBLE
            or type(fresh_dense.get('threshold')) is not int
            or not 0 <= fresh_dense['threshold'] < 1 << 21
            or type(threshold) is not int or not 0 <= threshold < 1 << 21
            or type(precision) is not int or precision < 128):
        raise ValueError('fresh same-ensemble replay, valid cutoff, and precision required')
    rows = fresh_dense.get('checked')
    if (not isinstance(rows, list) or not rows or len(rows) != fresh_dense.get('leaves')
            or any(not isinstance(v,dict) or not isinstance(v.get('path'),str) for v in rows)
            or len({v['path'] for v in rows}) != len(rows)):
        raise ValueError('all distinct freshly replayed leaf bounds required')
    sc.partition(SimpleNamespace(root=(Q(0), Q(1))),
                 {v['path']:dict(witness={}) for v in rows}, {})
    if sum((dyadic(v['upper']) for v in rows), Q(0)) > dyadic(fresh_dense['upper']):
        raise ValueError('fresh dense aggregate must enclose the exact sum of checked endpoints')
    ctx.prec = precision
    total = Q(0)
    delta = threshold-fresh_dense['threshold']
    for row in rows:
        original = dyadic(row['upper'])
        tilt = Q(row['output_tilt'])
        if original <= 0 or tilt <= 0:
            raise ValueError('positive fresh endpoint and output tilt required')
        # Arb receives exact integer/rational arguments, never binary floats.
        value = arb(original.numerator)/arb(original.denominator)
        value *= (arb(tilt.numerator)*delta/arb(tilt.denominator)).exp()
        total += dyadic([int(v) for v in value.upper().man_exp()])
    return compact_endpoint(total,precision+32)


def fixed_witness_frontier(fresh_dense, bits=20, precision=384):
    """Largest integer cutoff certified by these fixed dense witnesses.

    This is a dense-only screening frontier. Sparse coverage must be
    regenerated or validly retargeted and added before a whole-code claim.
    """
    if type(bits) is not int or bits < 1:
        raise ValueError('positive integer dense margin required')
    limit = Q(2)**-bits
    if dyadic(fixed_witness_upper(fresh_dense, 0, precision)) >= limit:
        return None
    lo, hi = 0, (1 << 21)-1
    while lo < hi:
        mid = (lo+hi+1)//2
        if dyadic(fixed_witness_upper(fresh_dense, mid, precision)) < limit:
            lo = mid
        else:
            hi = mid-1
    return dict(threshold=lo, minimum_distance=lo+1, dense_bits=bits,
                upper=fixed_witness_upper(fresh_dense, lo, precision),
                note='Fixed-witness dense-only frontier; sparse contribution is not included.')


def sum_sparse(rows, through):
    if (type(through) is not int or not 0 <= through < 2048
            or not isinstance(rows, list) or len(rows) != through
            or any(not isinstance(row, dict) or type(row.get('occupancy')) is not int
                   or type(row.get('updates')) is not int or row['updates'] != 2 for row in rows)
            or sorted(row['occupancy'] for row in rows) != list(range(1, through+1))):
        raise ValueError('exactly one fresh two-update bound for every sparse occupancy required')
    total = Q(0)
    for row in rows:
        upper = dyadic(row.get('upper'))
        if upper <= 0:
            raise ValueError('strictly positive sparse upper endpoint required')
        total += upper
    return total


def combine(dense_result, sparse_rows, through, bits=20):
    if type(bits) is not int or bits < 1:
        raise ValueError('positive integer margin target required')
    if dense_result.get('complete') is not True:
        raise ValueError('fresh complete dense replay required')
    dense_upper = dyadic(dense_result['upper'])
    if dense_upper <= 0:
        raise ValueError('strictly positive dense endpoint required')
    sparse_upper = sum_sparse(sparse_rows, through)
    # Store compact outward component endpoints, then enclose THEIR sum.
    # This preserves receipt consistency even after the final rounding.
    dense_endpoint=compact_endpoint(dense_upper)
    sparse_endpoint=compact_endpoint(sparse_upper)
    total_endpoint=compact_endpoint(dyadic(dense_endpoint)+dyadic(sparse_endpoint))
    if not dyadic(total_endpoint) < Q(2)**-bits:
        raise ValueError('fresh aggregate does not meet the strict margin target')
    return dict(dense_upper=dense_endpoint, sparse_upper=sparse_endpoint,
                total_upper=total_endpoint,
                aggregation='Exact dyadic component sums, each rounded upward to precision+32 '
                            'mantissa bits; total encloses the stored component endpoints.')


def verify(record, precision=384, bits=20, index=0, max_splits=128,
           sparse_output=None, count_witnesses=(), sparse_target_bits=32, dense_replayer=None):
    require_uncached_verification()
    row = validate_record(record, index, complete=True)
    if (type(precision) is not int or precision < 256 or type(bits) is not int or bits < 1
            or type(max_splits) is not int or max_splits < 0
            or type(sparse_target_bits) is not int or sparse_target_bits < 1
            or (dense_replayer is not None and not callable(dense_replayer))):
        raise ValueError('precision >=256, positive margin targets, and nonnegative search budget required')
    dense = (replay_dense if dense_replayer is None else dense_replayer)(record, precision, index)
    if ctx.prec != precision:
        raise ArithmeticError('precision changed during dense replay')
    if dyadic(dense['upper']) >= Q(2)**-bits:
        raise ValueError('fresh dense contribution alone exceeds the margin budget')
    through = record['minimum_groups']-1
    sparse = [] if not through else shared_sparse.run(
        [2], list(range(1, through+1)), row['threshold'], list(TILTS), precision,
        max_splits, sparse_target_bits, output=sparse_output, refined_counts=True,
        coupled_counts=True, joint_return_through=3, lazy_density_through=6,
        count_witnesses=count_witnesses)
    if ctx.prec != precision:
        raise ArithmeticError('precision changed during sparse regeneration')
    result = combine(dense, sparse, through, bits)
    result.update(schema=SCHEMA, ensemble=ENSEMBLE, complete=True,
        precision=precision, requested_bits=bits, message_length=1 << 20,
        output_length=1 << 21, updates=2, threshold=row['threshold'],
        minimum_distance=row['threshold']+1, sparse_occupancies=[1, through],
        dense_occupancies=[through+1, 2048], dense_leaves=dense['leaves'],
        dense_checked=dense['checked'], sparse_checked=sparse,
        count_witnesses=[dict(path=str(p), sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest())
                         for p in count_witnesses],
        note='Fresh shared-route dense replay and sparse regeneration. Saved numerical endpoints '
             'and saved acceptance labels were not proof inputs.')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('dense', type=Path)
    parser.add_argument('--index', type=int, default=0)
    parser.add_argument('--precision', type=int, default=384)
    parser.add_argument('--workers', type=int, choices=range(1,9), default=1,
        help='Fresh dense replay workers; sparse regeneration remains serial')
    parser.add_argument('--bits', type=int, default=20)
    parser.add_argument('--max-splits', type=int, default=128)
    parser.add_argument('--sparse-target-bits', type=int, default=32)
    parser.add_argument('--count-witnesses', nargs='+', type=Path, default=[])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sparse_output = args.output.with_suffix('.sparse.json')
    dense_logs = args.output.with_suffix('.dense-workers')
    if (args.output.exists() or sparse_output.exists()
            or args.dense.resolve() in (args.output.resolve(), sparse_output.resolve())
            or (args.workers > 1 and dense_logs.exists())):
        parser.error('use new output paths; existing proofs are preserved')
    raw = args.dense.read_bytes()
    dense_replayer = None
    if args.workers > 1:
        import shared_relaxed_replay
        def dense_replayer(record, precision, index):
            return shared_relaxed_replay.replay_dense(record, precision, index, args.workers, dense_logs)
    result = verify(json.loads(raw), args.precision, args.bits, args.index, args.max_splits,
                    sparse_output, args.count_witnesses, args.sparse_target_bits, dense_replayer)
    result['dense_replay_workers'] = args.workers
    result['dense_sha256'] = hashlib.sha256(raw).hexdigest()
    if sparse_output.exists():
        result['sparse_receipt_sha256'] = hashlib.sha256(sparse_output.read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    m, e = result['total_upper']
    margin = -(arb(m).log()/arb(2).log()+e)
    print('SHARED R2 COMPLETE: minimum distance', result['minimum_distance'],
          'over', result['output_length'], '; aggregate margin', margin, flush=True)


if __name__ == '__main__':
    main()
