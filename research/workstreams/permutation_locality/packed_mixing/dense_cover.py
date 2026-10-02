"""Resumable dense coverage for the canonical GL32 packet ensemble.

Every invocation authenticates BCH premises and reconstructs the expected
outer CDF. Saved coefficients are merely positive-majorant candidates and
saved cell witnesses are replayed; saved numeric bounds are never accepted
as proof inputs. A complete dense cover still excludes the sparse prefix.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

from canonical_counts import authenticated_bch_cdf

PARENT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PARENT))
sys.path.insert(0, str(PARENT / 'gf16_packets'))

SCHEMA = 'packed-canonical-gl32-dense-cover-1'
ENSEMBLE = 'canonical-gl32-width8-shared4-r2'


def cutoff(distance):
    distance = Q(distance)
    if not 0 < distance < Q(1, 2):
        raise ValueError('distance must be in (0,1/2)')
    return int(distance * (1 << 21))


def fingerprint(caps):
    return hashlib.sha256(json.dumps([str(v) for v in caps]).encode()).hexdigest()


def validate_record(record):
    if (record.get('schema') != SCHEMA or record.get('ensemble') != ENSEMBLE
            or record.get('K') != 1 << 20 or record.get('N') != 1 << 21
            or record.get('updates') != 2 or record.get('block_width') != 8
            or type(record.get('minimum_groups')) is not int
            or not 1 <= record['minimum_groups'] <= 2048
            or record.get('maximum_groups') != 2048
            or record.get('threshold') != cutoff(record['distance'])
            or record.get('comparison') not in ('expected-cdf-shell-majorant', 'direct-expected-shell-majorant')
            or record.get('last_lp') != 104 or record.get('refined') is not True
            or Q(record.get('base_tilt')) != Q(3, 16)
            or record.get('variance_bins') != 16
            or record.get('regional_count') is not True):
        raise ValueError('matching canonical GL32 construction and proof scope required')
    cover = record.get('cover')
    if not isinstance(cover, dict) or not isinstance(cover.get('leaves'), dict) or not isinstance(cover.get('unresolved'), dict):
        raise ValueError('saved dense partition required')
    return cover


def exact_mixture(caps, saved=None):
    import shared_mixture
    from positive_prune import prune
    if saved is not None:
        if not isinstance(saved, list) or not saved:
            raise ValueError('positive rational comparison candidate required')
        mixture = [(Q(row['mass']), Q(row['activity'])) for row in saved]
        if len({p for _, p in mixture}) != len(mixture):
            raise ValueError('distinct comparison centers required')
        shared_mixture.verify(caps, mixture)
        return mixture, dict(source='saved candidate; every shell freshly verified')
    centers = sorted({Q(u, 256) for u in range(1, 257, 4)} | {Q(1)})
    mixture = shared_mixture.envelope(caps, centers, zero_bits=64, cost_tilt=Q(1, 4))
    mixture, info = prune(caps, mixture, Q(1, 4))
    shared_mixture.verify(caps, mixture)
    return mixture, info


def dyadic(pair):
    if (not isinstance(pair, list) or len(pair) != 2
            or any(type(v) is not int for v in pair) or pair[0] < 0):
        raise ValueError('nonnegative dyadic endpoint required')
    return Q(pair[0]) * Q(2)**pair[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance', default='.095')
    parser.add_argument('--minimum-groups', type=int, default=33)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--target-bits', type=int, default=36)
    parser.add_argument('--max-cells', type=int, default=120)
    parser.add_argument('--max-depth', type=int, default=22)
    parser.add_argument('--means', nargs='*', default=[])
    parser.add_argument('--comparison', choices=('cdf', 'shell'), default='shell')
    parser.add_argument('--regional-presplit', action='store_true',
        help='Geometrically split unresolved cells to width 1/1024 before numerical search')
    parser.add_argument('--point-seeds', type=Path, nargs='+', default=[],
        help='Use old point witnesses as proposals; never transfer saved bounds')
    parser.add_argument('--region-cache', action='store_true',
        help='Reuse a freshly built cell-independent inner operator during seeded search')
    parser.add_argument('--seed-width', type=Q,
        help='Retile unresolved cells to this width; failed wide hints split without optimization')
    parser.add_argument('--replay-workers', type=int, default=1,
        help='Use 1..4 fresh uncached processes for independent replay only')
    prior = parser.add_mutually_exclusive_group()
    prior.add_argument('--resume', type=Path)
    prior.add_argument('--replay', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if (args.output.exists() or args.precision < 128 or args.target_bits < 1
            or args.max_cells < 0 or args.max_depth < 1
            or not 1 <= args.minimum_groups <= 2048
            or any(not 0 <= Q(x) <= 1 for x in args.means)):
        parser.error('fresh output path, valid scope, precision, and work limits required')
    if args.replay and (args.regional_presplit or args.point_seeds or args.region_cache
            or args.seed_width is not None):
        parser.error('search scheduling and caches are not enabled during independent replay')
    if args.region_cache and not args.point_seeds:
        parser.error('the current operator cache requires point-seeded search')
    if not 1 <= args.replay_workers <= 4 or (args.replay_workers != 1 and not args.replay):
        parser.error('replay workers must be 1..4 and parallel workers require --replay')
    if args.seed_width is not None and (not args.point_seeds or args.regional_presplit
            or not Q(1, 1024) <= args.seed_width <= 1):
        parser.error('seed width requires point seeds, width 1/1024..1, and no presplit flag')

    saved = None
    source_bytes = None
    source = args.resume or args.replay
    if source:
        source_bytes = source.read_bytes()
        saved = json.loads(source_bytes)
        validate_record(saved)
        args.distance, args.minimum_groups = saved['distance'], saved['minimum_groups']
        args.comparison = 'shell' if saved['comparison'] == 'direct-expected-shell-majorant' else 'cdf'
    threshold = cutoff(args.distance)
    caps, premises = authenticated_bch_cdf()
    digest = fingerprint(caps)
    if saved and saved.get('expected_cdf_sha256') != digest:
        raise ValueError('fresh outer bounds differ from saved comparison premises')
    comparison_caps = caps
    if args.comparison == 'shell':
        from monotone import transport_shells
        from local_models import full_block
        comparison_caps = transport_shells(premises['canonical_cdf'], full_block(8))
    comparison_digest = fingerprint(comparison_caps)
    if saved and saved.get('comparison_caps_sha256') != comparison_digest:
        raise ValueError('fresh comparison caps differ from saved premises')
    mixture, info = exact_mixture(comparison_caps, saved.get('mixture') if saved else None)

    import birth_classes
    import shared_mixture
    import scalar_cover as sc
    from flint import arb, ctx
    data = birth_classes.actual(2)
    ctx.prec = args.precision
    model = sc.Model(shared_mixture.as_components(mixture), data, threshold,
        args.minimum_groups, Q(3, 16), inner=birth_classes,
        variance_shuffle=True, variance_bins=16, regional_count=True)
    model.proposal_stop_bits = args.target_bits + 2
    point_model, regional_cache, hint_sources = None, None, []
    if args.point_seeds:
        from point_reuse import load_points, PointReuseModel
        probes, hint_sources = load_points(args.point_seeds)
        evaluator, regional_warm = None, None
        if args.region_cache:
            from region_cache import RegionalCache
            regional_cache = RegionalCache()
            base_model = model
            evaluator = lambda cell, witness: regional_cache.bound(base_model, cell, witness)
            regional_warm = lambda witness: regional_cache.is_warm(base_model, witness)
        point_model = PointReuseModel(model, probes, args.target_bits, outward=evaluator,
            regional_warm=regional_warm,
            maximum_optimization_width=Q(1, 1024) if args.seed_width is not None else None)
        model = point_model
    record = dict(schema=SCHEMA, ensemble=ENSEMBLE, K=1 << 20, N=1 << 21,
        updates=2, block_width=8, distance=str(Q(args.distance)), threshold=threshold,
        minimum_groups=args.minimum_groups, maximum_groups=2048,
        root=list(map(str, model.root)), precision=args.precision,
        target_bits=args.target_bits, max_cells=args.max_cells, max_depth=args.max_depth,
        comparison=('direct-expected-shell-majorant' if args.comparison == 'shell' else 'expected-cdf-shell-majorant'),
        expected_cdf_sha256=digest, comparison_caps_sha256=comparison_digest,
        last_lp=104, refined=True, base_tilt='3/16', variance_bins=16,
        regional_count=True, outer_premises=premises,
        mixture=[dict(mass=str(c), activity=str(p)) for c, p in mixture],
        mixture_verification=info, probes=[], cover=dict(leaves={}, unresolved={'': {}}),
        proof_status='Dense-only attempt; sparse prefix excluded. No whole-code claim.')
    if hint_sources:
        record['point_hint_sources'] = hint_sources
    if args.seed_width is not None:
        record['seed_width'] = str(args.seed_width)
        record['maximum_optimization_width'] = '1/1024'
    if source:
        if saved.get('root') != record['root']:
            raise ValueError('saved and freshly derived comparison domains differ')
        record['source'] = dict(path=str(source.resolve()), sha256=hashlib.sha256(source_bytes).hexdigest())

    resume_cover = saved['cover'] if saved else None
    if args.regional_presplit:
        if args.replay:
            parser.error('geometric preprocessing is a search option, not a replay option')
        from search_strategy import presplit
        resume_cover, info = presplit(model.root,
            resume_cover or record['cover'], Q(1, 1024), args.max_depth)
        record['search_preprocessing'] = info
        print('PACKED PRESPLIT', info, flush=True)
    if args.seed_width is not None:
        from search_strategy import retile
        resume_cover, info = retile(model.root,
            resume_cover or record['cover'], args.seed_width, args.max_depth)
        record['search_preprocessing'] = info
        print('PACKED RETILE', info, flush=True)
    if resume_cover is not None:
        record['cover'] = resume_cover

    def save(state=None):
        if state is not None:
            record['cover'] = state
        if point_model is not None:
            record['point_reuse_stats'] = dict(point_model.stats)
        if regional_cache is not None:
            record['regional_cache_stats'] = dict(hits=regional_cache.hits, misses=regional_cache.misses)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2) + '\n')

    save()
    print('PACKED DENSE domain', record['root'], 'q', args.minimum_groups,
          '..2048; exact comparison components', len(mixture), flush=True)
    if args.replay:
        cells = sc.partition(model, saved['cover']['leaves'], saved['cover']['unresolved'])
        total, checked = Q(0), []
        if args.replay_workers > 1 and saved['cover']['leaves']:
            from replay_pool import replay_dense
            tasks = [(path, cells[path], row['witness'])
                     for path, row in saved['cover']['leaves'].items()]
            logdir = args.output.with_name(args.output.stem+'-worker-logs')
            print('PACKED REPLAY starting', args.replay_workers, 'fresh workers for', len(tasks), 'cells', flush=True)
            checked = replay_dense(record['mixture'], threshold, args.minimum_groups,
                args.precision, tasks, workers=args.replay_workers, logdir=logdir)
            total = sum((dyadic(row['upper']) for row in checked), Q(0))
        else:
            for path, row in saved['cover']['leaves'].items():
                upper = model.outward(cells[path], row['witness'])
                if not upper.is_finite() or not upper > 0:
                    raise ArithmeticError('finite positive fresh cell endpoint required')
                endpoint = [int(v) for v in upper.upper().man_exp()]
                total += dyadic(endpoint)
                checked.append(dict(path=path, cell=list(map(str, cells[path])), upper=endpoint))
                print('PACKED REPLAY', len(checked), '/', len(saved['cover']['leaves']), flush=True)
        bounded = arb(total.numerator) / arb(total.denominator)
        record['cover'] = saved['cover']
        record['fresh_replay'] = dict(checked=checked, unresolved=len(saved['cover']['unresolved']),
            upper=[int(v) for v in bounded.upper().man_exp()],
            log2_upper=str(bounded.log()/arb(2).log()) if total else '-inf',
            complete_dense=not saved['cover']['unresolved'], maximum_workers=args.replay_workers)
        save()
        return
    for mean in map(Q, args.means):
        cell = (mean, mean)
        if model.empty(cell):
            record['probes'].append(dict(mean=str(mean), empty=True))
        else:
            score, witness = model.proposal(cell)
            upper = model.outward(cell, witness)
            row = dict(mean=str(mean), proposal=float(score), witness=witness,
                upper=[int(v) for v in upper.upper().man_exp()],
                log2_upper=str(upper.log()/arb(2).log()))
            record['probes'].append(row)
            print('PACKED POINT', mean, row['log2_upper'], flush=True)
        save()
    if args.max_cells:
        result = sc.run(model, args.max_cells, args.max_depth, args.target_bits,
            checkpoint=save, resume=resume_cover, coalesce=False)
        sc.partition(model, result['leaves'], result['unresolved'])
        save(result)
        print('PACKED DENSE covered', len(result['leaves']), 'unresolved', len(result['unresolved']), flush=True)
    print('No whole-code certificate; fresh dense replay and sparse aggregation are required.', flush=True)


if __name__ == '__main__':
    main()
