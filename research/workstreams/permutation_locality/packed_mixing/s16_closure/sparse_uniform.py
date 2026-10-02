"""Fresh sparse support bounds for the declared S16 uniform-GL inner.

Only requested occupancies are covered. This file neither imports an S19
bound nor claims a whole-code certificate. The outer count premises are
authenticated afresh; map preparation is delegated to kernel_maps.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from flint import arb, ctx
import sparse_kernel

ROOT = Path(__file__).resolve().parents[2]
PACKED = ROOT/'packed_mixing'
for directory in (ROOT, ROOT/'gf16_packets', PACKED):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

from occupancy_model import placement
from occupancy_memory import rounded
import sparse_hill_variant as old_variant

INNER = dict(t=128, s=16, distribution='fresh independent uniform GL16 per epoch')
KERNEL_VARIANT = 'uniform-all-birth-classes-v1'
SCOPE = dict(
    mixing='independent uniform GL32 per group/block; no additional GF16 stage',
    routing='independent shared column shuffle per group and independent regional shuffles',
    outer_count_refinement='packed-canonical-full32-h5-incidence-ekr-expected-cdf-1')


def legacy_sparse():
    # Use a unique module name for the retained packed sparse interface.
    name = 's16_retained_packed_sparse'
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, PACKED/'sparse.py')
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def build_operators(data, tilts, precision, maximum_groups):
    if (data.get('bits') != 16 or data.get('windows') != 32
            or data.get('distribution') != 'uniform_gl'
            or not tilts or any(not isinstance(t, str) or Q(t) <= 0 for t in tilts)
            or len(set(tilts)) != len(tilts)
            or type(precision) is not int or precision < 128
            or type(maximum_groups) is not int or not 1 <= maximum_groups <= 32):
        raise ValueError('distinct positive tilts, precision and sparse maximum required')
    import kernel_maps
    kernel_maps.authenticate(data)
    result, computed = {}, {}
    for tilt in tilts:
        # Preserve each witness's spelling, but evaluate rational aliases once.
        rational = Q(tilt)
        if rational in computed:
            result[tilt, '1'] = computed[rational]
            continue
        ctx.prec = precision
        local = sparse_kernel.outward(data, tilt)
        exact = placement(local, rounding=rounded, maximum_groups=maximum_groups)
        n = local[0].nrows()
        arrays = [np.array([[float(m[i,j]) for j in range(n)] for i in range(n)]) for m in exact]
        result[tilt, '1'] = exact, arrays
        computed[rational] = result[tilt, '1']
        print('S16 UNIFORM fixed-occupancy operators', tilt, 'states', n,
              'degree', maximum_groups, flush=True)
    return result


def run(occupancies, tilts, *, feedback='weight5', distance='.10', precision=256,
        target_bits=46, max_splits=64, output=None):
    if (not occupancies or len(set(occupancies)) != len(occupancies)
            or any(type(q) is not int or not 1 <= q <= 32 for q in occupancies)
            or not 0 < Q(distance) < Q(1,2)
            or type(target_bits) is not int or target_bits < 1
            or type(max_splits) is not int or max_splits < 0
            or output is not None and output.exists()):
        raise ValueError('distinct sparse occupancies, valid limits and fresh output required')
    import kernel_maps
    data, map_record = kernel_maps.prepare(feedback=feedback, seed=0, refresh='uniform', birth_density='classes')
    map_record = json.loads(json.dumps(map_record))
    sparse = legacy_sparse()
    counts, floor, premises = old_variant.authenticated_counts(sparse, 'ekr')
    threshold = (Q(distance)*(1 << 21)).__floor__()
    record = dict(schema='s16-uniform-gl-sparse-search-2', K=1 << 20, N=1 << 21,
        group_count=2048, outer='BCH256128', block_rows=4, block_columns=8,
        inner=INNER.copy(), kernel_variant=KERNEL_VARIANT,
        feedback=feedback, map_record=map_record, distance=str(Q(distance)), threshold=threshold,
        precision=precision, target_bits=target_bits, max_splits=max_splits,
        count_sha256=sparse.count_hash(counts), count_premises=premises, support_min=floor,
        requested=list(occupancies), tilts=list(tilts), results=[], complete_requested=False,
        complete_sparse_prefix=False, whole_code_certificate=False, **SCOPE)
    sparse.save(output, record)
    operators = build_operators(data, tilts, precision, max(occupancies))
    terminal = np.ones(next(iter(operators.values()))[0][0].nrows())
    args = SimpleNamespace(groups=0, precision=precision, target_bits=target_bits,
        max_splits=max_splits, joint_witness=True, joint_top=2)
    total = arb(0)
    for q in occupancies:
        args.groups = q
        ctx.prec = precision
        if q == 1:
            upper, details = old_variant.one_group(sparse, counts, operators, threshold, precision)
        else:
            from support_cover import cover
            upper, details = cover(args, operators, counts, terminal, cutoff=threshold, support_min=floor)
        passed = bool(upper is not None and upper.is_finite() and 0 < upper < arb(2)**(-target_bits))
        record['results'].append(dict(occupancy=q, passed=passed,
            upper=sparse.endpoint(upper) if upper is not None else None, details=details))
        if passed:
            total = sparse.up(total+upper)
        record['complete_requested'] = (len(record['results']) == len(occupancies)
                                       and all(r['passed'] for r in record['results']))
        record['complete_sparse_prefix'] = (record['complete_requested']
                                            and set(occupancies) == set(range(1,33)))
        record['successful_aggregate_upper'] = sparse.endpoint(total) if total > 0 else None
        record['aggregate_upper'] = sparse.endpoint(total) if record['complete_requested'] else None
        record['aggregate_target_passed'] = bool(record['complete_requested']
                                                and total < arb(2)**(-target_bits))
        sparse.save(output, record)
        print('S16 UNIFORM sparse', feedback, 'q', q, 'passed', passed, 'log2 upper',
              upper.log()/arb(2).log() if upper is not None else 'split budget exhausted', flush=True)
    return record


def validate_record(record):
    """Select search witnesses, never numerical endpoints, in the fixed scope."""
    schema = record.get('schema')
    # Version1's only implementation fixed these three fields implicitly.
    # Its witnesses may be replayed, but all new receipts state them explicitly.
    scope = {name: record.get(name, value if schema == 's16-uniform-gl-sparse-search-1' else None)
             for name, value in SCOPE.items()}
    if (schema not in ('s16-uniform-gl-sparse-search-1', 's16-uniform-gl-sparse-search-2')
            or scope != SCOPE
            or record.get('K') != 1 << 20 or record.get('N') != 1 << 21
            or record.get('group_count') != 2048 or record.get('outer') != 'BCH256128'
            or record.get('block_rows') != 4 or record.get('block_columns') != 8
            or record.get('inner') != INNER or record.get('feedback') not in ('weight5', 'bch16')
            or record.get('kernel_variant', KERNEL_VARIANT) != KERNEL_VARIANT
            or not 0 < Q(record['distance']) < Q(1,2)
            or type(record.get('threshold')) is not int
            or record['threshold'] != (Q(record['distance'])*(1 << 21)).__floor__()):
        raise ValueError('matching S16 uniform-GL sparse search scope required')
    requested, rows, tilts = record.get('requested'), record.get('results'), record.get('tilts')
    if (not isinstance(requested, list) or not requested
            or len(set(requested)) != len(requested)
            or any(type(q) is not int or not 1 <= q <= 32 for q in requested)
            or not isinstance(rows, list) or len({r['occupancy'] for r in rows}) != len(rows)
            or any(type(r['occupancy']) is not int or r['occupancy'] not in requested
                   or type(r.get('passed')) is not bool for r in rows)
            or not isinstance(tilts, list) or not tilts or len(set(tilts)) != len(tilts)
            or any(not isinstance(t, str) or Q(t) <= 0 for t in tilts)):
        raise ValueError('distinct occupancy witnesses and exact positive tilts required')
    selected = [row for row in rows if row['passed']]
    for row in selected:
        details = row['details']
        if details.get('support_min') != record.get('support_min'):
            raise ValueError('complete support domain required')
        if row['occupancy'] == 1:
            choices = details.get('support_choices')
            if (details.get('method') != 'all-support exact placement'
                    or not isinstance(choices, list) or len(choices) != 257
                    or any(t is not None and t not in tilts for t in choices)):
                raise ValueError('full q1 support witnesses required')
        elif (details.get('method') != 'full support-box CDF cover' or not details.get('leaves')
              or any(leaf['tilt'] not in tilts for leaf in details['leaves'])):
            raise ValueError('declared complete support-box witness required')
    return selected


def replay(sources, output, *, precision=384, target_bits=46, distance=None):
    """Authenticate maps/counts, rebuild operators, and check support witnesses."""
    if (not sources or output.exists() or type(precision) is not int or precision < 128
            or type(target_bits) is not int or target_bits < 1):
        raise ValueError('source receipts, fresh output and valid exact limits required')
    records, entries, metadata, seen, tilts = [], [], [], set(), set()
    for source in sources:
        payload = source.read_bytes()
        record = json.loads(payload)
        selected = validate_record(record)
        for row in selected:
            q = row['occupancy']
            if q in seen:
                raise ValueError('duplicate successful occupancy across receipts')
            seen.add(q)
            entries.append(row)
        tilts.update(record['tilts'])
        records.append(record)
        metadata.append(dict(path=str(source.resolve()), sha256=hashlib.sha256(payload).hexdigest()))
    if not entries or len({record['feedback'] for record in records}) != 1:
        raise ValueError('successful witnesses for exactly one feedback map required')
    distance = min(Q(r['distance']) for r in records) if distance is None else Q(distance)
    if not 0 < distance < Q(1,2) or any(distance > Q(r['distance']) for r in records):
        raise ValueError('only the same or lower distance may be replayed')
    import kernel_maps
    feedback = records[0]['feedback']
    data, map_record = kernel_maps.prepare(feedback=feedback, seed=0, refresh='uniform', birth_density='classes')
    map_record = json.loads(json.dumps(map_record))
    sparse = legacy_sparse()
    counts, floor, premises = old_variant.authenticated_counts(sparse, 'ekr')
    digest = sparse.count_hash(counts)
    if any(r['map_record'] != map_record or r['count_premises'] != premises
           or r['count_sha256'] != digest or r['support_min'] != floor for r in records):
        raise ValueError('saved map/count metadata differ from freshly authenticated premises')
    threshold = (distance*(1 << 21)).__floor__()
    tilts = sorted(tilts, key=Q)
    operators = build_operators(data, tilts, precision, max(seen))
    terminal = np.ones(next(iter(operators.values()))[0][0].nrows())
    result = dict(schema='s16-uniform-gl-sparse-replay-2', K=1 << 20, N=1 << 21,
        inner=INNER.copy(), kernel_variant=KERNEL_VARIANT,
        feedback=feedback, map_record=map_record, precision=precision,
        distance=str(distance), threshold=threshold, target_bits=target_bits,
        count_sha256=digest, count_premises=premises, support_min=floor,
        sources=metadata, requested=sorted(seen), tilts=tilts, results=[],
        complete_requested=False, complete_sparse_prefix=False, whole_code_certificate=False, **SCOPE)
    sparse.save(output, result)
    ctx.prec = precision
    total = arb(0)
    for entry in sorted(entries, key=lambda row: row['occupancy']):
        q = entry['occupancy']
        if q == 1:
            upper, _ = old_variant.one_group(sparse, counts, operators, threshold, precision)
        else:
            from support_cover import replay as support_replay
            upper = support_replay(entry['details'], operators, counts, terminal,
                groups=q, cutoff=threshold, precision=precision, target_bits=target_bits)
        if not upper.is_finite() or not 0 < upper < arb(2)**(-target_bits) or ctx.prec != precision:
            raise ArithmeticError('fresh sparse S16 replay failed')
        total = sparse.up(total+upper)
        result['results'].append(dict(occupancy=q, passed=True, upper=sparse.endpoint(upper)))
        result['aggregate_upper'] = sparse.endpoint(total)
        result['complete_requested'] = len(result['results']) == len(entries)
        result['complete_sparse_prefix'] = result['complete_requested'] and seen == set(range(1,33))
        result['aggregate_target_passed'] = bool(result['complete_requested'] and total < arb(2)**(-target_bits))
        if result['complete_requested']:
            for item in metadata:
                if hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() != item['sha256']:
                    raise ValueError('source receipt changed during replay')
        sparse.save(output, result)
        print('FRESH S16 SPARSE', feedback, 'q', q, 'margin', -upper.log()/arb(2).log(), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occupancies', type=int, nargs='+')
    parser.add_argument('--tilts', nargs='+')
    parser.add_argument('--feedback', choices=('weight5', 'bch16'))
    parser.add_argument('--distance')
    parser.add_argument('--precision', type=int)
    parser.add_argument('--target-bits', type=int, default=46)
    parser.add_argument('--max-splits', type=int)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--replay', type=Path, nargs='+')
    args = parser.parse_args()
    if args.replay:
        if any(getattr(args, option) is not None for option in ('occupancies', 'tilts', 'feedback', 'max_splits')):
            parser.error('replay rejects search occupancy, tilt, feedback and split options')
        replay(args.replay, args.output, precision=384 if args.precision is None else args.precision,
               target_bits=args.target_bits, distance=args.distance)
    else:
        run([1] if args.occupancies is None else args.occupancies,
            ['.00032', '.0016', '.008', '.032'] if args.tilts is None else args.tilts,
            feedback='weight5' if args.feedback is None else args.feedback,
            distance='.10' if args.distance is None else args.distance,
            precision=256 if args.precision is None else args.precision,
            target_bits=args.target_bits, max_splits=64 if args.max_splits is None else args.max_splits,
            output=args.output)


if __name__ == '__main__':
    main()
