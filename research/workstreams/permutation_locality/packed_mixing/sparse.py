"""Fresh sparse first moments for canonical GL32 mixing; never a full certificate.

For each fixed nonzero four-row outer input, independent local GL32 maps
produce the expected support CDF authenticated here. The later independent
uniform column permutation and region permutations supply the fixed-support
inner interface. Counts are rational expectations, not integer shell counts.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'gf16_packets'))

import numpy as np
from flint import arb, ctx
import occupancy_birth_classes
import occupancy_rank
import single_group
import sparse_cover as old_args
from canonical_counts import authenticated_bch_cdf

aq, up = occupancy_rank.aq, occupancy_rank.up
N = 1 << 21


def checked_counts(counts, support_min=None):
    """Do not omit any positive expected mass, including sub-BCH-distance tails."""
    counts = tuple(map(Q, counts))
    if (len(counts) != 257 or counts[0] != 0 or counts[-1] != (1 << 512)-1
            or any(a < 0 or a > b for a, b in zip(counts, counts[1:]))):
        raise ValueError('complete nonzero-message cumulative expected counts required')
    first = next(u for u, c in enumerate(counts) if c)
    if support_min is not None and (type(support_min) is not int or support_min != first):
        raise ValueError('support floor must equal the first positive expected CDF entry')
    return counts, first


def fold_cdf(counts, weights):
    """Directed Abel fold with rational cumulative caps and a suffix maximum."""
    if (len(counts) != len(weights) or not counts or counts[0] != 0
            or any(c < 0 for c in counts)
            or any(a > b for a, b in zip(counts, counts[1:]))
            or any(w < 0 for w in weights)):
        raise ValueError('nonnegative cumulative counts and weights required')
    tails = [arb(0)]*len(weights)
    for u in range(len(weights)-1, -1, -1):
        tails[u] = max(weights[u], tails[u+1] if u+1 < len(weights) else arb(0))
    return up(sum((aq(Q(counts[u]-counts[u-1]))*tails[u]
                   for u in range(1, len(counts))), arb(0)))


def endpoint(value):
    return [int(x) for x in up(value).man_exp()]


def count_hash(counts):
    return hashlib.sha256(json.dumps([str(c) for c in counts], separators=(',', ':')).encode()).hexdigest()


def save(path, record):
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(record, indent=2)+'\n')


def one_group(counts, operators, threshold, precision):
    """Exact placement averaging for every support 0..256, including 5..37."""
    checked_counts(counts, 5)
    best = [arb(1)]*257
    choices = [None]*257
    for (tilt, penalty), (exact, _) in operators.items():
        if penalty != '1':
            raise ValueError('this proof uses unweighted counts')
        ctx.prec = precision
        moments = single_group.support_moments(exact[0], exact[1])
        factor = (aq(Q(tilt))*threshold).exp()
        for u, moment in enumerate(moments):
            bound = up(factor*moment)
            if bound < best[u]:
                best[u], choices[u] = bound, tilt
        upper = up(2048*fold_cdf(counts, best))
        print('PACKED q1 tilt', tilt, 'log2 upper', upper.log()/arb(2).log(), flush=True)
    return upper, dict(method='all-support exact placement', support_min=5,
        support_choices=choices, support_probability_uppers=[endpoint(v) for v in best])


def validate_record(record):
    if (record.get('schema') != 'packed-gl32-sparse-1'
            or record.get('K') != 1 << 20 or record.get('N') != N
            or record.get('group_count') != 2048 or record.get('outer') != 'BCH256128'
            or record.get('block_rows') != 4 or record.get('block_columns') != 8
            or record.get('mixing') != 'independent uniform GL32 per group/block; no additional GF16 stage'
            or record.get('routing') != 'independent shared column shuffle per group and independent regional shuffles'
            or record.get('inner') != dict(t=128, s=19, updates=2)
            or record.get('support_min') != 5
            or not 0 < Q(record['distance']) < Q(1, 2)
            or record.get('threshold') != (Q(record['distance'])*N).__floor__()):
        raise ValueError('matching canonical GL32 construction and exact threshold required')
    entries = record['results']
    if (not isinstance(entries, list) or len({x['occupancy'] for x in entries}) != len(entries)
            or any(type(x['occupancy']) is not int or not 1 <= x['occupancy'] <= 32 for x in entries)):
        raise ValueError('distinct sparse occupancy results required')
    return [x for x in entries if x['passed'] is True and x['upper'] is not None]


def replay_records(sources, output, precision=384, target_bits=24, distance=None):
    """Freshly authenticate counts and operators; saved numbers are not proof inputs."""
    if output.exists() or type(precision) is not int or precision < 128 or type(target_bits) is not int or target_bits < 1:
        raise ValueError('fresh output and valid precision/target required')
    records, entries, seen, tilts, source_metadata = [], [], set(), set(), []
    for source in sources:
        source_bytes = source.read_bytes()
        record = json.loads(source_bytes)
        source_metadata.append(dict(path=str(source.resolve()), sha256=hashlib.sha256(source_bytes).hexdigest()))
        candidates = validate_record(record)
        for entry in candidates:
            q = entry['occupancy']
            if q in seen:
                raise ValueError('duplicate verified occupancy across input records')
            seen.add(q)
            if q == 1:
                tilts.update(x for x in entry['details']['support_choices'] if x is not None)
            else:
                tilts.update(x['tilt'] for x in entry['details']['leaves'])
            entries.append(entry)
        records.append(record)
    if not entries:
        raise ValueError('at least one previously covered sparse occupancy required')
    distance = Q(distance) if distance is not None else min(Q(r['distance']) for r in records)
    if not 0 < distance < Q(1, 2) or any(distance > Q(r['distance']) for r in records):
        raise ValueError('replay supports only the same or a lower claimed distance')
    threshold = (distance*N).__floor__()
    counts, premises = authenticated_bch_cdf(104, True, 8)
    counts, _ = checked_counts(counts, 5)
    digest = count_hash(counts)
    if any(r['count_sha256'] != digest or r['count_premises'] != premises for r in records):
        raise ValueError('saved count metadata does not match fresh authenticated premises')
    tilts = sorted(tilts, key=Q)
    args = old_args.build_args(max(seen), tilts, precision, 0, target_bits, None, 2)
    args.exact_feedback = True
    args.joint_return_through = max(r['joint_return_through'] for r in records)
    args.lazy_density_through = max(r['lazy_density_through'] for r in records)
    operators = occupancy_birth_classes.build_operators(args)
    terminal = np.ones(next(iter(operators.values()))[0][0].nrows())
    result = dict(schema='packed-gl32-sparse-replay-1', precision=precision,
        distance=str(distance), threshold=threshold, target_bits=target_bits,
        count_sha256=digest, count_premises=premises,
        scope={k: records[0][k] for k in ('K', 'N', 'outer', 'group_count', 'block_rows',
            'block_columns', 'mixing', 'routing', 'inner', 'support_min')},
        sources=source_metadata,
        tilts=tilts, joint_return_through=args.joint_return_through,
        lazy_density_through=args.lazy_density_through, requested=sorted(seen),
        uncovered_sparse=sorted(set(range(1, 33))-seen), results=[], complete_requested=False,
        complete_sparse_prefix=False,
        note='Aggregate endpoint, not per-occupancy targets, determines the total margin. Dense occupancies remain outside this receipt.')
    save(output, result)
    total = arb(0)
    for entry in sorted(entries, key=lambda x: x['occupancy']):
        q = entry['occupancy']
        ctx.prec = precision
        if q == 1:
            upper, _ = one_group(counts, operators, threshold, precision)
        else:
            from support_cover import replay
            upper = replay(entry['details'], operators, counts, terminal, groups=q,
                cutoff=threshold, precision=precision, target_bits=target_bits)
        if not 0 < upper < arb(2)**(-target_bits):
            raise ArithmeticError('fresh sparse replay failed')
        total = up(total+upper)
        result['results'].append(dict(occupancy=q, upper=endpoint(upper)))
        result['aggregate_upper'] = endpoint(total)
        result['complete_requested'] = len(result['results']) == len(entries)
        result['complete_sparse_prefix'] = result['complete_requested'] and not result['uncovered_sparse']
        save(output, result)
        print('FRESH packed sparse replay q', q, 'margin', -upper.log()/arb(2).log(), flush=True)
    return result


def run(occupancies, distance, tilts, precision=256, target_bits=32,
        max_splits=100, joint_return_through=1, lazy_density_through=1, output=None):
    if output is not None and output.exists():
        raise ValueError('new output path required; earlier receipts are preserved')
    if (not occupancies or len(set(occupancies)) != len(occupancies)
            or any(type(q) is not int or not 1 <= q <= 32 for q in occupancies)
            or not 0 < Q(distance) < Q(1, 2) or not tilts or any(Q(t) <= 0 for t in tilts)
            or type(precision) is not int or precision < 128
            or type(target_bits) is not int or target_bits < 1
            or type(max_splits) is not int or max_splits < 0):
        raise ValueError('valid sparse occupancies, exact distance, precision and limits required')
    threshold = (Q(distance)*N).__floor__()
    counts, premises = authenticated_bch_cdf(104, True, 8)
    counts, floor = checked_counts(counts, 5)
    record = dict(schema='packed-gl32-sparse-1', K=1 << 20, N=N, group_count=2048,
        outer='BCH256128', block_rows=4, block_columns=8,
        mixing='independent uniform GL32 per group/block; no additional GF16 stage',
        routing='independent shared column shuffle per group and independent regional shuffles',
        inner=dict(t=128, s=19, updates=2), distance=str(Q(distance)), threshold=threshold,
        precision=precision, target_bits=target_bits, tilts=list(tilts), max_splits=max_splits,
        joint_return_through=joint_return_through, lazy_density_through=lazy_density_through,
        count_sha256=count_hash(counts), count_premises=premises, support_min=floor,
        requested=list(occupancies), results=[],
        note='Only listed sparse occupancies. Freshly regenerate counts and operators to replay; no full-code certificate.')
    save(output, record)
    args = old_args.build_args(max(occupancies), tilts, precision, max_splits, target_bits, None, 2)
    args.exact_feedback = True
    args.analytic_gradient = True
    args.joint_return_through = joint_return_through
    args.lazy_density_through = lazy_density_through
    operators = occupancy_birth_classes.build_operators(args)
    terminal = np.ones(next(iter(operators.values()))[0][0].nrows())
    total = arb(0)
    for q in occupancies:
        args.groups = q
        ctx.prec = precision
        if q == 1:
            upper, details = one_group(counts, operators, threshold, precision)
        else:
            from support_cover import cover
            upper, details = cover(args, operators, counts, terminal, cutoff=threshold, support_min=floor)
        passed = upper is not None and upper < arb(2)**(-target_bits)
        entry = dict(occupancy=q, passed=passed, upper=endpoint(upper) if upper is not None else None, details=details)
        record['results'].append(entry)
        if upper is not None:
            total = up(total+upper)
            print('PACKED sparse q', q, 'log2 upper', upper.log()/arb(2).log(), 'target passed', passed, flush=True)
        record['aggregate_upper'] = endpoint(total) if total > 0 else None
        record['complete_requested'] = len(record['results']) == len(occupancies) and all(x['passed'] for x in record['results'])
        save(output, record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--occupancies', type=int, nargs='+', default=[1])
    parser.add_argument('--distance', default='.095')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--target-bits', type=int, default=32)
    parser.add_argument('--max-splits', type=int, default=100)
    parser.add_argument('--joint-return-through', type=int, default=1)
    parser.add_argument('--lazy-density-through', type=int, default=1)
    parser.add_argument('--tilts', nargs='+', default=['.00032', '.001', '.0032', '.008', '.016', '.032', '.064', '.096'])
    parser.add_argument('--output', type=Path)
    parser.add_argument('--replay', type=Path, nargs='+', help='Freshly replay all passed entries from these disjoint-scope receipts')
    args = parser.parse_args()
    if args.replay:
        if args.output is None:
            parser.error('replay requires a fresh output path')
        replay_records(args.replay, args.output, args.precision, args.target_bits, args.distance)
    else:
        run(args.occupancies, args.distance, args.tilts, args.precision, args.target_bits,
            args.max_splits, args.joint_return_through, args.lazy_density_through, args.output)


if __name__ == '__main__':
    main()
