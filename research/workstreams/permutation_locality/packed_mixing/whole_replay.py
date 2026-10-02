"""Fresh whole-code replay for canonical GL32, shared routing, and R2.

Inputs are one complete dense witness cover and sparse witness receipts
whose passed occupancies partition 1..32. Their saved numerical bounds are
not proof inputs. The two existing replay programs run serially in fresh
processes and write into a new directory. Dense replay may itself use
independent uncached workers. Only their new per-cell and
per-occupancy dyadic endpoints enter the exact final sum.

Success certifies an expected number of nonzero messages of weight at most
floor(distance*2^21) below 2^-target_bits for the declared ideal setup law.
Markov's inequality gives the same bound on failure of this distance claim.
This does not authenticate a particular seeded implementation's random law.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import dense_cover

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SCHEMA = 'packed-canonical-gl32-whole-replay-1'
SPARSE_SCOPE = dict(K=1 << 20, N=1 << 21, group_count=2048, outer='BCH256128',
    block_rows=4, block_columns=8, support_min=5,
    mixing='independent uniform GL32 per group/block; no additional GF16 stage',
    routing='independent shared column shuffle per group and independent regional shuffles',
    inner=dict(t=128, s=19, updates=2))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_source(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    record = json.loads(raw)
    require(isinstance(record, dict), 'JSON object receipt required')
    return record, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())


def claim(record):
    require(type(record.get('distance')) in (str, int), 'exact rational distance required')
    distance = Q(record['distance'])
    require(0 < distance < Q(1, 2) and type(record.get('threshold')) is int
            and record['threshold'] == int(distance*(1 << 21)), 'inconsistent distance and cutoff')
    return distance, record['threshold']


def partition_cells(root, leaves):
    """Reconstruct a complete binary interval partition using exact rationals."""
    require(isinstance(root, list) and len(root) == 2, 'two rational root endpoints required')
    lo, hi = map(Q, root)
    require(0 <= lo < hi <= 1 and isinstance(leaves, dict) and bool(leaves), 'nonempty dense domain required')
    require(all(isinstance(p, str) and not set(p)-{'0', '1'} for p in leaves), 'binary dense paths required')
    paths = sorted(leaves)
    require(all(not b.startswith(a) for a, b in zip(paths, paths[1:])), 'overlapping dense cells')
    require(sum((Q(1, 1 << len(p)) for p in paths), Q(0)) == 1, 'dense partition has a gap')
    cells = {}
    for path in paths:
        left, right = lo, hi
        for digit in path:
            middle = (left+right)/2
            if digit == '0':
                right = middle
            else:
                left = middle
        cells[path] = (left, right)
    return cells


def complete_dense(record):
    cover = dense_cover.validate_record(record)
    claim(record)
    require(record['minimum_groups'] == 33 and record['maximum_groups'] == 2048,
            'dense occupancy scope must be exactly33..2048')
    require(cover['unresolved'] == {}, 'unresolved dense cells remain')
    cells = partition_cells(record['root'], cover['leaves'])
    require(all(isinstance(row, dict) and isinstance(row.get('witness'), dict)
                for row in cover['leaves'].values()), 'one witness per dense cell required')
    return cells


def prepare(dense_path, sparse_paths):
    """Validate scope before starting either expensive replay; ignore old bounds."""
    dense, dense_source = read_source(dense_path)
    complete_dense(dense)
    distance, threshold = claim(dense)
    premises = dense.get('outer_premises')
    require(isinstance(premises, dict) and bool(premises), 'dense outer premises required')
    sources, seen, paths, count_hashes = [], set(), set(), set()
    for path in sparse_paths:
        record, source = read_source(path)
        require(source['path'] not in paths, 'duplicate sparse source file')
        paths.add(source['path'])
        require(record.get('schema') == 'packed-gl32-sparse-1'
                and all(record.get(k) == v for k, v in SPARSE_SCOPE.items()),
                'sparse construction does not match canonical GL32/R2')
        source_distance, _ = claim(record)
        require(distance <= source_distance, 'sparse replay cannot raise the source distance')
        require(record.get('count_premises') == premises, 'sparse and dense count premises differ')
        digest = record.get('count_sha256')
        require(isinstance(digest, str) and len(digest) == 64
                and not set(digest)-set('0123456789abcdef'), 'sparse count fingerprint required')
        count_hashes.add(digest)
        entries, local_seen = record.get('results'), set()
        require(isinstance(entries, list), 'sparse result list required')
        for row in entries:
            require(isinstance(row, dict) and type(row.get('occupancy')) is int
                    and 1 <= row['occupancy'] <= 32 and row['occupancy'] not in local_seen,
                    'invalid or duplicate sparse occupancy')
            q = row['occupancy']
            local_seen.add(q)
            # The existing sparse replay selects these entries. The old
            # numeric endpoint is deliberately not parsed or compared.
            if row.get('passed') is not True or row.get('upper') is None:
                continue
            require(q not in seen and isinstance(row.get('details'), dict),
                    'overlapping sparse scopes or missing witness details')
            seen.add(q)
        sources.append(source)
    require(seen == set(range(1, 33)), 'passed sparse occupancies must cover exactly1..32')
    require(len(count_hashes) == 1, 'sparse count fingerprints differ')
    return dict(dense=dense, dense_source=dense_source, sparse_sources=sources,
        distance=distance, threshold=threshold, premises=premises, sparse_count_sha256=count_hashes.pop())


def dyadic(pair):
    require(isinstance(pair, list) and len(pair) == 2
            and all(type(v) is int for v in pair) and pair[0] > 0,
            'positive exact dyadic endpoint required')
    return Q(pair[0])*Q(2)**pair[1]


def compact_endpoint(value, bits):
    """Round an exact positive dyadic sum upward once for compact storage."""
    value = Q(value)
    require(value > 0 and value.denominator & (value.denominator-1) == 0,
            'positive dyadic sum required')
    numerator, exponent = value.numerator, 1-value.denominator.bit_length()
    shift = max(0, numerator.bit_length()-bits)
    if shift:
        numerator = (numerator+(1 << shift)-1) >> shift
        exponent += shift
    shift = (numerator & -numerator).bit_length()-1
    return [numerator >> shift, exponent+shift]


def _assemble_fresh(prepared, sparse, dense, precision, target_bits):
    """Internal post-replay check, not an interface for trusting saved receipts."""
    require(type(precision) is int and precision >= 128 and type(target_bits) is int
            and target_bits >= 20, 'precision>=128 and whole-code target>=20 required')
    expected_claim = prepared['distance'], prepared['threshold']
    require(claim(sparse) == claim(dense) == expected_claim, 'fresh replay cutoffs differ')
    require(sparse.get('precision') == dense.get('precision') == precision, 'fresh precision mismatch')
    require(sparse.get('schema') == 'packed-gl32-sparse-replay-1'
            and sparse.get('scope') == SPARSE_SCOPE, 'fresh sparse construction mismatch')
    require(sparse.get('sources') == prepared['sparse_sources']
            and dense.get('source') == prepared['dense_source'], 'fresh source hashes differ')
    require(sparse.get('count_premises') == dense.get('outer_premises') == prepared['premises'],
            'fresh count premises differ')
    require(sparse.get('count_sha256') == prepared['sparse_count_sha256']
            and all(dense.get(k) == prepared['dense'].get(k) for k in
                    ('expected_cdf_sha256', 'comparison_caps_sha256', 'comparison')),
            'fresh comparison fingerprints differ from the authenticated sources')
    require(sparse.get('requested') == list(range(1, 33)) and sparse.get('uncovered_sparse') == []
            and sparse.get('complete_requested') is True and sparse.get('complete_sparse_prefix') is True,
            'fresh sparse prefix is incomplete')
    rows = sparse.get('results')
    require(isinstance(rows, list) and len(rows) == 32
            and all(isinstance(row, dict) and type(row.get('occupancy')) is int for row in rows)
            and sorted(row['occupancy'] for row in rows) == list(range(1, 33)),
            'fresh sparse rows must cover each occupancy once')
    sparse_sum = sum((dyadic(row['upper']) for row in rows), Q(0))
    require(sparse_sum <= dyadic(sparse['aggregate_upper']), 'sparse aggregate omits row mass')

    cells = complete_dense(dense)
    require(dense['cover'] == prepared['dense']['cover'] and dense['root'] == prepared['dense']['root'],
            'fresh dense witness partition differs from its source')
    replay = dense.get('fresh_replay')
    require(isinstance(replay, dict) and replay.get('complete_dense') is True
            and replay.get('unresolved') == 0, 'fresh dense replay is incomplete')
    checked = replay.get('checked')
    require(isinstance(checked, list) and len(checked) == len(cells)
            and all(isinstance(row, dict) and isinstance(row.get('path'), str) for row in checked)
            and sorted(row['path'] for row in checked) == sorted(cells),
            'fresh dense rows omit, duplicate, or add cells')
    for row in checked:
        require(isinstance(row.get('cell'), list) and tuple(map(Q, row['cell'])) == cells[row['path']],
                'fresh dense cell geometry differs')
    dense_sum = sum((dyadic(row['upper']) for row in checked), Q(0))
    require(dense_sum <= dyadic(replay['upper']), 'dense aggregate omits cell mass')
    total = sparse_sum+dense_sum
    limit = Q(1, 1 << target_bits)
    require(total < limit, 'fresh whole-code sum misses the requested margin')
    upper = compact_endpoint(total, precision+32)
    require(dyadic(upper) < limit, 'compact outward endpoint misses the requested margin')
    return dict(schema=SCHEMA, ensemble=dense_cover.ENSEMBLE, verified=True, complete=True, target_bits=target_bits,
        precision=precision, distance=str(prepared['distance']), threshold=prepared['threshold'],
        construction=SPARSE_SCOPE, sparse_occupancies=[1, 32], dense_occupancies=[33, 2048],
        count_premises=prepared['premises'], sparse_count_sha256=sparse['count_sha256'],
        dense_expected_cdf_sha256=dense['expected_cdf_sha256'],
        sparse_upper=compact_endpoint(sparse_sum, precision+32),
        dense_upper=compact_endpoint(dense_sum, precision+32), upper=upper,
        aggregation='Exact sum of all fresh row dyadics, followed by one upward rounding.',
        probability_space='Independent ideal uniform GL32 maps per group/block, independent shared column shuffles per group, independent regional shuffles, and the declared R2 inner setup.',
        claim='Expected number of nonzero messages with output weight at most threshold is below 2^-target_bits; the same bound holds for distance failure probability.')


def run(dense_path, sparse_paths, output_dir, precision=384, target_bits=20, dense_workers=1):
    """Run both authenticating replays serially; never accept preexisting outputs."""
    require(type(precision) is int and precision >= 128 and type(target_bits) is int
            and target_bits >= 20, 'precision>=128 and whole-code target>=20 required')
    require(type(dense_workers) is int and 1 <= dense_workers <= 4,
            'dense worker count must be an integer in 1..4')
    output_dir = Path(output_dir).resolve()
    require(not output_dir.exists(), 'a new output directory is required')
    prepared = prepare(dense_path, sparse_paths)
    output_dir.mkdir(parents=True, exist_ok=False)
    sparse_path, dense_path = output_dir/'sparse.json', output_dir/'dense.json'
    common = ['--precision', str(precision), '--target-bits', str(target_bits)]
    commands = [
        [sys.executable, str(HERE/'sparse.py'), '--replay',
         *(v['path'] for v in prepared['sparse_sources']), '--distance', str(prepared['distance']),
         *common, '--output', str(sparse_path)],
        [sys.executable, str(HERE/'dense_cover.py'), '--replay', prepared['dense_source']['path'],
         *common, '--max-cells', '0', '--output', str(dense_path)],
    ]
    if dense_workers > 1:
        commands[1].extend(['--replay-workers', str(dense_workers)])
    for command in commands:
        print('WHOLE REPLAY starting', Path(command[1]).name, flush=True)
        subprocess.run(command, cwd=REPO, check=True)
    for source in [prepared['dense_source'], *prepared['sparse_sources']]:
        require(hashlib.sha256(Path(source['path']).read_bytes()).hexdigest() == source['sha256'],
                'an input witness receipt changed during replay')
    sparse, sparse_receipt = read_source(sparse_path)
    dense, dense_receipt = read_source(dense_path)
    result = _assemble_fresh(prepared, sparse, dense, precision, target_bits)
    result.update(witness_sources=dict(dense=prepared['dense_source'], sparse=prepared['sparse_sources']),
        fresh_receipts=dict(sparse=sparse_receipt, dense=dense_receipt), replay_commands=commands,
        dense_replay_workers=dense_workers)
    with (output_dir/'whole.json').open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('VERIFIED complete canonical GL32/R2 whole-code bound <2^-'+str(target_bits), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dense', type=Path, required=True)
    parser.add_argument('--sparse', type=Path, nargs='+', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=384)
    parser.add_argument('--target-bits', type=int, default=20)
    parser.add_argument('--dense-workers', type=int, default=1)
    args = parser.parse_args()
    run(args.dense, args.sparse, args.output_dir, args.precision, args.target_bits, args.dense_workers)


if __name__ == '__main__':
    main()
