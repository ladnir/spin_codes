"""Fresh whole-code replay for one explicit canonical GL32 hill variant.

The complete dense and sparse searches must describe the same R2/R3/R4
construction and exactly the same claim. Their old bounds are not inputs.
Fresh replay programs run serially into a new directory. Only their new
per-cell and per-occupancy dyadic endpoints enter the exact final sum.

Success bounds the expected number of nonzero messages of weight at most
the cutoff, and hence the setup failure probability, under the ideal setup
law. It does not certify an implementation's seed generator or running time.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import hill_cover
import sparse_hill_variant

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SCHEMA = 'packed-canonical-gl32-whole-hill-replay-1'
SPARSE_KEYS = ('K', 'N', 'outer', 'group_count', 'block_rows', 'block_columns',
               'mixing', 'routing', 'inner', 'support_min')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_source(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    record = json.loads(raw)
    require(isinstance(record, dict), 'JSON object required')
    return record, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())


def proof_source_manifest():
    """Record source bytes for replay provenance, not as proof premises.

    Rescanning also detects added or removed Python sources. The production
    headers fix the BCH matrix and selected inner maps used by authentication.
    """
    python_root = HERE.parent.resolve()
    require(python_root.is_dir(), 'proof source directory is missing')
    python_paths = sorted(python_root.rglob('*.py'), key=lambda path: str(path))
    require(python_paths, 'proof source directory has no Python sources')
    headers = [REPO/'spin/src/kernels/generated'/name
               for name in ('BchCircuit.h', 'SelectedMaps.h')]
    require(all(path.is_file() for path in headers), 'production proof header is missing')
    files = [dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
             for path in [*python_paths, *headers]]
    require(len({row['path'] for row in files}) == len(files), 'duplicate proof source path')
    return dict(schema='packed-hill-proof-source-provenance-1', python_root=str(python_root),
        files=files, purpose='Source-byte provenance and unchanged-source guard only; '
                            'saved hashes are not mathematical premises or proof authority.')


def check_proof_sources(manifest):
    require(proof_source_manifest() == manifest, 'proof implementation sources changed during replay')


def claim(record):
    require(type(record.get('distance')) in (str, int), 'exact rational distance required')
    distance = Q(record['distance'])
    require(0 < distance < Q(1, 2) and type(record.get('threshold')) is int
        and record['threshold'] == int(distance*(1 << 21)), 'matching distance and cutoff required')
    return distance, record['threshold']


def count_fingerprints(premises):
    """Check both serialization conventions on the same exact transported CDF.

    This is a consistency check, not authentication of saved count premises.
    Each child later regenerates and authenticates those premises itself.
    """
    from canonical_counts import transport_cdf
    from local_models import full_block
    caps = premises.get('canonical_cdf')
    require(isinstance(caps, list) and len(caps) == 33 and caps[-1] == (1 << 512)-1,
        'complete four-row canonical cumulative caps required')
    counts = transport_cdf(caps, full_block(8))
    dense_hash = hill_cover.cell_search.dense.fingerprint(counts)
    sparse_hash = hashlib.sha256(json.dumps([str(v) for v in counts], separators=(',', ':')).encode()).hexdigest()
    floor = next(j for j, value in enumerate(counts) if value)
    return dense_hash, sparse_hash, floor


def prepare(dense_path, sparse_path):
    dense, dense_source = read_source(dense_path)
    sparse, sparse_source = read_source(sparse_path)
    hill_cover.validate_record(dense, complete=True)
    scope = dense['scope']
    require(scope['minimum_groups'] == 33 and scope['maximum_groups'] == 2048,
        'dense occupancies must be exactly33..2048')
    selected = sparse_hill_variant.validate_record(sparse)
    require(sparse.get('requested') == list(range(1, 33))
        and sparse.get('complete_requested') is True and sparse.get('complete_sparse_prefix') is True
        and sparse.get('uncovered_sparse') == []
        and sorted(row['occupancy'] for row in selected) == list(range(1, 33))
        and len(sparse['results']) == 32, 'complete successful sparse occupancies1..32 required')
    require(claim(scope) == claim(sparse), 'dense and sparse claims must match exactly')
    require(sparse['inner']['updates'] == scope['updates'], 'actual inner variants differ')
    require(sparse['count_premises'] == scope['outer_premises'], 'outer count premises differ')
    refinement = 'ekr' if hill_cover.validate_scope(scope) else 'incidence'
    require(sparse['count_refinement'] == refinement, 'count refinement selectors differ')
    dense_hash, sparse_hash, floor = count_fingerprints(scope['outer_premises'])
    require(scope['expected_cdf_sha256'] == dense_hash and sparse['count_sha256'] == sparse_hash
        and sparse['support_min'] == floor, 'count fingerprints or support floor differ from common premises')
    return dict(dense=dense, sparse=sparse, dense_source=dense_source, sparse_source=sparse_source,
        scope=scope, construction={key: sparse[key] for key in SPARSE_KEYS}, refinement=refinement,
        distance=claim(scope)[0], threshold=scope['threshold'], count_premises=scope['outer_premises'],
        dense_count_sha256=dense_hash, sparse_count_sha256=sparse_hash)


def dyadic(pair):
    require(isinstance(pair, list) and len(pair) == 2
        and all(type(value) is int for value in pair) and pair[0] > 0, 'positive dyadic endpoint required')
    return Q(pair[0])*Q(2)**pair[1]


def compact_endpoint(value, bits):
    value = Q(value)
    require(type(bits) is int and bits >= 1 and value > 0
        and value.denominator & (value.denominator-1) == 0, 'positive dyadic and integer precision required')
    numerator, exponent = value.numerator, 1-value.denominator.bit_length()
    shift = max(0, numerator.bit_length()-bits)
    if shift:
        numerator = (numerator+(1 << shift)-1) >> shift
        exponent += shift
    shift = (numerator & -numerator).bit_length()-1
    return [numerator >> shift, exponent+shift]


def _assemble_fresh(prepared, sparse, dense, precision, target_bits):
    """Internal post-child validation; never a saved-receipt acceptance mode."""
    require(type(precision) is int and precision >= 128
        and type(target_bits) is int and target_bits >= 40, 'precision>=128 and whole target>=40 required')
    require(sparse.get('schema') == 'packed-gl32-sparse-hill-variant-replay-1'
        and sparse.get('precision') == precision and sparse.get('target_bits') == target_bits+1,
        'fresh sparse replay schema, precision, and target required')
    require(sparse.get('scope') == prepared['construction']
        and sparse.get('sources') == [prepared['sparse_source']]
        and claim(sparse) == (prepared['distance'], prepared['threshold'])
        and sparse.get('count_premises') == prepared['count_premises']
        and sparse.get('count_sha256') == prepared['sparse_count_sha256']
        and sparse.get('count_refinement') == prepared['refinement']
        and sparse.get('joint_return_through') == 2 and sparse.get('lazy_density_through') == 4,
        'fresh sparse construction, count premises, or provenance mismatch')
    require(sparse.get('requested') == list(range(1, 33)) and sparse.get('uncovered_sparse') == []
        and sparse.get('complete_requested') is True and sparse.get('complete_sparse_prefix') is True,
        'fresh sparse coverage is incomplete')
    rows = sparse.get('results')
    require(isinstance(rows, list) and len(rows) == 32
        and all(isinstance(row, dict) and type(row.get('occupancy')) is int
                and row.get('passed') is True for row in rows)
        and sorted(row['occupancy'] for row in rows) == list(range(1, 33)),
        'fresh sparse rows must cover each occupancy exactly once')
    sparse_sum = sum((dyadic(row['upper']) for row in rows), Q(0))
    require(all(dyadic(row['upper']) < Q(2)**-(target_bits+1) for row in rows)
        and sparse_sum <= dyadic(sparse['aggregate_upper']), 'fresh sparse row/aggregate inconsistency')
    cells = hill_cover.validate_record(dense, complete=True)
    require(dense.get('precision') == precision and dense.get('scope') == prepared['scope']
        and dense.get('source') == prepared['dense_source']
        and dense.get('cover') == prepared['dense']['cover'], 'fresh dense source or construction mismatch')
    replay = dense.get('fresh_replay')
    require(isinstance(replay, dict) and replay.get('precision') == precision
        and replay.get('target_bits') == target_bits+1 and replay.get('complete_dense') is True
        and replay.get('passed') is True and replay.get('whole_code_certificate') is False
        and replay.get('scope') == prepared['scope'], 'fresh dense replay is incomplete or has wrong scope')
    checked = replay.get('checked')
    require(isinstance(checked, list) and len(checked) == len(cells)
        and all(isinstance(row, dict) and isinstance(row.get('path'), str) for row in checked)
        and sorted(row['path'] for row in checked) == sorted(cells),
        'fresh dense rows omit, duplicate, or add cells')
    for row in checked:
        require(tuple(map(Q, row['cell'])) == cells[row['path']], 'fresh dense interval geometry differs')
    dense_sum = sum((dyadic(row['upper']) for row in checked), Q(0))
    require(dense_sum <= dyadic(replay['aggregate_upper']) < Q(2)**-(target_bits+1),
        'fresh dense aggregate does not match its successful replay')
    total = sparse_sum+dense_sum
    upper = compact_endpoint(total, precision+32)
    require(total < Q(2)**-target_bits and dyadic(upper) < Q(2)**-target_bits,
        'fresh whole-code sum misses the requested strict margin')
    return dict(schema=SCHEMA, verified=True, complete=True, target_bits=target_bits,
        precision=precision, ensemble=prepared['scope']['ensemble'], construction=prepared['construction'],
        distance=str(prepared['distance']), threshold=prepared['threshold'],
        sparse_occupancies=[1, 32], dense_occupancies=[33, 2048],
        count_premises=prepared['count_premises'], sparse_count_sha256=prepared['sparse_count_sha256'],
        dense_expected_cdf_sha256=prepared['dense_count_sha256'],
        dense_comparison_caps_sha256=prepared['scope']['comparison_caps_sha256'],
        sparse_upper=compact_endpoint(sparse_sum, precision+32),
        dense_upper=compact_endpoint(dense_sum, precision+32), upper=upper,
        aggregation='Exact sum of every fresh per-occupancy and per-cell dyadic, then one upward rounding.',
        probability_space='Independent ideal uniform GL32 per group/block; independent shared column shuffles per group; independent regional shuffles; the declared actual inner update count.',
        claim='Expected number of nonzero messages of output weight at most threshold is below 2^-target_bits; this also bounds distance failure probability.')


def run(dense_path, sparse_path, output_dir, precision=384, target_bits=40, dense_workers=1):
    require(type(precision) is int and precision >= 128 and type(target_bits) is int
        and target_bits >= 40, 'precision>=128 and whole target>=40 required')
    require(type(dense_workers) is int and 1 <= dense_workers <= 4,
        'dense worker count must be an integer in 1..4')
    output_dir = Path(output_dir).resolve()
    require(not output_dir.exists(), 'new output directory required')
    proof_sources = proof_source_manifest()
    prepared = prepare(dense_path, sparse_path)
    output_dir.mkdir(parents=True, exist_ok=False)
    sparse_output, dense_output = output_dir/'sparse.json', output_dir/'dense.json'
    common = ['--precision', str(precision), '--target-bits', str(target_bits+1)]
    commands = [
        [sys.executable, str(HERE/'sparse_hill_variant.py'), '--replay', prepared['sparse_source']['path'],
         '--updates', str(prepared['scope']['updates']), '--refinement', prepared['refinement'],
         '--distance', str(prepared['distance']), *common, '--output', str(sparse_output)],
        [sys.executable, str(HERE/'hill_cover.py'), prepared['dense_source']['path'], '--replay',
         *common, '--output', str(dense_output)],
    ]
    if dense_workers > 1:
        commands[1].extend(['--replay-workers', str(dense_workers)])
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       PYTHONDONTWRITEBYTECODE='1')
    for command, output in zip(commands, (sparse_output, dense_output)):
        check_proof_sources(proof_sources)
        require(not output.exists(), 'replay output unexpectedly exists')
        print('WHOLE HILL REPLAY starting', Path(command[1]).name, flush=True)
        subprocess.run(command, cwd=REPO, env=environment, check=True)
        require(output.is_file(), 'fresh child did not write its receipt')
    for source in (prepared['dense_source'], prepared['sparse_source']):
        require(hashlib.sha256(Path(source['path']).read_bytes()).hexdigest() == source['sha256'],
            'source witness receipt changed during replay')
    sparse, sparse_metadata = read_source(sparse_output)
    dense, dense_metadata = read_source(dense_output)
    require(dense.get('fresh_replay', {}).get('workers', 1) == dense_workers,
        'fresh dense worker configuration differs from the requested replay')
    result = _assemble_fresh(prepared, sparse, dense, precision, target_bits)
    result.update(witness_sources=dict(dense=prepared['dense_source'], sparse=prepared['sparse_source']),
        fresh_receipts=dict(dense=dense_metadata, sparse=sparse_metadata), replay_commands=commands,
        proof_sources=proof_sources, dense_replay_workers=dense_workers)
    check_proof_sources(proof_sources)
    with (output_dir/'whole.json').open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('VERIFIED fresh whole-code', prepared['scope']['ensemble'], 'bound <2^-'+str(target_bits), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dense', type=Path, required=True)
    parser.add_argument('--sparse', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=384)
    parser.add_argument('--target-bits', type=int, default=40)
    parser.add_argument('--dense-workers', type=int, default=1)
    args = parser.parse_args()
    run(args.dense, args.sparse, args.output_dir, args.precision, args.target_bits, args.dense_workers)


if __name__ == '__main__':
    main()
