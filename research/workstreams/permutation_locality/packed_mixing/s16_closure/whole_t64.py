"""Fresh whole-code replay for the selected t64/S16 uniform-GL construction.

The inputs are search witnesses, not certificates or numerical premises.
This program requires a complete dense partition and sparse occupancies
1..32. It runs both replays in new processes and a new output directory.
Only their newly computed per-cell and per-occupancy dyadics enter the sum.

The claim concerns independent ideal uniform setup randomness. Each physical
64-symbol step emits y=x+Aq and updates q'=Mq+Cx with uniform M in GL16.
State starts at zero, persists between steps, and is not flushed. Two such
steps form each 128-symbol proof macro. This does not certify a seeded
implementation's randomness, binary, or running time.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import dense_t64
import sparse_t64

HERE = Path(__file__).resolve().parent
PERMUTATION = HERE.parents[1]
REPO = HERE.parents[4]
SCHEMA = 'packed-gl32-t64-s16-whole-replay-1'
DENSE_SCOPE = dict(schema='packed-gl32-t64-s16-dense-context-1',
    K=1 << 20, N=1 << 21, physical_t=64, state_bits=16,
    physical_steps=32768, macro_t=128, macro_windows=32, macro_steps=16384,
    physical_steps_per_macro=2, regions=256, macros_per_region=64,
    outer='fixed BCH[256,128]', block_rows=4, block_width=8, groups=2048,
    mixing='independent uniform GL32 on each canonical four-row/eight-column block',
    routing='shared uniform column shuffle per four-row group; independent uniform packet shuffle per region',
    inner_recurrence='y=x+Aq; next_q=Mq+Cx; independent uniform GL16 each physical step; zero initial state; no flush',
    birth_density='capped', minimum_groups=33, maximum_groups=2048,
    comparison='direct-expected-shell-majorant', outer_authenticated=True)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON object key: '+key)
        result[key] = value
    return result


def read_source(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    record = json.loads(raw, object_pairs_hook=_object)
    require(isinstance(record, dict), 'JSON object required')
    return record, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest())


def source_metadata(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def unchanged(source):
    require(source_metadata(source['path']) == source, 'source path or bytes changed')


def proof_source_manifest(inputs):
    """Byte provenance only; both children independently rebuild premises."""
    paths = sorted(PERMUTATION.rglob('*.py'), key=str)
    require(paths, 'proof Python sources are missing')
    paths += [REPO/'spin/src/kernels/generated'/name for name in ('BchCircuit.h', 'SelectedMaps.h')]
    paths += [dense_t64.kernel_t64.SELECTED_MAP, *[Path(row['path']) for row in inputs]]
    unique = sorted({path.resolve() for path in paths}, key=str)
    require(all(path.is_file() for path in unique), 'proof source or input is missing')
    return dict(schema='t64-proof-source-manifest-1', python_root=str(PERMUTATION.resolve()),
        files=[source_metadata(path) for path in unique],
        purpose='Source-byte provenance and unchanged-source guard; hashes are not mathematical premises.')


def claim(record):
    require(type(record.get('distance')) in (str, int), 'exact distance required')
    distance = Q(record['distance'])
    require(0 < distance < Q(1, 2) and type(record.get('threshold')) is int
        and record['threshold'] == int(distance*(1 << 21)), 'matching exact distance and cutoff required')
    return distance, record['threshold']


def count_fingerprints(premises):
    """Cross-check saved metadata on common exact CDF and shell transforms.

    This does not authenticate the saved premises. Both child replays call
    authenticated_bch_cdf and reject a mismatch against fresh BCH premises.
    """
    from canonical_counts import transport_cdf
    from local_models import full_block
    from monotone import transport_shells
    caps = premises.get('canonical_cdf')
    require(premises.get('schema') == sparse_t64.SCOPE['outer_count_refinement']
        and premises.get('block_width') == 8 and premises.get('block_rows') == 4
        and premises.get('groups_per_outer_word') == 32
        and isinstance(caps, list) and len(caps) == 33
        and all(type(v) is int for v in caps) and caps[-1] == (1 << 512)-1,
        'canonical four-row width-eight EKR count premises required')
    counts = transport_cdf(caps, full_block(8))
    shells = transport_shells(caps, full_block(8))
    fingerprint = dense_t64.prior.cell_search.dense.fingerprint
    sparse_hash = hashlib.sha256(json.dumps([str(v) for v in counts], separators=(',', ':')).encode()).hexdigest()
    return fingerprint(counts), fingerprint(shells), sparse_hash, next(j for j, value in enumerate(counts) if value)


def fresh_map_record():
    # Local map preparation is bounded. No outer numerical cover is run here.
    return dense_t64.kernel_t64.prepare(birth_density='capped')[1]


def prepare(dense_path, sparse_path):
    dense, dense_source = read_source(dense_path)
    sparse, sparse_source = read_source(sparse_path)
    cells = dense_t64.validate_cover_record(dense, complete=True)
    scope = dense['scope']
    require(all(scope.get(key) == value for key, value in DENSE_SCOPE.items()),
        'the complete explicit physical64/macro128 dense scope is required')
    require(type(scope.get('variance_bins')) is int and 1 <= scope['variance_bins'] <= 64
        and Q(scope['base_tilt']) == Q(3, 16), 'supported dense proof parameters required')
    selected = sparse_t64.validate_record(sparse)
    require(sorted(sparse['requested']) == list(range(1, 33))
        and len(sparse['results']) == 32
        and sorted(row['occupancy'] for row in selected) == list(range(1, 33)),
        'every sparse occupancy1..32 must have one successful full-support witness')
    require(claim(scope) == claim(sparse), 'dense and sparse claims differ')
    require(sparse['map_record'] == scope['maps'] == fresh_map_record()
        and scope.get('maps_sha256') == dense_t64.prior.fingerprint(scope['maps']),
        'search maps or map source differ from fresh selected-map preparation')
    require(sparse['count_premises'] == scope['outer_premises'], 'outer count premises differ')
    cdf_hash, shell_hash, sparse_hash, floor = count_fingerprints(scope['outer_premises'])
    require(scope.get('expected_cdf_sha256') == cdf_hash
        and scope.get('comparison_caps_sha256') == shell_hash
        and sparse.get('count_sha256') == sparse_hash and sparse.get('support_min') == floor,
        'count transforms, fingerprints, or support floors differ')
    proposal, proposal_source = read_source(dense['source']['path'])
    require(dense['source'] == proposal_source, 'dense mixture proposal source changed')
    old = proposal['scope']
    for key in ('outer_premises', 'expected_cdf_sha256', 'comparison_caps_sha256', 'mixture'):
        require(old.get(key) == scope.get(key), 'dense mixture proposal metadata differ: '+key)
    # Reconstruct the complete mean domain from the actual mixture, not flags.
    mixture = scope['mixture']
    require(isinstance(mixture, list) and mixture
        and all(Q(row['mass']) > 0 and 0 < Q(row['activity']) <= 1 for row in mixture),
        'positive exact mixture components required')
    tilt = Q(scope['base_tilt'])
    minimum = Q(33, 2048)*min(Q(row['activity'])*tilt/(1-Q(row['activity'])+Q(row['activity'])*tilt)
                            for row in mixture)
    require(tuple(map(Q, scope['root'])) == (minimum, Q(1)), 'dense cover root omits part of the full mean domain')
    inputs = [dense_source, sparse_source, proposal_source, scope['maps']['source']]
    for source in inputs:
        unchanged(source)
    return dict(dense=dense, sparse=sparse, dense_source=dense_source, sparse_source=sparse_source,
        proposal_source=proposal_source, inputs=inputs, cells=cells, scope=scope,
        distance=claim(scope)[0], threshold=scope['threshold'], count_premises=scope['outer_premises'],
        dense_count_sha256=cdf_hash, comparison_sha256=shell_hash, sparse_count_sha256=sparse_hash)


def dyadic(pair):
    require(isinstance(pair, list) and len(pair) == 2
        and all(type(v) is int for v in pair) and pair[0] > 0, 'positive exact dyadic endpoint required')
    return Q(pair[0])*Q(2)**pair[1]


def exact_endpoint(value):
    value = Q(value)
    require(value > 0 and value.denominator & (value.denominator-1) == 0, 'positive dyadic sum required')
    return [value.numerator, 1-value.denominator.bit_length()]


def compact_endpoint(value, bits):
    """Round a positive dyadic upward to at most ``bits`` mantissa bits.

    Tiny cells can make an exact sum's denominator enormous. The exact sum
    remains the acceptance input; only its stored upper bound is compacted.
    """
    value = Q(value)
    require(type(bits) is int and bits >= 1 and value > 0
        and value.denominator & (value.denominator-1) == 0,
        'positive dyadic and integer mantissa precision required')
    numerator, exponent = value.numerator, 1-value.denominator.bit_length()
    shift = max(0, numerator.bit_length()-bits)
    if shift:
        numerator = (numerator+(1 << shift)-1) >> shift
        exponent += shift
    # A rounding carry can produce 2**bits, whose normalized mantissa is 1.
    shift = (numerator & -numerator).bit_length()-1
    return [numerator >> shift, exponent+shift]


def _assemble_fresh(prepared, sparse, dense, precision, target_bits):
    """Internal post-process check, never an accept-saved-receipts command."""
    scope = prepared['scope']
    require(sparse.get('schema') == 't64-s16-uniform-gl-sparse-replay-1'
        and sparse.get('precision') == precision and sparse.get('target_bits') == target_bits+1
        and sparse.get('K') == 1 << 20 and sparse.get('N') == 1 << 21
        and sparse.get('inner') == sparse_t64.INNER
        and sparse.get('kernel_variant') == sparse_t64.KERNEL_VARIANT
        and all(sparse.get(key) == value for key, value in sparse_t64.SCOPE.items())
        and sparse.get('sources') == [prepared['sparse_source']]
        and sparse.get('map_record') == scope['maps']
        and sparse.get('count_premises') == prepared['count_premises']
        and sparse.get('count_sha256') == prepared['sparse_count_sha256']
        and sparse.get('support_min') == prepared['sparse']['support_min']
        and claim(sparse) == (prepared['distance'], prepared['threshold']),
        'fresh sparse scope, count premises, or provenance mismatch')
    rows = sparse.get('results')
    require(sparse.get('requested') == list(range(1, 33))
        and sparse.get('complete_requested') is True and sparse.get('complete_sparse_prefix') is True
        and sparse.get('whole_code_certificate') is False
        and isinstance(rows, list) and len(rows) == 32
        and all(isinstance(row, dict) and type(row.get('occupancy')) is int and row.get('passed') is True for row in rows)
        and sorted(row['occupancy'] for row in rows) == list(range(1, 33)),
        'fresh sparse replay omits, duplicates, or fails an occupancy')
    sparse_sum = sum((dyadic(row['upper']) for row in rows), Q(0))
    require(all(dyadic(row['upper']) < Q(2)**-(target_bits+1) for row in rows)
        and sparse_sum <= dyadic(sparse['aggregate_upper']), 'fresh sparse row/aggregate inconsistency')
    require(dense.get('schema') == 'packed-gl32-t64-s16-dense-replay-1'
        and dense.get('precision') == precision and dense.get('scope') == scope
        and dense.get('source') == prepared['proposal_source']
        and dense.get('input_source') == prepared['dense_source']
        and dense.get('complete_dense') is True and dense.get('whole_code_certificate') is False,
        'fresh dense replay scope, completeness, or provenance mismatch')
    cells, checked = prepared['cells'], dense.get('cells')
    require(isinstance(checked, list) and len(checked) == len(cells)
        and all(isinstance(row, dict) and isinstance(row.get('path'), str) for row in checked)
        and sorted(row['path'] for row in checked) == sorted(cells),
        'fresh dense replay omits, duplicates, or adds a cell')
    for row in checked:
        require(tuple(map(Q, row['cell'])) == cells[row['path']], 'fresh dense interval differs from its exact path')
    dense_sum = sum((dyadic(row['upper']) for row in checked), Q(0))
    require(dense_sum <= dyadic(dense['aggregate_upper']), 'fresh dense row/aggregate inconsistency')
    total = dense_sum+sparse_sum
    require(total < Q(2)**-target_bits, 'fresh exact whole-code sum misses the strict target')
    sparse_upper = compact_endpoint(sparse_sum, precision+32)
    dense_upper = compact_endpoint(dense_sum, precision+32)
    total_upper = compact_endpoint(total, precision+32)
    require(dyadic(total_upper) < Q(2)**-target_bits,
        'upward-rounded whole-code endpoint misses the strict target')
    return dict(schema=SCHEMA, verified=True, complete=True, whole_code_certificate=True,
        precision=precision, target_bits=target_bits, distance=str(prepared['distance']),
        threshold=prepared['threshold'], K=1 << 20, N=1 << 21, rate='1/2',
        inner=dict(sparse_t64.INNER), maps=scope['maps'], dense_scope=scope,
        sparse_occupancies=[1, 32], dense_occupancies=[33, 2048],
        count_premises=prepared['count_premises'], sparse_count_sha256=prepared['sparse_count_sha256'],
        dense_expected_cdf_sha256=prepared['dense_count_sha256'], comparison_sha256=prepared['comparison_sha256'],
        sparse_upper=sparse_upper, dense_upper=dense_upper, upper=total_upper,
        aggregation='Exact Fraction sum of all fresh per-occupancy and per-cell dyadic endpoints; '
            'each reported aggregate is then rounded upward once to precision+32 mantissa bits. '
            'Both exact and rounded whole sums must pass the strict target.',
        probability_space='Independent ideal uniform GL32 per canonical four-row/eight-column block; '
            'independent shared column shuffles per group and independent regional packet shuffles; '
            'fresh independent uniform GL16 at every physical 64-symbol step.',
        claim='The expected number of nonzero messages with output weight at most threshold is below '
            '2^-target_bits. Therefore the probability over ideal setup that such a message exists is below this bound.',
        implementation_claim=False)


def run(dense_path, sparse_path, output_dir, precision=384, target_bits=40, dense_workers=1):
    require(type(precision) is int and precision >= 128 and type(target_bits) is int
        and target_bits >= 40, 'precision>=128 and whole target>=40 required')
    require(type(dense_workers) is int and 1 <= dense_workers <= 4, 'dense workers must be an integer1..4')
    output_dir = Path(output_dir).resolve()
    require(not output_dir.exists(), 'new output directory required')
    # Capture implementation bytes before even the bounded map preparation.
    _, dense_meta = read_source(dense_path)
    _, sparse_meta = read_source(sparse_path)
    initial = proof_source_manifest([dense_meta, sparse_meta])
    prepared = prepare(dense_path, sparse_path)
    require(proof_source_manifest([dense_meta, sparse_meta]) == initial, 'sources changed during preflight')
    manifest = proof_source_manifest(prepared['inputs'])
    output_dir.mkdir(parents=True, exist_ok=False)
    sparse_output, dense_output = output_dir/'sparse.json', output_dir/'dense.json'
    commands = [
        [sys.executable, '-B', str(HERE/'sparse_t64.py'), '--replay', prepared['sparse_source']['path'],
         '--precision', str(precision), '--target-bits', str(target_bits+1),
         '--distance', str(prepared['distance']), '--output', str(sparse_output)],
        [sys.executable, '-B', str(HERE/'dense_t64.py'), 'replay', '--input', prepared['dense_source']['path'],
         '--source', prepared['proposal_source']['path'], '--precision', str(precision),
         '--variance-bins', str(prepared['scope']['variance_bins']), '--birth-density', 'capped',
         '--distance', str(prepared['distance']), '--output', str(dense_output)],
    ]
    if dense_workers > 1:
        commands[1] = [sys.executable, '-B', str(HERE/'replay_t64_parallel.py'),
            '--input', prepared['dense_source']['path'], '--source', prepared['proposal_source']['path'],
            '--precision', str(precision), '--workers', str(dense_workers), '--output', str(dense_output)]
    environment = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       PYTHONDONTWRITEBYTECODE='1')
    receipts = []
    for command, output in zip(commands, (sparse_output, dense_output)):
        require(proof_source_manifest(prepared['inputs']) == manifest, 'proof sources or inputs changed before replay')
        require(not output.exists(), 'fresh child output unexpectedly exists')
        print('T64 WHOLE starting', output.stem, 'replay', flush=True)
        with output.with_suffix('.log').open('x', encoding='utf-8') as log:
            subprocess.run(command, cwd=REPO, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
        require(output.is_file() and not output.is_symlink(), 'fresh child did not write a regular receipt')
        require(proof_source_manifest(prepared['inputs']) == manifest, 'proof sources or inputs changed during replay')
        receipts.append(read_source(output))
    (sparse, sparse_meta), (dense, dense_meta) = receipts
    require(dense.get('workers', 1) == dense_workers, 'fresh dense worker count differs')
    if dense_workers > 1:
        metadata = dense.get('worker_metadata')
        require(dense.get('worker_start_method') == 'spawn' and isinstance(metadata, list)
            and 1 <= len(metadata) <= dense_workers
            and all(isinstance(row, dict) and type(row.get('pid')) is int and row['pid'] > 0 for row in metadata)
            and len({row['pid'] for row in metadata}) == len(metadata)
            and all(row.get('precision') == precision and row.get('source') == prepared['proposal_source']
                and row.get('scope_sha256') == dense_t64.prior.fingerprint(prepared['scope'])
                and row.get('initialization') == 'fresh dense_t64.fresh_model in spawned process' for row in metadata),
            'fresh spawned dense worker provenance is missing or inconsistent')
    result = _assemble_fresh(prepared, sparse, dense, precision, target_bits)
    after = proof_source_manifest(prepared['inputs'])
    require(after == manifest, 'proof sources or inputs changed after replay')
    for source in [*prepared['inputs'], sparse_meta, dense_meta]:
        unchanged(source)
    result.update(witness_sources=dict(dense=prepared['dense_source'], sparse=prepared['sparse_source'],
        mixture_proposal=prepared['proposal_source']), fresh_receipts=dict(dense=dense_meta, sparse=sparse_meta),
        replay_commands=commands, source_manifest_before=manifest, source_manifest_after=after,
        dense_replay_workers=dense_workers, dense_worker_metadata=dense.get('worker_metadata', []))
    with (output_dir/'whole.json').open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print('VERIFIED fresh t64/S16 whole-code bound <2^-'+str(target_bits), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dense', type=Path, required=True)
    parser.add_argument('--sparse', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--precision', type=int, default=384)
    parser.add_argument('--target-bits', type=int, default=40)
    parser.add_argument('--dense-workers', type=int, choices=range(1, 5), default=1)
    args = parser.parse_args()
    run(args.dense, args.sparse, args.output_dir, args.precision, args.target_bits, args.dense_workers)


if __name__ == '__main__':
    main()
