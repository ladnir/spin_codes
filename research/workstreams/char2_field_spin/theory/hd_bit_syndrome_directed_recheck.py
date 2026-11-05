"""Independent final-stage audit of a pinned selected-class directed receipt.

The trust boundary is the authenticated upward COLUMN cache, not the saved
word or aggregate. For every labelled A, this verifier imports each binary64
column entry as an exact dyadic and performs sixteen successive row-vector
products in Arb512. It independently takes the eighteenth-root sum and
applies the source count, Chernoff factor, RS density and projective factor.

This does not rebuild the regional caps or certify the column construction.
It checks only r18, homogeneous arity8, Fmax<=1 for the bound construction.
No root arithmetic, aggregation, geometry or regional modules are imported.
"""
from fractions import Fraction as F
from hashlib import sha256
from io import BytesIO
from itertools import combinations
import json
from math import comb, isfinite, nextafter
from pathlib import Path
import re
from time import perf_counter
from zipfile import ZipFile

import numpy as np
from flint import arb, ctx

DIRECTORY = Path(__file__).resolve().parent
PRECISION, ROWS, ARITY, GROUPS, COLUMNS = 512, 12870, 8, 18, 16
MARK, LIMIT, DENSITY_EXPONENT = F(993, 1000), 734003, 1153
CONSTRUCTION = 'G16-RS16x8-D2-bit-syndrome-D1-research'
DRIVER_SHA = '0ea0597506e167bc9da77c31e6e2f5679d9cbfd8e7ee65d133e56453cce00086'
SOURCE_PINS = {
    'hd_bit_syndrome_atoms_directed.py': '8dad6f3a541e91cd6b14b925745b514f0539eb07c7391a7eee4ff098c6b2ceb6',
    'hd_bit_syndrome_regions_directed.py': '431442dbc6c3772a9f818ff624720ecaab0e6208c67601f66cf23ec0613e129b',
    'hd_bit_syndrome_outer.py': 'a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d',
    'hd_bit_syndrome_kernel.py': 'a27e5557aceaf04273c3e67ce0a01a5d5746c6d01a2ecda8881205ba23541c77',
    'hd_bit_syndrome_gate.py': '3f8fa49998b1aacc00a964e5c6b83611d62510424afeca98875979aea8063505',
    'hd_bit_syndrome_orbit_gate.py': '1aad0e0c7b88397f18f5c81806ba04b4a2aab64015396e740af4155ea7636023',
    'hd_post_refresh_quiet_rows.py': '75e38819d225e9771ade1ed90509b7d4582a2637b24b6522f5b62b0200eb49e7',
    'hd_positive_upward.py': 'eaa44c7f3976220d1a088925cf51f2f77917b905353fc122e64ec1fdb0a0133c',
    'hd_bit_syndrome_columns_directed.py': 'c21eb93b5ceb7d59fb5f4188553b5007c82273c374185ce667c9e9669f757c88',
}
EXPECTED_SCOPE = dict(
    construction=CONSTRUCTION,
    inner='original one-state h4; packet64, state16, four RS[16,4] maps; D2 H16 D1 update',
    route='alternating even/odd adjacent pair matchings; independent pair permutations and fair swaps',
    scope='r18; homogeneous arity8; Fmax<=1; all labelled row sets and boundary paths',
    complete_requested_class=True, whole_code_certificate=False, directed=True,
    message_length=1048576, codeword_length=2097152, record_bytes=16,
    distance='35/100', forbidden_weight=LIMIT, mark=str(MARK), activity_marker='1',
    groups=GROUPS, arity=ARITY, columns=COLUMNS, rowset_count=ROWS,
    source_rowsets_fixed_across_columns=True, source_count_once=True,
    rs_normalizer_exponent=1152, projective_normalizations=1,
    young_weights=['1/16', '15/16'], young_before_word_power=True,
    information_cover=True, per_orbit_caps=True, orbits_used_for_counts_only=True,
    all_boundary_state_paths=True, initial_state=0, terminal_vector=[1, 1],
    original_one_state_inner=True, post_refresh_bank_used=False,
    saved_float_bounds_used=False, all_76_moments_verified=True)


def _poly_product(left, right, maximum):
    answer = {}
    for (w, a, b), x in left.items():
        for (v, c, d), y in right.items():
            key = w+v, a+c, b+d
            if key[0] <= maximum and key[1] <= 1 and key[2] <= 1:
                answer[key] = answer.get(key, 0)+x*y
    return answer


def _count_product(left, right, maximum):
    answer = []
    for i in range(2):
        for j in range(2):
            entry = {}
            for k in range(2):
                for key, value in _poly_product(left[2*i+k], right[2*k+j], maximum).items():
                    entry[key] = entry.get(key, 0)+value
            answer.append(entry)
    return answer


def independent_source_counts(size=8192, weight=18):
    """Exact trace((T_even T_odd)^(size/2)), truncated only above F=1.

    T_even[i,j]=z^j a^(ij), T_odd[i,j]=z^j b^(ij). The trace
    counts indexed cyclic bit strings, with each vertex weight once.
    All exponents are nonnegative, so discarded degrees cannot return.
    """
    if (type(size) is not int or size < 4 or size % 2 or type(weight) is not int
            or not 0 <= weight <= min(32, size)):
        raise ValueError('Even cycle and bounded exact support weight required')
    even = [{(j, i*j, 0): 1} for i in range(2) for j in range(2)]
    odd = [{(j, 0, i*j): 1} for i in range(2) for j in range(2)]
    base = _count_product(even, odd, weight)
    result = [{(0, 0, 0): 1}, {}, {}, {(0, 0, 0): 1}]
    exponent = size//2
    while exponent:
        if exponent & 1:
            result = _count_product(result, base, weight)
        exponent >>= 1
        if exponent:
            base = _count_product(base, base, weight)
    counts = {(a, b): result[0].get((weight, a, b), 0)+result[3].get((weight, a, b), 0)
              for a in range(2) for b in range(2)}
    n0 = counts[0, 0]
    n1 = counts[0, 1]+counts[1, 0]+counts[1, 1]
    independent = (size*comb(size-weight, weight)//(size-weight)
                   if weight <= size//2 else 0)
    if n0 != independent or counts[0, 1] != counts[1, 0]:
        raise ArithmeticError('Independent cycle-count identities disagree')
    total = comb(size, weight)
    return dict(F0=n0, F1=n1, F_le1=n0+n1, F_ge2=total-n0-n1, all_sources=total)


def exact_scaled(value, exponent):
    """Exact binary64 value times an integer power of two; no float scaling."""
    if type(exponent) is not int or abs(exponent) > 100000:
        raise ValueError('Bounded exact integer scale required')
    value = float(value)
    if not isfinite(value) or value < 0:
        raise ValueError('Finite nonnegative binary64 required')
    value = F(value)
    return value*(1 << exponent) if exponent >= 0 else value/(1 << -exponent)


def _arb(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def _endpoint(value, upper=True):
    value = value if value.is_exact() else (value.upper() if upper else value.lower())
    if not value.is_finite():
        raise ArithmeticError('Finite Arb endpoint required')
    q = value.fmpq()
    return F(int(q.p), int(q.q))


def row_word(matrices):
    """Successive ordered row-vector products, at the caller's precision."""
    left, right = arb(1), arb(0)
    used = 0
    for matrix in matrices:
        if len(matrix) != 4:
            raise ValueError('Flat row-major 2x2 matrix required')
        a, b, c, d = map(_arb, matrix)
        if not all(x >= 0 for x in (a, b, c, d)):
            raise ValueError('Nonnegative exact column inputs required')
        left, right = left*a+right*c, left*b+right*d
        used += 1
    if not used:
        raise ValueError('A nonempty ordered column sequence is required')
    return left+right


def assert_dominated(value, upper, label):
    """Require domination of the independent UPPER, not merely its lower."""
    if _endpoint(value) > F(upper):
        raise ArithmeticError('Independent upper is not dominated: '+label)


def minkowski(words, groups):
    if type(groups) is not int or groups < 1 or not words:
        raise ValueError('Nonempty words and positive exact group count required')
    total = arb(0)
    for word in words:
        if not word >= 0:
            raise ValueError('Nonnegative word interval required')
        total += word.root(groups)
    return total**groups


def independent_prefix(aggregate, count):
    if type(count) is not int or count < 0:
        raise ValueError('Exact nonnegative source count required')
    return aggregate*count*(_arb(MARK)**(-LIMIT))/(arb(255)**DENSITY_EXPONENT)


def _sha(value):
    if not isinstance(value, str) or not re.fullmatch('[0-9a-f]{64}', value):
        raise ValueError('Lowercase SHA256 required')
    return value


def _unique_object(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError('Duplicate JSON key: '+key)
        answer[key] = value
    return answer


def _checked_scope(receipt):
    for key, expected in EXPECTED_SCOPE.items():
        actual = receipt.get(key)
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError('Unsupported receipt scope: '+key)
    if type(receipt.get('precision_bits')) is not int or receipt['precision_bits'] not in (256, 384):
        raise ValueError('Expected a selected256/384 receipt')
    if receipt.get('source_sha256') != DRIVER_SHA or receipt.get('reviewed_sources') != SOURCE_PINS:
        raise ValueError('Unexpected directed source identity')


def _check_source_files(receipt):
    pins = dict(SOURCE_PINS, **{'hd_bit_syndrome_directed_gate.py': DRIVER_SHA})
    metadata = receipt.get('region_metadata', {})
    extra = metadata.get('reviewed_dependency_sha256')
    if not isinstance(extra, dict) or not extra:
        raise ValueError('Missing regional dependency pins')
    for name, digest in extra.items():
        if not isinstance(name, str) or Path(name).name != name or not name.endswith('.py'):
            raise ValueError('Invalid dependency path')
        _sha(digest)
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting dependency pins')
        pins[name] = digest
    for name, digest in pins.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Pinned source changed: '+name)
    return pins


def _cache_arrays(raw, receipt):
    with ZipFile(BytesIO(raw)) as archive:
        if sum(item.file_size for item in archive.infolist()) > 32*1024*1024:
            raise ValueError('Unexpectedly large cache')
    with np.load(BytesIO(raw), allow_pickle=False) as cache:
        expected = {'quiet', 'quiet_scale', 'nonquiet', 'nonquiet_scale', 'combined',
                    'combined_scale', 'words', 'words_scale', 'rowsets', 'mark', 'precision'}
        if len(cache.files) != len(expected) or set(cache.files) != expected:
            raise ValueError('Unexpected cache fields')
        rows = cache['rowsets']
        labels = np.array(list(combinations(range(16), 8)), dtype=np.int8)
        if rows.dtype.kind not in 'iu' or rows.shape != (ROWS, ARITY) or not np.array_equal(rows, labels):
            raise ValueError('Missing, reordered or duplicate labelled row sets')
        if cache['mark'].shape != () or str(cache['mark'].item()) != str(MARK):
            raise ValueError('Cache mark mismatch')
        if cache['precision'].shape != () or cache['precision'].dtype.kind not in 'iu' or int(cache['precision']) != receipt['precision_bits']:
            raise ValueError('Cache precision mismatch')
        arrays, scales = {}, {}
        for name in ('quiet', 'nonquiet', 'combined', 'words'):
            values = cache[name]
            shape = (ROWS,) if name == 'words' else (ROWS, 2, 2)
            if (values.dtype != np.dtype('float64') or values.shape != shape
                    or not np.all(np.isfinite(values)) or np.any(values < 0)):
                raise ValueError('Malformed finite positive cache: '+name)
            scale = cache[name+'_scale']
            if scale.shape != () or scale.dtype.kind not in 'iu' or abs(int(scale)) > 100000:
                raise ValueError('Malformed exact cache scale')
            arrays[name], scales[name] = values.copy(), int(scale)
    return arrays, scales


def load_inputs(receipt_path, expected_sha256, cache_path=None):
    path = Path(receipt_path)
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != _sha(expected_sha256):
        raise ValueError('Receipt SHA256 mismatch')
    receipt = json.loads(raw, object_pairs_hook=_unique_object,
                         parse_constant=lambda x: (_ for _ in ()).throw(ValueError('Nonfinite JSON')))
    _checked_scope(receipt)
    pins = _check_source_files(receipt)
    if cache_path is None:
        cache_path = Path(receipt['cache'])
        if not cache_path.is_absolute():
            cache_path = path.parent/cache_path
    cache_path = Path(cache_path)
    raw_cache = cache_path.read_bytes()
    if len(raw_cache) > 32*1024*1024 or sha256(raw_cache).hexdigest() != _sha(receipt['cache_sha256']):
        raise ValueError('Column cache SHA256 mismatch')
    arrays, scales = _cache_arrays(raw_cache, receipt)
    return receipt, arrays, scales, pins


def run(receipt_path, expected_sha256, *, cache_path=None, max_seconds=120, progress=None):
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('Explicit positive budget at most120 seconds required')
    started = perf_counter()
    receipt, arrays, scales, pins = load_inputs(receipt_path, expected_sha256, cache_path)
    counts = independent_source_counts()
    if receipt.get('source_counts') != {key: str(value) for key, value in counts.items()}:
        raise ValueError('Independent exact source counts disagree')
    words = []
    with ctx.workprec(PRECISION):
        for index, matrix in enumerate(arrays['combined']):
            if perf_counter()-started >= max_seconds:
                raise TimeoutError('Independent final-stage budget exhausted')
            exact = tuple(exact_scaled(x, scales['combined']) for x in matrix.flat)
            word = row_word((exact for _ in range(COLUMNS)))
            cached = exact_scaled(arrays['words'][index], scales['words'])
            assert_dominated(word, cached, 'labelled word '+str(index))
            words.append(word)
            if progress and (index+1) % 2048 == 0:
                progress(dict(checked_words=index+1, elapsed_seconds=perf_counter()-started))
        aggregate = minkowski(words, GROUPS)
        assert_dominated(aggregate, F(receipt['aggregate_upper_dyadic']), 'Minkowski aggregate')
        result = {}
        for key in ('F0', 'F1', 'F_le1'):
            bound = independent_prefix(aggregate, counts[key])
            upper = _endpoint(bound)
            reported = receipt['covered_bounds'][key]
            assert_dominated(bound, F(reported['bound_upper_dyadic']), key+' probability')
            logarithm = bound.log()/arb(2).log()
            log_upper = _endpoint(logarithm)
            reported_log = reported.get('log2_upper')
            if (type(reported_log) not in (int, float) or not isfinite(reported_log)
                    or log_upper > F(reported_log)):
                raise ArithmeticError('Independent logarithm exceeds reported upper')
            if upper >= F(1, 2**40) or reported.get('below_2_minus_40') is not True:
                raise ArithmeticError('Selected class does not verify forty bits')
            display = float(log_upper)
            if F(display) < log_upper:
                display = nextafter(display, float('inf'))
            result[key] = dict(bound_upper_dyadic=str(upper), log2_upper=display,
                               below_reported_upper=True, below_2_minus_40=True)
        aggregate_upper = str(_endpoint(aggregate))
    return dict(construction=CONSTRUCTION, scope=EXPECTED_SCOPE['scope'],
        independent_final_stage=True, earlier_column_proof_replayed=False,
        whole_code_certificate=False, precision_bits=PRECISION,
        word_algorithm='sixteen successive Arb512 row-vector products; no squaring',
        checked_labelled_words=len(words), all_cached_words_dominate_independent_uppers=True,
        source_count_method='independent exact alternating-edge cycle transfer',
        source_counts={key: str(value) for key, value in counts.items()},
        aggregate_upper_dyadic=aggregate_upper, covered_bounds=result,
        receipt_sha256=expected_sha256, cache_sha256=receipt['cache_sha256'],
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), verified_sources=pins,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--receipt-sha256', required=True)
    parser.add_argument('--cache', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=120)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen independent recheck')
    answer = run(args.receipt, args.receipt_sha256, cache_path=args.cache,
                 max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(covered_bounds=answer['covered_bounds'],
                         checked_words=answer['checked_labelled_words']), allow_nan=False), flush=True)
