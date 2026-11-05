"""Fresh directed high-total-arity bounds for three g16 active-group classes.

This reuses the uniform-amplitude swapped-pair operator, not a saved outer
receipt. The unique-preimage comparison has density 256^(256*r)/255^(8*A0)
on a fixed r-group support and total constituent arity at least A0. Source
supports are unioned outside alpha; row sets and cohorts receive no extra
count factor. The low-total-arity complements are not covered here.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_shifted_swap_sparse as uniform
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS


DIRECTORY = Path(__file__).resolve().parent
ANCHORS = ((18, 247, F(691, 700)), (24, 317, F(172, 175)), (32, 463, F(171, 175)))
ALPHA = F(1, 5)
PACKET_BITS = 60
HELPER_PINS = {
    'hd_shifted_swap_sparse.py': '9b8875ec1433a2bc79dbd1e42d3b90815528a70ce0234069cfd03d61b030c58c',
    'hd_grouped_packet_screen.py': '99b0a558494c26da727f30ba5589ee3213df0617cb658b69de1c80dc0ebb3fa5',
    'hd_shifted_good_pair_screen.py': '900e8791c50be12ca3b1d82e023a004a7763891d3e786e0666efb78e8a5af045',
    'hd_shifted_pairs_count.py': 'bb1700278361343325007cf558751af2db6df19dbf6c2090502941f1a95397be',
    'hd_common_refresh_packets.py': '07c3e1c24663e34364cf8fd372d877ee39a17217112956f9f6fa3230b461e9b2',
    'hd_common_refresh_principal.py': '1ecd8f4243b81ca095f045d6c6b9c28529900bcfff1c5fa4382a2cd8ca259cb5',
    'hd_bucket_feedback.py': '66a1e97608046ae9e1fa6e37f6b6e5634de586fd03ae819278014d6a8ae16d38',
    'hd_blockwise_shared_transfer.py': '0245af670491b5aaab89cb6d9e620a2ebc490907eef4dea8dc1ee804085a260b',
    'hd_grouped_quiet.py': '4d84222bcaaeab1d134bb3ab8235d4082f39ebb9970e47efa110f52f9bf5a369',
    'hd_shifted_alias_certificate.py': 'c6363a12239a7573ee99b2365f4a101d2e51fc4e7043110e2e6dad779929b631',
    'hd_shifted_coverage.py': 'b5e7378a6019ad702ce66dc7c6046bdbb94eaf7f5e36cad528f724038f5a803e',
    'hd_sumrank_paired_certificate.py': '7277bd61663670fd52233919b539320d47b229e158811c4c06931951cdb69372',
    'hd_two_sum_gf16_sparse.py': 'e225f2fccb2f39d95dbfd7c6896c48957c3fd29898f4f574f406a090f6b95c17',
    'hd_systematic_mds.py': 'fe04df5613156b79215784e960e8a85e3cf6d4e17f8dc0c268c0eb911e66c5e2',
    'hd_blockwise_maps.py': '1f80f1f8c6f89509c69043b4542ba9598d56d6e1afe2805159c2c50848becab4',
}
EXPECTED_SCOPE = dict(
    construction=CONSTRUCTION, labels=list(LABELS), field_order=256, field_polynomial='0x11b',
    source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
    outer_groups=8192, constituent_length=16, constituent_dimension=8,
    outer_measure='unique-preimage-uniform-total-arity',
    regional_amplitude_measure='uniform-field-including-zero',
    density_per_fixed_source='256^(256*r)/255^(8*A0)',
    pair_swaps=True, independent_regional_pair_permutations=True,
    alternating_matchings=True, ordered_initial_zero_state=True, regional_pairs=128,
    coefficient_method='carry-free-Arb-truncated-bivariate', packet_mantissa_bits=PACKET_BITS,
    alpha=str(ALPHA), projective_normalizations=1,
    source_support_union_outside_alpha=True, source_count_method='profile_count_upper',
    rowset_count_factor=1, cohort_count_factor=1, extra_source_binomial_factor=False,
    all_source_supports_in_each_selected_high_arity_class=True,
    all_even_odd_source_splits=True, all_fixed_rowsets_above_threshold=True,
    uniform_operator_reused=True, saved_receipts_used=False, saved_float_estimates_used=False,
    sharper_toeplitz_measure_used=False, old_outer_certificates_transferred=False,
    directed=True, low_total_arity_complements_covered=False,
    entire_active_group_weight_classes_covered=False, continuous_weight_band_covered=False,
    whole_code_certificate=False)


def checked_dependencies():
    for name, digest in HELPER_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Reviewed dependency changed: {name}')
    if uniform.LIMIT != 734003:
        raise ValueError('Uniform operator target changed')


def require_fields(data, expected):
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError(f'Unexpected receipt field: {key}')


def log2(value):
    return value.log()/arb(2).log()


def high_prefix(r, a0, mark, alpha=ALPHA):
    if (type(r) is not int or type(a0) is not int or not 1 <= r <= 8192
            or not r <= a0 <= 16*r or not 0 < F(mark) < 1 or not 0 < F(alpha) <= 1):
        raise ValueError('Valid source size, total-arity threshold and rational witnesses required')
    return uniform.rational(F(alpha))*(2048*r-8*a0*log2(arb(255))
        -log2(arb(255))-734003*log2(uniform.rational(F(mark))))


def profile_counts(r, pairs=4096):
    if type(r) is not int or type(pairs) is not int or not 1 <= r <= 2*pairs:
        raise ValueError('Valid positive source weight and pair count required')
    return {(a, b): count for a in range(min(r//2, pairs)+1)
            for b in range(min(r//2, pairs)+1)
            if (count := uniform.profile_count_upper(pairs, r, a, b))}


def source_profile_terms(r, regions, *, pairs=4096, cycles=128, check_time=None):
    """Indexed-source upper counts and chronological (R_a R_b)^cycles.

    The profile count already unions source supports. No additional
    binomial source factor or ordering symmetry factor belongs here.
    """
    if type(cycles) is not int or cycles < 1:
        raise ValueError('Positive number of ordered matching cycles required')
    counts = profile_counts(r, pairs)
    for f in {x for pair in counts for x in pair}:
        values = regions.get((r, f))
        if (values is None or len(values) != 4
                or any(not x.is_finite() or not x >= 0 for x in values)):
            raise ValueError('Complete finite nonnegative regional matrix table required')
    terms = []
    for (a, b), count in counts.items():
        if check_time:
            check_time()
        left, right = regions[r, a], regions[r, b]
        pair = tuple(left[2*i]*right[j]+left[2*i+1]*right[2+j]
                     for i in range(2) for j in range(2))
        run = uniform.scalar_power(pair, cycles)
        value = count*(run[0]+run[1])
        if not value.is_finite() or not value > 0:
            raise ArithmeticError('Positive finite source-profile contribution required')
        terms.append((a, b, count, value))
    return terms


def anchor_metadata(r, a0, mark):
    return dict(active_groups=r, total_arity_minimum=a0, total_arity_maximum=16*r,
        mark=str(mark), alpha=str(ALPHA), density_q_exponent=256*r,
        density_nonzero_exponent=8*a0, original_indexed_source_supports=comb(8192, r),
        profiles=len(profile_counts(r)), coefficient_maximum=r,
        coefficient_marker=str(F(r, 4096)), packet_mantissa_bits=PACKET_BITS)


def receipt(precision=256, *, max_seconds_per_anchor=60, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Fresh replay precision must be256 or384')
    if (type(max_seconds_per_anchor) not in (int, float)
            or not isfinite(max_seconds_per_anchor) or not 0 < max_seconds_per_anchor <= 60):
        raise ValueError('A finite positive per-anchor budget at most60 seconds is required')
    checked_dependencies()
    started = perf_counter()
    entries = []
    with ctx.workprec(precision):
        total = arb(0)
        for r, a0, mark in ANCHORS:
            anchor_started = perf_counter()
            def check_time():
                if perf_counter()-anchor_started >= max_seconds_per_anchor:
                    raise TimeoutError('High-arity replay budget exceeded; no complete receipt produced')
            check_time()
            if progress:
                progress(dict(stage='regional-bank-start', active_groups=r, precision_bits=precision))
            regions = uniform.regional_directed(mark, ALPHA, r, PACKET_BITS)
            check_time()
            terms = source_profile_terms(r, regions, check_time=check_time)
            prefix = high_prefix(r, a0, mark)
            union = sum((term[3] for term in terms), arb(0))
            value_log = prefix+log2(union)
            rows = [dict(even_full=a, odd_full=b, indexed_source_count_upper=count,
                failure_log2_upper=str(uniform.rational_ceiling(prefix+log2(value))))
                    for a, b, count, value in terms]
            ceiling = uniform.rational_ceiling(value_log)
            total += arb(2)**value_log
            check_time()
            row = dict(anchor_metadata(r, a0, mark), profile_entries=rows,
                failure_log2_upper=str(ceiling), failure_log2_enclosure=str(value_log),
                scoped_class_reaches_40_bits=ceiling < -40,
                elapsed_seconds=perf_counter()-anchor_started)
            entries.append(row)
            if progress:
                progress(dict(stage='anchor-complete', **{k: v for k, v in row.items() if k != 'profile_entries'}))
        aggregate = uniform.rational_ceiling(log2(total))
    return dict(EXPECTED_SCOPE, precision_bits=precision, reviewed_helper_sha256=HELPER_PINS,
        entries=entries, selected_high_arity_classes_complete=True, fresh_regional_banks=3,
        failure_log2_upper=str(aggregate), selected_union_reaches_40_bits=aggregate < -40,
        compute_budget_seconds_per_anchor=max_seconds_per_anchor,
        elapsed_seconds=perf_counter()-started)


def checked_rows(data, precision):
    checked_dependencies()
    require_fields(data, dict(EXPECTED_SCOPE, precision_bits=precision,
        reviewed_helper_sha256=HELPER_PINS, selected_high_arity_classes_complete=True,
        fresh_regional_banks=3, source_sha256=sha256(Path(__file__).read_bytes()).hexdigest()))
    expected = {r: (a0, mark) for r, a0, mark in ANCHORS}
    anchors, profiles = {}, {}
    for row in data.get('entries', []):
        r = row.get('active_groups')
        if type(r) is not int or r not in expected or r in anchors:
            raise ValueError('Unexpected or repeated source weight')
        a0, mark = expected[r]
        require_fields(row, anchor_metadata(r, a0, mark))
        counts, seen = profile_counts(r), set()
        for cell in row.get('profile_entries', []):
            a, b = cell.get('even_full'), cell.get('odd_full')
            if type(a) is not int or type(b) is not int or (a, b) not in counts or (a, b) in seen:
                raise ValueError('Unexpected or repeated ordered matching profile')
            seen.add((a, b))
            require_fields(cell, dict(indexed_source_count_upper=counts[a, b]))
            if type(cell.get('failure_log2_upper')) is not str:
                raise ValueError('Exact rational profile ceiling required')
            profiles[r, a, b] = F(cell['failure_log2_upper'])
        if seen != set(counts) or type(row.get('failure_log2_upper')) is not str:
            raise ValueError('Incomplete ordered matching-profile receipt')
        anchors[r] = F(row['failure_log2_upper'])
    if set(anchors) != set(expected) or type(data.get('failure_log2_upper')) is not str:
        raise ValueError('Incomplete three-anchor receipt')
    return anchors, profiles


def rounded_union(ceilings, precision=512):
    with ctx.workprec(precision):
        value = sum((arb(2)**uniform.rational(F(x)) for x in ceilings), arb(0))
        return uniform.rational_ceiling(log2(value))


def verify(directory=DIRECTORY):
    rows, receipts, inputs = [], [], []
    for precision in (256, 384):
        path = Path(directory)/f'hd_two_sum_gf16_high_arity_receipt{precision}.json'
        raw = path.read_bytes()
        data = json.loads(raw)
        rows.append(checked_rows(data, precision))
        receipts.append(data)
        inputs.append(dict(filename=path.name, sha256=sha256(raw).hexdigest(),
            precision_bits=precision, elapsed_seconds=data['elapsed_seconds']))
    if rows[0] != rows[1] or receipts[0]['failure_log2_upper'] != receipts[1]['failure_log2_upper']:
        raise ValueError('Dual precision anchor, profile or aggregate ceilings disagree')
    anchors, profiles = rows[0]
    rounded = rounded_union(anchors.values())
    raw_ceiling = F(receipts[0]['failure_log2_upper'])
    if rounded < raw_ceiling:
        raise ArithmeticError('Rounded-anchor union lies below raw union ceiling')
    return dict(EXPECTED_SCOPE, reviewed_helper_sha256=HELPER_PINS,
        receipt_precision_bits=[256, 384], independent_sum_precision_bits=512,
        selected_high_arity_classes_complete=True, all3_anchor_ceilings_agree=True,
        all_ordered_profile_ceilings_agree=True, ordered_profile_ceilings_checked=len(profiles),
        aggregate_ceilings_agree=True, raw_aggregate_log2_upper=str(raw_ceiling),
        rounded_anchor_union_log2_upper=str(rounded),
        anchor_ceilings={str(r): str(value) for r, value in anchors.items()},
        selected_union_reaches_40_bits=rounded < -40, inputs=inputs,
        receipt_arithmetic_trust='observed-fresh-directed-replays-not-proved-by-byte-hashes')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--precision', type=int, choices=(256, 384))
    mode.add_argument('--verify', action='store_true')
    parser.add_argument('--directory', type=Path, default=DIRECTORY)
    parser.add_argument('--max-seconds-per-anchor', type=float, default=60)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen directed artifact')
    result = verify(args.directory) if args.verify else receipt(args.precision,
        max_seconds_per_anchor=args.max_seconds_per_anchor,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False), flush=True)
