"""Exact, bounded pair-only screen for two quiet transvections.

The full search space consists of all9730 unordered pairs among140 binary
quiet planes. Retain only commuting pairs independent modulo RM(1,4).
For each retained pair, test all120 distinct input-coordinate pairs.
No symmetry quotient, support-moment table, or regional operator is used.

For packed binary parity columns a,b, the three exceptional output masks
are a XOR b, FULL XOR(a AND NOT b), and FULL XOR(b AND NOT a). These are
the ratios1, lambda/(lambda+1), and its inverse. Remaining ratios and
singletons have full output support. The result is therefore an exact
local a2 minimum for every lambda outside{0,1}; it is not a code-distance
certificate. Incomplete runs retain witnesses but do not exclude unseen
configurations. All old construction receipts remain separate.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import isfinite
from pathlib import Path
from time import perf_counter

import hd_lambda_multi_transvection as core


DIRECTORY = Path(__file__).resolve().parent
CORE_SHA256 = '42491167749463d4570e5685e7b69ff606746589dce03f4a5030312b215c7ab9'
FULL = 65535
BASE_ROWS = core.binary_rows(())
BASE_COLUMNS = tuple(sum(((row >> j) & 1) << i for i, row in enumerate(BASE_ROWS)) for j in range(16))
INPUT_PAIRS = tuple(combinations(range(16), 2))


def dependencies():
    if sha256((DIRECTORY/'hd_lambda_multi_transvection.py').read_bytes()).hexdigest() != CORE_SHA256:
        raise ValueError('Frozen multi-transvection helper changed')
    return {**core.checked_dependencies(), Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def binary_columns(supports):
    """Transpose the full B+sum vvT without any field multiplication."""
    supports = core.configuration(supports)
    columns = list(BASE_COLUMNS)
    for support in supports:
        mask = sum(1 << i for i in support)
        for i in support:
            columns[i] ^= mask
    return tuple(columns)


def exceptional_masks(column_a, column_b):
    if any(type(x) is not int or not 0 <= x <= FULL for x in (column_a, column_b)):
        raise ValueError('Two packed16-bit binary parity columns required')
    return (column_a ^ column_b,
            FULL ^ (column_a & (FULL ^ column_b)),
            FULL ^ (column_b & (FULL ^ column_a)))


def exact_pair_minimum(supports):
    """All120 source pairs, all three exceptional ratios; no early rejection."""
    columns = binary_columns(supports)
    minimum, witness = 16, None
    for a, b in INPUT_PAIRS:
        ca, cb = columns[a], columns[b]
        masks = (ca ^ cb, FULL ^ (ca & (FULL ^ cb)), FULL ^ (cb & (FULL ^ ca)))
        for index, mask in enumerate(masks):
            weight = mask.bit_count()
            if weight < minimum:
                minimum = weight
                witness = dict(source=(a, b), exceptional_ratio_index=index,
                               output_mask=mask, output_labels=tuple(i for i in range(16) if mask & (1 << i)))
    return dict(minimum_output_weight=minimum, first_worst=witness,
                input_pairs_checked=120, exceptional_ratios_per_pair=3,
                singleton_output_weight=16, collapsed_columns=minimum == 0)


def example_receipt():
    """Only the four prescribed structural examples; no family-wide inference."""
    pins = dependencies()
    rows = [dict(supports=s, construction=core.construction(s),
                 circuit=core.primitive_counts(s), pair_profile=exact_pair_minimum(s))
            for s in core.concrete_candidates()]
    if dependencies() != pins:
        raise ValueError('Sources changed')
    return dict(mode='four-prescribed-examples', rows=rows, source_sha256=pins,
                complete_family_search=False, moment_tables=0, regional_provider_runs=0,
                native_benchmarks=0, old_receipts_transferred=False, whole_code_certificate=False)


def complete_pair_screen(*, max_seconds=10):
    """Optional root-run screen of all two-transvection configurations, <=30s."""
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 30):
        raise ValueError('Positive wall budget at most30 seconds required')
    pins, started = dependencies(), perf_counter()
    planes = core.quiet_supports()
    if len(planes) != 140:
        raise ArithmeticError('The quiet-plane list must have140 elements')
    masks = tuple(sum(1 << i for i in support) for support in planes)
    considered = noncommuting = dependent = 0
    rows, survivors, histogram, stop = [], [], Counter(), None
    for i, j in combinations(range(140), 2):
        if perf_counter()-started >= max_seconds:
            stop = 'Wall budget; unseen configurations are not excluded'
            break
        considered += 1
        if (masks[i] & masks[j]).bit_count() & 1:
            noncommuting += 1
            continue
        if core.binary_rank(core.BASE_SYNDROMES+(masks[i], masks[j])) != 7:
            dependent += 1
            continue
        supports = (planes[i], planes[j])
        result = exact_pair_minimum(supports)
        histogram[result['minimum_output_weight']] += 1
        row = dict(supports=supports, **result)
        rows.append(row)
        if result['minimum_output_weight'] >= 6:
            survivors.append(dict(row, circuit=core.primitive_counts(supports),
                                  construction=core.construction(supports)))
    complete = considered == 9730
    if considered != noncommuting+dependent+len(rows):
        raise ArithmeticError('Configuration accounting failed')
    if dependencies() != pins:
        raise ValueError('Sources changed during pair screen')
    return dict(mode='all-two-transvections-pair-only', field_order=256, lambda_aes_byte=188,
        quiet_plane_count=140, expected_unordered_pairs=9730, considered_pairs=considered,
        excluded_noncommuting=noncommuting, excluded_dependent_modulo_RM1=dependent,
        independently_extending_commuting_pairs=len(rows), desired_correction_rank=7,
        complete=complete, stop_reason=stop, all_retained_configurations_test_all120_input_pairs=True,
        no_symmetry_quotient=True, exact_minimum_histogram=dict(sorted(histogram.items())),
        local_minimum6_survivors=survivors, survivor_count=len(survivors), rows=rows,
        source_sha256=pins, source_pins_unchanged=True, elapsed_seconds=perf_counter()-started,
        moment_tables=0, regional_provider_runs=0, native_benchmarks=0,
        old_receipts_transferred=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--examples', action='store_true')
    mode.add_argument('--all-pairs', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=10)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = example_receipt() if args.examples else complete_pair_screen(max_seconds=args.max_seconds)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'local_minimum6_survivors', 'source_sha256')}))
