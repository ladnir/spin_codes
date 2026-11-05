"""Three output marks for the same post-refresh r18/a8/F1 FLOAT screen.

This retunes the output-weight witness, not the construction or scope. Each
raw eight-family bank and every derived moment parameter is preserved.
No selected-class result is a whole-code distance certificate.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

from hd_post_refresh_regions import raw_region_bank, serialize_bank
from hd_post_refresh_compress import compress_bank, entrywise_curves
from hd_post_refresh_moment_gate import check_source_pins, SOURCE_COUNT
from hd_coupled_state_moment_gate import MomentFitError, entry_measure, polynomial_control
from hd_two_sum_gf16_mean_column import raw_holder_weights
from hd_xor_multigroup import LIMIT


MARKS = (Fraction(99, 100), Fraction(497, 500), Fraction(499, 500))


def run():
    check_source_pins()
    weights = raw_holder_weights()
    result = dict(marks=list(map(str, MARKS)), parameter_search='three prescribed output marks',
        construction='proposed-g16-post-refresh-coupled-inner',
        scope='homogeneous r18/arity8/source-Fmax-exactly1 only',
        active_groups=18, arity=8, source_count=str(SOURCE_COUNT),
        output_weight_limit=LIMIT, input_surrogate_floor=1296,
        directed=False, whole_code_certificate=False, controls=[])
    started = perf_counter()
    for mark in MARKS:
        begin = perf_counter()
        def check_time():
            if perf_counter()-begin > 30:
                raise TimeoutError('One mark exceeded thirty seconds')
        row = dict(mark=str(mark), complete=False)
        try:
            bank = raw_region_bank(float(mark), 18, full_max=1)
            check_time()
            row.update(provider_metadata=bank['metadata'],
                regional_log2_matrices=serialize_bank(bank['regions']))
            compressed = compress_bank(bank['regions'])
            curves = entrywise_curves(compressed['regions'])
            measure = entry_measure(curves)
            prefix = log2(SOURCE_COUNT)-LIMIT*log2(float(mark))-1153*log2(255)
            row.update(state_measure=measure, prefix_log2=prefix,
                compression=compressed['metadata'])
            row.update(polynomial_control(weights, measure, prefix, check_time=check_time))
            row['complete'] = True
        except (TimeoutError, MomentFitError) as error:
            row['incomplete_reason'] = str(error)
        row['elapsed_seconds'] = perf_counter()-begin
        result['controls'].append(row)
        print(json.dumps({key:row[key] for key in ('mark', 'complete', 'elapsed_seconds')}
            | {key:row[key] for key in ('log2_screen', 'incomplete_reason') if key in row}), flush=True)
    result.update(elapsed_seconds=perf_counter()-started,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen probe')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
