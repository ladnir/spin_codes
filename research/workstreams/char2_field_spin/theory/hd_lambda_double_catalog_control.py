"""Eight prescribed exact-a6 gates from one reusable dense-form catalog.

This is selected labelled coverage, not all8008 source sets. The original
set is held fixed across each full sixteen-column word. Counts retain
total input weight, not individual input-row activity. Different groups
may choose among these eight sets in the explicitly reported subset union.
The cached count table never identifies source sets by an unproved symmetry.
"""
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_lambda_double_categories as categories
import hd_lambda_double_exact6_control as exact_control


DIRECTORY = Path(__file__).resolve().parent
CATALOG = 'hd_lambda_double_categories_a6_q256.npz'
CATALOG_SHA = '4d33d1b18dccfbfb76ae250da2fabae9c9d32e04047eaf5fa20dbcc074663692'
REFERENCE = 'hd_lambda_double_exact6_9927_float.json'
REFERENCE_SHA = 'ff18bc87635c925190869c8f12cabce6e9c07c694c70f0643d6691c82000eb26'
PINS = {
    'hd_lambda_double_categories.py': 'b4d59e12e99e112599261a7e40d3f2ccaaf9833f0cb23f73993a711c0676bbaa',
    'hd_lambda_double_exact6_control.py': '07ad009cde4608750bcce12ebe0cc04407f9b554352307736a7eafbc21e54571',
}
SOURCES = (tuple(range(6)), (0,1,2,3,4,8), (0,1,2,4,8,15), (0,1,4,5,8,9),
           (0,2,4,6,8,10), (0,3,5,6,9,12), (0,1,2,3,14,15), tuple(range(10,16)))


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen exact helper changed: '+name)
    return {**categories.dependencies(), **exact_control.dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def run(*, max_seconds=120, progress=None):
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('Use a positive wall budget at most120 seconds')
    pins, started = dependencies(), perf_counter()
    raw = (DIRECTORY/REFERENCE).read_bytes()
    if sha256(raw).hexdigest() != REFERENCE_SHA:
        raise ValueError('Exact uncompressed reference changed')
    reference = json.loads(raw)
    if reference['source_pins'] != exact_control.dependencies() or not reference['complete']:
        raise ValueError('Reference source graph or completeness differs')
    catalog = categories.load_catalog(DIRECTORY/CATALOG, CATALOG_SHA)
    if (catalog.arity, catalog.q) != (6, 256):
        raise ValueError('This gate requires the exact GF256/a6 catalog')
    cache = categories.MappedHistogramCache(catalog)
    bank, measure, prefix = exact_control.tail_backend.load_measure()
    if (prefix != reference['common_prefix_log2'] or bank['source_count'] != reference['source_count']
            or bank['mark'] != reference['mark']):
        raise ValueError('Reference and cached calculation use different common measures')
    table = exact_control.tail_backend.moment_table(measure, 18)
    controls, error = [], None
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Bounded catalog-backed gate exhausted its budget')
    try:
        for source in SOURCES:
            check_time()
            counts, diagnostic = cache.scheduled_columns(source, check_time=check_time)
            digest = sha256(json.dumps([[[w, mask, count] for (w, mask), count in sorted(column.items())]
                for column in counts], separators=(',', ':')).encode()).hexdigest()
            columns = np.stack([exact_control.coefficients(column, table) for column in counts])
            tail = float(exact_control.tail_backend.positive_deficit_tail(columns, 18, 972, check_time=check_time))
            root = tail/18-48*log2(255)
            value = prefix+18*root
            replay = source == SOURCES[0]
            if replay and (digest != reference['enumeration']['ordered_joint_counts_sha256']
                    or abs(value-reference['result']['fixed_A_repeated18_bound_log2']) > 1e-9):
                raise ArithmeticError('Compression changed raw counts or completed reference bound')
            row = dict(original_source=source, normalized_fixed_A_root_log2=root,
                fixed_A_repeated18_bound_log2=value, ordered_joint_counts_sha256=digest,
                uncompressed_reference_reproduced=replay, cache=diagnostic)
            controls.append(row)
            if progress: progress(dict(original_source=source, bound_log2=value, elapsed_seconds=perf_counter()-started))
    except TimeoutError as exception:
        error = str(exception)
    complete = len(controls) == len(SOURCES)
    union = prefix+18*float(np.logaddexp2.reduce([row['normalized_fixed_A_root_log2'] for row in controls])) if complete else None
    if (dependencies() != pins or sha256((DIRECTORY/CATALOG).read_bytes()).hexdigest() != CATALOG_SHA
            or sha256((DIRECTORY/REFERENCE).read_bytes()).hexdigest() != REFERENCE_SHA):
        raise ValueError('Frozen sources or evidence changed during gate')
    return dict(construction=categories.exact.outer.CONSTRUCTION, groups=18, arity=6,
        source_class='Fmax<=1', source_count=bank['source_count'], mark=bank['mark'],
        source_sets=SOURCES, complete=complete, error=error, controls=controls,
        eight_set_mixed_group_bound_log2=union, all8008_original_rowsets_covered=False,
        common_prefix_log2=prefix, RS_factor_whole_word='255^(-864)',
        catalog=CATALOG, catalog_sha256=CATALOG_SHA, reference=REFERENCE, reference_sha256=REFERENCE_SHA,
        source_pins=pins, original_source_fixed_across_columns=True,
        directed=False, whole_code_certificate=False, native_operations=0,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    result = run(max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=result['complete'], subset_union=result['eight_set_mixed_group_bound_log2'],
        elapsed_seconds=result['elapsed_seconds']), allow_nan=False))
